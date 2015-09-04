#!/usr/bin/perl

use warnings;
use strict;

use File::Basename ();
use File::Path ();
use File::Spec ();
use File::Temp ();

# Really not needed, but it enforces the architecture
use Net::FTP::AutoReconnect;
#use Net::SFTP::Foreign 1.76;

package BP::DCCLoader::WorkingDir;

use Log::Log4perl;

use constant {
	ANONYMOUS_USER	=> 'ftp',
	ANONYMOUS_PASS	=> 'guest@',
};

# _cachedGet parameters:
#	ftpServer: A Net::FTP or Net::FTP::AutoReconnect instance
#	remotePath: The resource to be fetched using the FTP server
#	cachingDir: The directory where the resource is going to be
#		downloaded, preserving the paths
#	LOG: a Log::Log4perl instance
sub _cachedGet($$$$) {
	my($ftpServer,$remotePath,$cachingDir,$LOG)=@_;
	
	my $filedate = $ftpServer->mdtm($remotePath);
	my $filesize = $ftpServer->size($remotePath);
	
	my $localPath = File::Spec->catfile($cachingDir,$remotePath);
	my $localBasePath = File::Basename::basename($remotePath);
	my $localRelDir = File::Basename::dirname($remotePath);
	my $localDir = File::Spec->catdir($cachingDir,$localRelDir);
	
	my $mirrored = undef;
	if(-f $localPath) {
		my($localsize,$localdate) = ( stat($localPath) )[7,9];
		$mirrored = $filedate == $localdate && $filesize == $localsize;
	}
	
	unless($mirrored) {
		$remotePath = '/'.$remotePath  unless(substr($remotePath,0,1) eq '/');
		File::Path::make_path($localDir);
		#print STDERR join(" -=- ",$remotePath,$cachingDir,$localPath,$localBasePath,$localRelDir,$localDir),"\n";
		my $targetLocalPath = $localPath;
		$localPath = $ftpServer->get($remotePath,$localPath);
		$LOG->debug("($remotePath -> $targetLocalPath) ".$ftpServer->message)  unless(defined($localPath));
		utime($filedate,$filedate,$localPath)  if(defined($localPath));
	}
	
	return $localPath;
}

sub new(;$) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = bless({},$class);
	
	# The optional parameter
	my $workingDir = shift;
	
	# First, explicitly create the caching directory
	if(defined($workingDir)) {
		File::Path::make_path($workingDir);
	} else {
		my $tempDir = File::Temp->newdir('bpCachedGetXXXXXX',TMPDIR => 1);
		$self->{_tempdir} = $tempDir;
		$workingDir = $tempDir->dirname;
	}
	$self->{LOG} = Log::Log4perl->get_logger(__PACKAGE__);
	$self->{workingDir} = $workingDir;
	
	return $self;
}

# cachedGet parameters:
#	ftpServer: A Net::FTP or Net::FTP::AutoReconnect instance
#	remotePath: The path to the resource to be fetched from the FTP server
# It returns the local path to the fetched resource
sub cachedGet($$) {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	return _cachedGet($_[0],$_[1],$self->{workingDir},$self->{LOG});
}

# mirror parameters:
#	resourceURI: An URI instance, pointing to the resource to be fetched
#	localPath: a relative name or path where the resource is going to be stored (optional)
#	ua: A LWP::UserAgent instance (optional)
# It returns the local path to the fetched resource
sub mirror($;$$) {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my $resourceURI = shift;
	my $localPath = shift;
	my $ua = shift;
	
	my @segments = ();
	if(defined($localPath)) {
		push(@segments,$localPath);
	} else {
		(undef,@segments) = $resourceURI->path_segments();
	}
	
	$ua = LWP::UserAgent->new()  unless(defined($ua));
	
	my $workpathLocal = File::Spec->catfile($self->{workingDir},@segments);
	my $workpathDir = File::Basename::dirname($workpathLocal);
	# Assuring the path exists
	File::Path::make_path($workpathDir);
	
	my $res = $ua->mirror($resourceURI->as_string,$workpathLocal);
	
	$self->{LOG}->logdie("FATAL ERROR: Unable to fetch ".$resourceURI->as_string.". Reason: ".$res->status_line)  unless($res->is_success || $res->is_redirect);
	
	return $workpathLocal;
}

1;
