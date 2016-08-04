#!/usr/bin/perl

use warnings;
use strict;
use 5.010;

use Digest::MD5 qw();
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
#	expectedSize: If defined, the expected file size
#	expectedMD5: If defined, the expected MD5 sum of the file, in hex format
#	cachingDir: The directory where the resource is going to be
#		downloaded, preserving the paths
#	LOG: a Log::Log4perl instance
#	cachedBoost:
sub _cachedGet($$$$$$$) {
	my($ftpServer,$remotePath,$expectedSize,$expectedMD5,$cachingDir,$LOG,$cachedBoost)=@_;
	
	my $localPath = File::Spec->catfile($cachingDir,$remotePath);
	my $filesize = undef;
	my $reason = undef;
	
	if(!$cachedBoost || ! -r $localPath) {
		my $filedate = $ftpServer->mdtm($remotePath);
		if(defined($filedate)) {
			$filesize = $ftpServer->size($remotePath);
			
			if(defined($filesize)) {
				my $localBasePath = File::Basename::basename($remotePath);
				my $localRelDir = File::Basename::dirname($remotePath);
				my $localDir = File::Spec->catdir($cachingDir,$localRelDir);
				
				$LOG->warn("Remote file $remotePath has size $filesize, but it was expected to have $expectedSize")  if(defined($expectedSize) && $expectedSize != $filesize);
				
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
					if(defined($localPath)) {
						utime($filedate,$filedate,$localPath);
					} else {
						$reason = $ftpServer->message;
						$LOG->debug("($remotePath -> $targetLocalPath) ".$reason);
					}
				}
			} else {
				$localPath = undef;
				$reason = $ftpServer->message();
			}
		} else {
			$localPath = undef;
			$reason = $ftpServer->message();
		}
	} else {
		$filesize = -s $localPath;
	}
	
	# Skipping heavy checks when there is nothing to check or when file size differs
	if(!defined($expectedSize) || $expectedSize == $filesize) {
		if(defined($localPath) && defined($expectedMD5)) {
			if(open(my $F,'<', $localPath)) {
				binmode($F);
				
				eval {
					my $md5 = Digest::MD5->new();
					
					$md5->addfile($F);
					
					my $hexdigest = $md5->hexdigest();
					$LOG->error("Wrong MD5: expected $expectedMD5 , got $hexdigest on local file $localPath (from $remotePath)")  unless($hexdigest eq $expectedMD5);
				};
				
				if($@) {
					$reason = "Error while checksumming $localPath. Reason: ".$@;
					$localPath = undef;
				}
				
				close($F);
			} else {
				$reason = "Unable to open local file $localPath for MD5 validation. Reason: ".$!;
				$localPath = undef;
			}
		}
	} else {
		$reason = "Local file $localPath size of $filesize differ from expected size $expectedSize";
		$localPath = undef;
	}
	
	return wantarray ? ($localPath,$reason) : $localPath;
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

# It returns a File::Temp instance of a temporary directory
sub createTempDir() {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my $tempDir = File::Temp->newdir('bpTempDirXXXXXX', 'DIR' => $self->{workingDir});
	
	return $tempDir;
}

# It creates a FTP or SFTP connection
sub doBPConnect($$$$) {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my($protocol,$host,$user,$pass) = @_;
	
	my $LOG = $self->{LOG};
	
	my $bpDataServer = undef;
	
	if($protocol eq 'ftp') {
		$bpDataServer = Net::FTP::AutoReconnect->new($host,Debug=>0) || $LOG->logdie("FTP connection to server $host failed: ".$@);
		$bpDataServer->login($user,$pass) || $LOG->logdie("FTP login to server $host failed: ".$bpDataServer->message());
		$bpDataServer->binary();
	} elsif($protocol eq 'sftp') {
		$LOG->logdie("Unfinished protocol $protocol. Ask the developers to finish it");
		
		$bpDataServer = Net::SFTP::Foreign->new('host' => $host,'user' => $user,'password' => $pass,'fs_encoding' => 'utf8');
		$bpDataServer->die_on_error("SSH connection to server $host failed");
	} else {
		$LOG->logdie("Unknown protocol $protocol");
	}
	
	return $bpDataServer;
}

# cachedGet parameters:
#	ftpServer: A Net::FTP or Net::FTP::AutoReconnect instance
#	remotePath: The path to the resource to be fetched from the FTP server
#	expectedSize: If defined, the expected file size
#	expectedMD5: If defined, the expected MD5 sum of the file, in hex format
#	cachedBoost:
# It returns the local path to the fetched resource
sub cachedGet($$;$$$) {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	return _cachedGet($_[0],$_[1],scalar(@_)>=3 ? $_[2]: undef,(scalar(@_)>=4 ? $_[3]: undef),$self->{workingDir},$self->{LOG},(scalar(@_)>=5 ? $_[4]: undef));
}

# listing parameters:
#	ftpServer: A Net::FTP or Net::FTP::AutoReconnect instance
#	remotePath: The path to the resource to be fetched from the FTP server
# It returns the list of files in that remote path, in absolute path format (as answered by the FTP server)
sub listing($$) {
	my $self = shift;
	
	$self->{LOG}->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my $ftpServer = shift;
	my $remotePath = shift;
	
	my $p_listing = undef;
	
	$remotePath = '/'.$remotePath  unless(substr($remotePath,0,1) eq '/');
	$p_listing = $ftpServer->ls($remotePath);
	unless(defined($p_listing)) {
		$self->{LOG}->warn("Cannot directory $remotePath listing. Reason: ". $ftpServer->message);
	}
	
	return $p_listing;
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
