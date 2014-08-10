#!/usr/bin/perl -w

use strict;

#use diagnostics;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";

use Carp;
use File::Basename;
use File::Path;
use File::Spec;
use Net::FTP;
#use Net::SFTP::Foreign 1.76;
use XML::LibXML::Reader;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Relational;

use constant DCC_LOADER_SECTION => 'dcc-loader';

sub cachedGet($$$) {
	my($bpDataServer,$remotePath,$cachingDir)=@_;
	my $filedate = $bpDataServer->mdtm($remotePath);
	my $filesize = $bpDataServer->size($remotePath);
	
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
		File::Path::make_path($localDir);
		$localPath = $bpDataServer->get($remotePath,$localPath);
		utime($filedate,$filedate,$localPath)  if(defined($localPath));
	}
	
	return $localPath;
}

if(scalar(@ARGV)>=2) {
	my $iniFile = shift(@ARGV);
	my $cachingDir = shift(@ARGV);
	my $modelDomain = shift(@ARGV);
	
	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
	
	# Check the needed parameters for the SSH connection
	my $protocol = undef;
	my $host = undef;
	my $user = undef;
	my $pass = undef;
	my $indexPath = undef;
	
	if($ini->exists(DCC_LOADER_SECTION,'protocol')) {
		$protocol = $ini->val(DCC_LOADER_SECTION,'protocol');
	} else {
		Carp::croak("Configuration file $iniFile must have 'protocol'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,'host')) {
		$host = $ini->val(DCC_LOADER_SECTION,'host');
	} else {
		Carp::croak("Configuration file $iniFile must have 'host'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,'user')) {
		$user = $ini->val(DCC_LOADER_SECTION,'user');
	} else {
		$user = 'ftp'  if($protocol eq 'ftp');
		Carp::croak("Configuration file $iniFile must have 'user'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,'pass')) {
		$pass = $ini->val(DCC_LOADER_SECTION,'pass');
	} else {
		$pass = ($user eq 'ftp')?'guest@':''  if($protocol eq 'ftp');
		Carp::croak("Configuration file $iniFile must have 'pass'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,'index-path')) {
		$indexPath = $ini->val(DCC_LOADER_SECTION,'index-path');
	} else {
		Carp::croak("Configuration file $iniFile must have 'index-path'");
	}
	
	# First, explicitly create the caching directory
	File::Path::make_path($cachingDir);
	
	print "Connecting to $host...\n";
	my $bpDataServer = undef;
	if($protocol eq 'ftp') {
		$bpDataServer = Net::FTP->new($host,Debug=>0) || Carp::croak("FTP connection to server $host failed: ".$@);
		$bpDataServer->login($user,$pass) || Carp::croak("FTP login to server $host failed: ".$bpDataServer->message());
		$bpDataServer->binary();
		
	} elsif($protocol eq 'sftp') {
		Carp::croak("Unfinished protocol $protocol. Ask the developers to finish it");
		
		$bpDataServer = Net::SFTP::Foreign->new('host' => $host,'user' => $user,'password' => $pass,'fs_encoding' => 'utf8');
		$bpDataServer->die_on_error("SSH connection to server $host failed");
	} else {
		Carp::croak("Unknown protocol $protocol");
	}
	
	my $localIndexPath = cachedGet($bpDataServer,$indexPath,$cachingDir);
	
	if(defined($localIndexPath)) {
		# Try getting a connection to 
		
		# Let's parse the model
		my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
		# Setting up the right path on relative cases
		$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

		print "Parsing model $modelFile...\n";
		my $model = undef;
		eval {
			$model = BP::Model->new($modelFile);
		};
		
		if($@) {
			Carp::croak('ERROR: Model parsing and validation failed. Reason: '.$@);
		}
		print "\tDONE!\n";
		
		
	} else {
		Carp::croak("FATAL ERROR: Unable to fetch index $indexPath from $host");
	}
	
	$bpDataServer->disconnect()  if($bpDataServer->can('disconnect'));
	$bpDataServer->quit()  if($bpDataServer->can('quit'));
	
} else {
	print STDERR "Usage: $0 iniFile cachingDir [sdata|pdna|rnaseq|dnase]\n"
}
