#!/usr/bin/perl -w

use strict;

use diagnostics;

use Config::IniFiles;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";
use lib $FindBin::Bin."/libs";

use Log::Log4perl;

# These are needed for self-registration
use BP::DCCLoader::Parsers;
use BP::DCCLoader::Parsers::AbstractInsertionParser;
use BP::DCCLoader::Parsers::CpGInsertionParser;
use BP::DCCLoader::Parsers::DNASEBedInsertionParser;
use BP::DCCLoader::Parsers::MACSBedInsertionParser;
use BP::DCCLoader::Parsers::MethRegionsBedInsertionParser;
use BP::DCCLoader::Parsers::RNASeqGFFInsertionParser;
use BP::DCCLoader::Parsers::RNASeqStarInsertionParser;
use BP::DCCLoader::Parsers::WigglerInsertionParser;

use BP::DCCLoader::Parsers::BlueprintDCCMetadataParser;

# Logging facility initialization
my $LOG;
BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::INFO, layout => "[%d{ISO8601}]%p %m%n" } );
	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
};

# Constants and static memory structures

my %DOMAIN2EXPANAL = (
	'pdna'	=>	['cs',['pdna']],
	'rnaseq'	=>	['mrna',['exp','jcn']],
	'dnase'	=>	['chro',['rreg']],
	'meth'	=>	['wgbs',['dlat']],
);

# The main program

my $testmode = undef;
my $skipmode = undef;
my $skipmodeText;
my $testmodeText;

my $doCheck = 1;
while($doCheck) {
	if(scalar(@ARGV)>0) {
		if(index($ARGV[0],'-s')==0) {
			if(index($ARGV[0],'-ss')==0) {
				$skipmode = 2;
				$skipmodeText = 'analysis metadata and primary data';
			} else {
				$skipmode = 1;
				$skipmodeText = 'only primary data';
			}
			shift(@ARGV);
		} elsif(index($ARGV[0],'-t')==0) {
			if(index($ARGV[0],'-tt')==0) {
				$testmode = 2;
				$testmodeText = 'only validating, both data and metadata';
			} else {
				$testmode = 1;
				$testmodeText = 'only validating metadata, and skipping parsing big data';
			}
			shift(@ARGV);
		} else {
			$doCheck = undef;
		}
	} else {
		$doCheck = undef;
	}
}

$LOG->info("* [SKIPMODE] Enabled skip mode, level $skipmode ($skipmodeText)")  if($skipmode);
$LOG->info("* [TESTMODE] Enabled test mode, level $testmode ($testmodeText)")  if($testmode);

if(scalar(@ARGV)>=2) {
	STDOUT->autoflush(1);
	STDERR->autoflush(1);
	my $iniFile = shift(@ARGV);
	
	# Defined outside
	my $cachingDir = shift(@ARGV);
	my $modelDomain = shift(@ARGV);
	
	$LOG->logdie('ERROR: Unknown knowledge domain '.$modelDomain)  if(defined($modelDomain) && $modelDomain ne 'sdata' && !exists($DOMAIN2EXPANAL{$modelDomain}));

	my $metadataParser = BP::DCCLoader::Parsers::BlueprintDCCMetadataParser->new($iniFile,$cachingDir);
	
	my $publicIndexPayload = $metadataParser->getPublicIndexPayload($testmode);
	my($p_loadModels,$p_storageModels) = $metadataParser->getLoadAndStorageModels();
	
	# Common objects
	my $ini = $metadataParser->getConfigIni();
	my $workingDir = $metadataParser->getWorkingDir();
	my $model = $metadataParser->getModel();
	
	my $p_donors = $publicIndexPayload->{donors};
	my $p_specimens = $publicIndexPayload->{specimens};
	my $p_samples = $publicIndexPayload->{samples};
	my $p_studies = $publicIndexPayload->{studies};
	my $p_lab = $publicIndexPayload->{lab};
	my $p_anal = $publicIndexPayload->{anal};
	my $p_primary_anal = $publicIndexPayload->{primary_anal};
	
	foreach my $loadModelName (@{$p_loadModels}) {
		$LOG->info("Storing data using $loadModelName mapper");
		my $mapper = $p_storageModels->{$loadModelName};
		
		# Now, do we need to push the metadata there?
		if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
			if($testmode) {
				$LOG->info("\t [TESTMODE]Skipping storage of metadata model");
			} else {
				$LOG->info("\t* Storing native model");
				$mapper->storeNativeModel();
			}
		}
		
		# Several hacks in a row... Yuck!
		if(!defined($modelDomain) || $modelDomain eq 'sdata') {
			my $conceptDomain = $model->getConceptDomain('sdata');
			$LOG->info("Storing ".$conceptDomain->fullname);

			my %corrConcepts = map { $_ => BP::Loader::CorrelatableConcept->new($conceptDomain->conceptHash->{$_}) } keys(%{$conceptDomain->conceptHash});
			
			my $destination = undef;
			my $bulkData = undef;
			my $entorp = undef;
			my @bulkArray = ();
			
			# donor
			$mapper->setDestination($corrConcepts{'donor'});
			
			@bulkArray = values(%{$p_donors});
			$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
			unless($testmode) {
				$destination = $mapper->getInternalDestination();
				$bulkData = $mapper->_bulkPrepare($entorp);
				$mapper->_bulkInsert($destination,$bulkData);
			} else {
				$LOG->info("[TESTMODE] Skipping storage of donors");
			}
			
			$destination = undef;
			$mapper->freeDestination();
			@bulkArray = ();
			$bulkData = undef;
			$entorp = undef;
			
			# specimen
			$mapper->setDestination($corrConcepts{'specimen'});
			
			@bulkArray = values(%{$p_specimens});
			$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
			unless($testmode) {
				$destination = $mapper->getInternalDestination();
				$bulkData = $mapper->_bulkPrepare($entorp);
				$mapper->_bulkInsert($destination,$bulkData);
			} else {
				$LOG->info("[TESTMODE] Skipping storage of specimens");
			}
			
			$destination = undef;
			$mapper->freeDestination();
			@bulkArray = ();
			$bulkData = undef;
			$entorp = undef;
			
			# sample
			$mapper->setDestination($corrConcepts{'sample'});
			
			@bulkArray = values(%{$p_samples});
			$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
			unless($testmode) {
				$destination = $mapper->getInternalDestination();
				$bulkData = $mapper->_bulkPrepare($entorp);
				$mapper->_bulkInsert($destination,$bulkData);
			} else {
				$LOG->info("[TESTMODE] Skipping storage of samples");
			}
			
			$destination = undef;
			$mapper->freeDestination();
			@bulkArray = ();
			$bulkData = undef;
			$entorp = undef;
			
		}
		
		if(!defined($modelDomain) || $modelDomain ne 'sdata') {
			my $labConceptDomain = $model->getConceptDomain('lab');
			my $labFullname = $labConceptDomain->fullname;
			
			$LOG->info("Storing $labFullname");

			if(scalar(keys(%{$p_studies}))>0) {
				my $studyConcept = 'study';
				$LOG->info("\t* ".$labConceptDomain->conceptHash->{$studyConcept}->fullname."...");
				$mapper->setDestination(BP::Loader::CorrelatableConcept->new($labConceptDomain->conceptHash->{$studyConcept}));
				
				my @bulkArray = values(%{$p_studies});
				my $entorp = $mapper->validateAndEnactEntry(\@bulkArray);
				unless($testmode) {
					my $destination = $mapper->getInternalDestination();
					my $bulkData = $mapper->_bulkPrepare($entorp);
					$mapper->_bulkInsert($destination,$bulkData);
				} else {
					$LOG->info("[TESTMODE] Skipping storage of studies");
				}
				
				$mapper->freeDestination();
			}

			my @modelDomains = defined($modelDomain)?($DOMAIN2EXPANAL{$modelDomain}) : values(%DOMAIN2EXPANAL);
			
			foreach my $p_modelDomain (@modelDomains) {
				my($expDomain,$p_analDomains) = @{$p_modelDomain};
				
				if(exists($p_lab->{$expDomain})) {
					my $destination = undef;
					my $bulkData = undef;
					my $entorp = undef;

					$LOG->info("\t* ".$labConceptDomain->conceptHash->{$expDomain}->fullname."...");
					$mapper->setDestination(BP::Loader::CorrelatableConcept->new($labConceptDomain->conceptHash->{$expDomain}));
					$entorp = $mapper->validateAndEnactEntry($p_lab->{$expDomain});
					unless($testmode) {
						$destination = $mapper->getInternalDestination();
						$bulkData = $mapper->_bulkPrepare($entorp);
						$mapper->_bulkInsert($destination,$bulkData);
					} else {
						$LOG->info("\t[TESTMODE] Skipping storage of IHEC experiment data ".$labConceptDomain->conceptHash->{$expDomain}->fullname);
					}
					
					$destination = undef;
					$mapper->freeDestination();
					$bulkData = undef;
					$entorp = undef;
					
					next  if(defined($skipmode) && $skipmode==2);
					
					foreach my $analDomain (@{$p_analDomains}) {
						if(exists($p_anal->{$analDomain})) {
							my $conceptDomain = $model->getConceptDomain($analDomain);
							my %corrConcepts = map { $_ => BP::Loader::CorrelatableConcept->new($conceptDomain->conceptHash->{$_}) } keys(%{$conceptDomain->conceptHash});
							
							$LOG->info("Storing ".$conceptDomain->fullname);
							
							if(exists($corrConcepts{'m'})) {
								$LOG->info("\t* ".$corrConcepts{'m'}->concept->fullname."...");
								$mapper->setDestination($corrConcepts{'m'});
								$entorp = $mapper->validateAndEnactEntry($p_anal->{$analDomain});
								unless($testmode) {
									$destination = $mapper->getInternalDestination();
									$bulkData = $mapper->_bulkPrepare($entorp);
									$mapper->_bulkInsert($destination,$bulkData);
								} else {
									$LOG->info("\t[TESTMODE] Skipping storage of analysis metadata ".$corrConcepts{'m'}->concept->fullname);
								}
								
								$destination = undef;
								$mapper->freeDestination();
								$bulkData = undef;
								$entorp = undef;
								
								next  if(defined($skipmode));
								
								# And here the different bulk load
								if(exists($p_primary_anal->{$analDomain})) {
									foreach my $p_primary (@{$p_primary_anal->{$analDomain}}) {
										my($analysis_id,$conceptName,$instance,$p_remote_files) = @{$p_primary};
										
										$p_remote_files = [ $p_remote_files ]  unless(ref($p_remote_files) eq 'ARRAY');
										
										my $conceptFullName = $corrConcepts{$conceptName}->concept->fullname;
										
										$mapper->setDestination($corrConcepts{$conceptName});
										
										foreach my $p_r_file (@{$p_remote_files}) {
											my($remote_file_path,$expectedSize,$expectedMD5) = @{$p_r_file}{('r_file','expectedSize','expectedMD5')};
											
											$LOG->info("\t* ".$conceptFullName." ($remote_file_path)...");
											my($local_file,$reason) = $metadataParser->dataServerGet($remote_file_path,$expectedSize,$expectedMD5);
											
											if(defined($local_file)) {
												my $f_mode = undef;
												my @f_params = ();
												
												# Compressed file detection
												if($local_file =~ /\.gz$/) {
													$f_mode = '-|';
													push(@f_params,BP::Loader::Tools::GUNZIP,'-c',$local_file);
												} else {
													$f_mode = '<';
													push(@f_params,$local_file);
												}
												
												if(open(my $F,$f_mode,@f_params)) {
													unless($testmode && $testmode==1) {
														eval {
															$instance->insert($F,$analysis_id,$mapper);
														};
														
														if($@) {
															$LOG->logwarn("Errors while processing $remote_file_path: ".$@);
														}
													} else {
														$LOG->info("\t[TESTMODE] Skipping storage of $conceptFullName ($local_file)");
													}
													close($F);
												} else {
													$LOG->logwarn("File $local_file (fetched from $remote_file_path) not processed. Reason: ".$!);
												}
												
												# At the end, free space of the huge downloaded file
												unlink($local_file);
											} else {
												$LOG->logwarn("File $remote_file_path not processed (unable to fetch it). Reason: ".$reason);
											}
										}
										
										$mapper->freeDestination();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	$metadataParser->disconnect();
	$LOG->info("Program has finished");
} else {
	print STDERR "Usage: $0 [-t|-tt] [-s|-ss] iniFile cachingDir [",join('|','sdata',sort(keys(%DOMAIN2EXPANAL))),"]\n"
}
