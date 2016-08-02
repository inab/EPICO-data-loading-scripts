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
use BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser;

use BP::DCCLoader::Parsers::BlueprintDCCMetadataParser;
use BP::DCCLoader::Parsers::NucleosomeDetectionParser;
use BP::DCCLoader::Parsers::ChromatinSegmentationParser;

# Logging facility initialization
my $LOG;
BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::INFO, layout => "[%d{ISO8601}]%p %m%n" } );
	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
};

# Constants and static memory structures

my %SECDOMAIN2EXPANAL = (
	BP::DCCLoader::Parsers::ChromatinSegmentationParser::ANALYSIS_CONCEPT_DOMAIN_NAME()	=>	[BP::DCCLoader::Parsers::ChromatinSegmentationParser::EXPERIMENT_CONCEPT_NAME(),[BP::DCCLoader::Parsers::ChromatinSegmentationParser::ANALYSIS_CONCEPT_DOMAIN_NAME()]],
	BP::DCCLoader::Parsers::NucleosomeDetectionParser::ANALYSIS_CONCEPT_DOMAIN_NAME()	=>	['cs',[BP::DCCLoader::Parsers::NucleosomeDetectionParser::ANALYSIS_CONCEPT_DOMAIN_NAME()]],
	'cosi'	=>	['mrna',['cosi']],
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
				$skipmodeText = 'analysis metadata and secondary data';
			} else {
				$skipmode = 1;
				$skipmodeText = 'only secondary data';
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
	
	$LOG->logdie('ERROR: Unknown secondary knowledge domain '.$modelDomain)  if(defined($modelDomain) && !exists($SECDOMAIN2EXPANAL{$modelDomain}));

	my $metadataParser = BP::DCCLoader::Parsers::BlueprintDCCMetadataParser->new($iniFile,$cachingDir);
	
	# With this, we are getting the CoSI related files
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
	
	# Now, we have to fill-in here the data structures related to
	
	# * nucleosome detection
	if(!defined($modelDomain) || $modelDomain eq BP::DCCLoader::Parsers::NucleosomeDetectionParser::ANALYSIS_CONCEPT_DOMAIN_NAME()) {
		BP::DCCLoader::Parsers::NucleosomeDetectionParser::GetAnalysisMetadata($metadataParser,%{$p_anal},%{$p_primary_anal});
	}
	
	# * chromatin segmentation
	if(!defined($modelDomain) || $modelDomain eq BP::DCCLoader::Parsers::ChromatinSegmentationParser::ANALYSIS_CONCEPT_DOMAIN_NAME()) {
		BP::DCCLoader::Parsers::ChromatinSegmentationParser::GetAnalysisMetadata($metadataParser,%{$p_anal},%{$p_primary_anal});
	}
	
	# Iterating over the loadModels
	foreach my $loadModelName (@{$p_loadModels}) {
		$LOG->info("Storing secondary data using $loadModelName mapper");
		my $mapper = $p_storageModels->{$loadModelName};
		
		# We need the concept domain, only to show labels
		my $labConceptDomain = $model->getConceptDomain('lab');
		my $labFullname = $labConceptDomain->fullname;
		
		my @modelDomains = defined($modelDomain)?($SECDOMAIN2EXPANAL{$modelDomain}) : values(%SECDOMAIN2EXPANAL);
		
		# Iterating over the model domains which are going to be processed
		foreach my $p_modelDomain (@modelDomains) {
			my($expDomain,$p_analDomains) = @{$p_modelDomain};
			
			if(exists($p_lab->{$expDomain})) {
				my $destination = undef;
				my $bulkData = undef;
				my $entorp = undef;

				$LOG->info("\t* ".$labConceptDomain->conceptHash->{$expDomain}->fullname."...");
				# This task is already done in blueprint-dcc-loader.pl, so skip it!
				#$mapper->setDestination(BP::Loader::CorrelatableConcept->new($labConceptDomain->conceptHash->{$expDomain}));
				#$entorp = $mapper->validateAndEnactEntry($p_lab->{$expDomain});
				#unless($testmode) {
				#	$destination = $mapper->getInternalDestination();
				#	$bulkData = $mapper->_bulkPrepare($entorp);
				#	$mapper->_bulkInsert($destination,$bulkData);
				#} else {
				#	$LOG->info("\t[TESTMODE] Skipping storage of IHEC experiment data ".$labConceptDomain->conceptHash->{$expDomain}->fullname);
				#}
				
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
									
									my $conceptFullName = $corrConcepts{$conceptName}->concept->fullname;
									
									$LOG->info("\t* Analysis $analysis_id, ".$conceptFullName);
									
									$mapper->setDestination($corrConcepts{$conceptName});
									
									my $singleFile = ref($p_remote_files) eq 'HASH';
									$p_remote_files = [ $p_remote_files ]  if($singleFile);
									
									# Mirroring all of them
									my @localFiles = ();
									my @unlinkableFiles = ();
									my $remote_file_path;
									my $reason = undef;
									
									# Preparing several files for a single batch insertion
									foreach my $p_r_file (@{$p_remote_files}) {
										my $local_file;
										
										if(exists($p_r_file->{'r_file'})) {
											# Remote resource
											my $expectedSize;
											my $expectedMD5;
											($remote_file_path,$expectedSize,$expectedMD5) = @{$p_r_file}{('r_file','expectedSize','expectedMD5')};
											$LOG->info("\t\t- Fetching remote $remote_file_path...");
											($local_file,$reason) = $metadataParser->dataServerGet($remote_file_path,$expectedSize,$expectedMD5);
										} elsif(exists($p_r_file->{'l_file'})) {
											# Already fetched resource
											$local_file = $p_r_file->{'l_file'};
											$LOG->info("\t\t- Local $local_file...");
										} else {
											$remote_file_path = '(unset)';
											$reason = 'Design error';
										}
										
										if(defined($local_file)) {
											push(@localFiles,$local_file);
											# Recording what to remove
											push(@unlinkableFiles,$local_file)  if(exists($p_r_file->{'r_file'}));
											$reason = undef;
										} else {
											last;
										}
									}
									
									if(!defined($reason) && scalar(@localFiles) > 0) {
										my $F = undef;
										
										if($singleFile) {
											$singleFile = $localFiles[0];
											my $f_mode = undef;
											my @f_params = ();
											
											# Compressed file detection
											if($singleFile =~ /\.gz$/) {
												$f_mode = '-|';
												push(@f_params,BP::Loader::Tools::GUNZIP,'-c',$singleFile);
											} else {
												$f_mode = '<';
												push(@f_params,$singleFile);
											}
											unless(open($F,$f_mode,@f_params)) {
												$LOG->logwarn("File $singleFile (fetched from $remote_file_path) not processed. Reason: ".$!);
											}
										} else {
											$F = \@localFiles;
										}
										
										if(defined($F)) {
											unless($testmode && $testmode==1) {
												eval {
													$LOG->info("\t* Processing data for $analysis_id");
													$instance->insert($F,$analysis_id,$mapper);
												};
												
												if($@) {
													$LOG->logwarn("Errors while processing remote files for analysis $analysis_id: ".$@);
												}
											} else {
												$LOG->info("\t[TESTMODE] Skipping storage of $conceptFullName for analysis $analysis_id");
											}
											
											close($F)  if($singleFile);
										#} else {
										#	$LOG->error("What?!?!?!?");
										}
										
										# At the end, free space of the huge downloaded file
										foreach my $unlinkableFile (@unlinkableFiles) {
											$LOG->info("\t\tRemoving $unlinkableFile");
											unlink($unlinkableFile);
										}
									} elsif(defined($reason)) {
										$LOG->logwarn("File $remote_file_path not processed (unable to fetch it). Reason: ".$reason);
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
	
	
	$metadataParser->disconnect();
	$LOG->info("Program has finished");
} else {
	print STDERR "Usage: $0 [-t|-tt] [-s|-ss] iniFile cachingDir [",join('|',sort(keys(%SECDOMAIN2EXPANAL))),"]\n"
}
