#!/usr/bin/perl -w

use v5.12;
use warnings qw(all);
no warnings qw(experimental);
use strict;
use diagnostics;

use Config::IniFiles;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";
use lib $FindBin::Bin."/libs";

use Carp;
use File::Basename;
use File::Spec;
use IO::Handle;
use Log::Log4perl;
use Net::FTP::AutoReconnect;
use URI;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Tools;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Autoload::Relational;
use BP::Loader::Mapper::Autoload::Elasticsearch;
use BP::Loader::Mapper::Autoload::MongoDB;

use BP::DCCLoader::Parsers;
use BP::DCCLoader::WorkingDir;
use BP::DCCLoader::Parsers::GencodeGTFParser;

use TabParser;

my $LOG;
BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::INFO, layout => "[%d{ISO8601}]%p %m%n" } );
	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
};

use constant {
	REACTOME_BASE_TAG	=> 'reactome-base-uri',
	REACTOME_BUNDLE_TAG	=> 'reactome-bundle',
	REACTOME_INTERACTIONS_FILE	=> 'homo_sapiens.interactions.txt.gz',
	REACTOME_PATHWAYS_FILE	=> 'Ensembl2Reactome_All_Levels.txt',
	REACTOME_COMPLEXES_FILE	=> 'curated_complexes.stid.txt',
	
	APPRIS_BASE_TAG	=> 'appris-base-uri',
	APPRIS_PRINCIPAL_FILE	=> 'appris_data.principal.txt',
};

sub parseAPPRISPrincipal($$$);
sub parseReactomeComplexDesc($$$);
sub parseReactomeInteractions($$$$$);
sub parseReactomePathways($$$$);

#####
# Method bodies
#####
use constant PRINCIPAL_PREFIX	=> 'PRINCIPAL:';
use constant PRINCIPAL_LENGTH	=> length(PRINCIPAL_PREFIX);

sub parseAPPRISPrincipal($$$) {
	my($p_APPRISannot,$transcript_id,$principal_status)=@_;
	
	if(substr($principal_status,0,PRINCIPAL_LENGTH) eq PRINCIPAL_PREFIX) {
		my $dotpos = index($transcript_id,'.');
		$transcript_id = substr($transcript_id,0,$dotpos)  if($dotpos!=-1);
		
		$p_APPRISannot->{$transcript_id} = substr($principal_status,PRINCIPAL_LENGTH);
		#push(@{$p_ENShash->{$transcript_id}{'attribute'}},{'domain' => 'APPRIS_PRINCIPAL', 'value' => [substr($principal_status,PRINCIPAL_LENGTH)]});
	}
	
	1;
}

my @TRANSKEYS = ('chromosome','chromosome_start','chromosome_end');

use constant REACTOME_NS => 'Reactome';

sub parseReactomePathways($$$$) {
	my($payload,$participantId,$stablePathwayId,$pathwayName)=@_;
	
	my($p_ENShash,$p_pathways)=@{$payload};
	
	if(exists($p_ENShash->{$participantId})) {
		$p_pathways->{$stablePathwayId} = {
			'feature_cluster_id'	=> [$stablePathwayId],
			'feature_ns'	=> REACTOME_NS,
			'feature_id'	=> $stablePathwayId,
			'feature'	=> 'pathway',
			'coordinates'	=> [],
			'symbol'	=> [
				{
					'domain' => REACTOME_NS,
					'value' => [$stablePathwayId]
				},
				{
					'domain' => 'description',
					'value' => [ $pathwayName ]
				}
			]
		}  unless(exists($p_pathways->{$stablePathwayId}));
		
		my $participant = { 'feature_id' => exists($p_ENShash->{$participantId}) ? $p_ENShash->{$participantId}{'feature_id'} : $participantId };
		
		@{$participant}{@TRANSKEYS} = @{$p_ENShash->{$participantId}{'coordinates'}[0]}{@TRANSKEYS}  if(exists($p_ENShash->{$participantId}));
		
		push(@{$p_pathways->{$stablePathwayId}{'coordinates'}},$participant);
	}
	
	1;
}

sub parseReactomeComplexDesc($$$) {
	my($p_complexDesc,$pathway_id,$desc)=@_;
	
	$p_complexDesc->{$pathway_id} = $desc  if(length($pathway_id) > 0);
	
	1;
}

sub parseReactomeInteractions($$$$$) {
	my($payload,$left,$right,$pathway_id,$type)=@_;
	
	my($p_ENShash,$p_pathways,$p_foundInter,$p_complexDesc)=@{$payload};
	
	$left = ''  unless(defined($left));
	$right= ''  unless(defined($right));
	foreach my $ensId (split(/\|/,$left),split(/\|/,$right)) {
		my $colon = index($ensId,':');
		if($colon != -1) {
			$ensId = substr($ensId,$colon+1);
			unless(exists($p_foundInter->{$pathway_id}{$ensId})) {
				$p_foundInter->{$pathway_id}{$ensId} = undef;
				
				unless(exists($p_pathways->{$pathway_id})) {
					my @symbols = (
						{
							'domain' => REACTOME_NS,
							'value' => [$pathway_id]
						}
					);
					push(@symbols,{
						'domain' => 'description',
						'value' => [ $p_complexDesc->{$pathway_id} ]
					})  if(exists($p_complexDesc->{$pathway_id}));
					$p_pathways->{$pathway_id} = {
						'feature_cluster_id'	=> [$pathway_id],
						'feature_ns'	=> REACTOME_NS,
						'feature_id'	=> $pathway_id,
						'feature'	=> $type,
						'coordinates'	=> [],
						'symbol'	=> \@symbols,
					};
				}
				
				my $participant = { 'feature_id' => exists($p_ENShash->{$ensId}) ? $p_ENShash->{$ensId}{'feature_id'} : $ensId };
				
				@{$participant}{@TRANSKEYS} = @{$p_ENShash->{$ensId}{'coordinates'}[0]}{@TRANSKEYS}  if(exists($p_ENShash->{$ensId}));
				
				push(@{$p_pathways->{$pathway_id}{'coordinates'}},$participant);
			}
		}
	}
	
	1;
}

sub postprocessFeatures(\@) {
	my $p_data = shift;
	
	# Post-processing
	foreach my $p_regionData (@{$p_data}) {
		my @keywords = ();
		foreach my $symbol (@{$p_regionData->{'symbol'}}) {
			push(@keywords,@{$symbol->{'value'}});
		}
		$p_regionData->{'keyword'} = \@keywords;
	}
}

#####
# main
#####
my $testmode = undef;
if(scalar(@ARGV)>0 && $ARGV[0] eq '-t') {
	$testmode = 1;
	shift(@ARGV);
	$LOG->info("* [TESTMODE] Enabled test mode (only validating data)");
}

my $skipReactome = undef;
if(scalar(@ARGV)>0 && $ARGV[0] eq '-skipReactome') {
	$skipReactome = 1;
	shift(@ARGV);
}

if(scalar(@ARGV)>=2) {
	STDOUT->autoflush(1);
	STDERR->autoflush(1);
	my $iniFile = shift(@ARGV);
	# Defined outside
	my $cachingDir = shift(@ARGV);
	
	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
	
	my $gencode_ftp_base = undef;
	my $reactome_http_base = undef;
	
	my $gencode_gtf_file = undef;
	my $reactome_bundle_file = undef;
	
	# Check the needed parameters for the construction
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,REACTOME_BASE_TAG)) {
		$reactome_http_base = URI->new($ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,REACTOME_BASE_TAG));
	} else {
		Carp::croak("Configuration file $iniFile must have '".REACTOME_BASE_TAG."'");
	}
	
	unless($skipReactome) {
		if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,REACTOME_BUNDLE_TAG)) {
			$reactome_bundle_file = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,REACTOME_BUNDLE_TAG);
		} else {
			Carp::croak("Configuration file $iniFile must have '".REACTOME_BUNDLE_TAG."'");
		}
	}
	
	my $appris_http_base = undef;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,APPRIS_BASE_TAG)) {
		$appris_http_base = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,APPRIS_BASE_TAG);
	} else {
		Carp::croak("Configuration file $iniFile must have '".APPRIS_BASE_TAG."'");
	}
	
	Carp::croak('ERROR: undefined destination storage model')  unless($ini->exists($BP::Loader::Mapper::SECTION,'loaders'));
	
	# Zeroth, load the data model
	
	# Let's parse the model
	my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
	# Setting up the right path on relative cases
	$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

	$LOG->info("Parsing model $modelFile...");
	my $model = undef;
	eval {
		$model = BP::Model->new($modelFile);
	};
	
	if($@) {
		Carp::croak('ERROR: Model parsing and validation failed. Reason: '.$@);
	}
	$LOG->info("\tDONE!");
	
	# Now, let's patch the properies of the different remote resources, using the properties inside the model
	eval {
		$model->annotations->applyAnnotations(\($appris_http_base));
		$appris_http_base = URI->new($appris_http_base);
		
		$model->annotations->applyAnnotations(\($reactome_bundle_file))  if(defined($reactome_bundle_file));
	};
	
	if($@) {
		Carp::croak("$@ (does not exist in model $modelFile, used in $iniFile)");
	}
	
	# Setting up the loader storage model(s)
	my %storageModels = ();
	
	my $loadModelNames = $ini->val($BP::Loader::Mapper::SECTION,'loaders');
	
	my @loadModels = ();
	foreach my $loadModelName (split(/,/,$loadModelNames)) {
		unless(exists($storageModels{$loadModelName})) {
			$storageModels{$loadModelName} = BP::Loader::Mapper->newInstance($loadModelName,$model,$ini);
			push(@loadModels,$loadModelName);
		}
	}
	
	my $concept = $model->getConceptDomain('external')->conceptHash->{'features'};
	my $corrConcept = BP::Loader::CorrelatableConcept->new($concept);
	
	my @mappers = @storageModels{@loadModels};
	# Now, do we need to push the metadata there?
	if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
		foreach my $mapper (@mappers) {
			if($testmode) {
				$LOG->info("\t [TESTMODE] Skipping storage of metadata model");
			} else {
				$LOG->info("\t* Storing native model");
				$mapper->storeNativeModel();
			}
			
			$mapper->setDestination($corrConcept);
		}
	}
	
	# First, explicitly create the caching directory
	my $workingDir = BP::DCCLoader::WorkingDir->new($cachingDir);
	
	# Fetching HTTP resources
	my $reactome_bundle_local = undef;
	if(defined($reactome_bundle_file)) {
		$LOG->info("Connecting to $reactome_http_base...");
		my $reactome_bundle_uri = $reactome_http_base->clone();
		my @reactSeg = $reactome_http_base->path_segments();
		pop(@reactSeg)  if($reactSeg[$#reactSeg] eq '');
		$reactome_bundle_uri->path_segments(@reactSeg,$reactome_bundle_file);
		
		$reactome_bundle_local = $workingDir->mirror($reactome_bundle_uri);
	} else {
		$LOG->info("Skipping Reactome");
	} 
	
	$LOG->info("Connecting to $appris_http_base...");
	my $appris_principal_uri = $appris_http_base->clone();
	my @apprisSeg = $appris_http_base->path_segments();
	pop(@apprisSeg)  if($apprisSeg[$#apprisSeg] eq '');
	$appris_principal_uri->path_segments(@apprisSeg,APPRIS_PRINCIPAL_FILE);
	
	my $appris_principal_local = $workingDir->mirror($appris_principal_uri);
	
	# Now, let's get the information to enrich the transcripts labelled as principal
	$LOG->info("Parsing ".$appris_principal_local);
	my %APPRISannot;
	if(open(my $APPRIS,'<',$appris_principal_local)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> \%APPRISannot,
			TabParser::TAG_CALLBACK => \&parseAPPRISPrincipal,
			TabParser::TAG_FOLLOW	=> 1,	# We had to add it in order to avoid crap from Reactome generators
			TabParser::TAG_FETCH_COLS => [2,4],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($APPRIS,%config);
		
		close($APPRIS);
	} else {
		Carp::croak("ERROR: Unable to open APPRIS principal isoforms file ".$appris_principal_local);
	}
	
	my $batchProcessor;
	if($testmode) {
		$LOG->info("[TESTMODE] Skipping storage of remaining gene and transcript coordinates");
		$batchProcessor = sub {
			my $p_data = shift;
			
			# Post-processing
			postprocessFeatures(@{$p_data});
			
			# Storing the final genes and transcripts data
			#use JSON;
			#my $j = JSON->new->pretty->allow_blessed();
			foreach my $mapper (@mappers) {
				$mapper->validateAndEnactEntry($p_data);
				#use Devel::Size qw();
				#print 'SIZE ',Devel::Size::size($p_ENShash),' TOTAL SIZE ',Devel::Size::total_size($p_ENShash),"\n";
				
				#my $before = `ps h -o vsz $$`;
				
				#foreach my $entry (values(%{$p_ENShash})) {
				#	print $j->encode($entry),"\n";
				#}
				
				#my $after = `ps h -o vsz $$`;
				#print "DEBUG: $before $after ",$after-$before,"\n";
				#exit(1);
			}
		};
	} else {
		$batchProcessor = sub {
			my $p_data = shift;
			
			# Post-processing
			postprocessFeatures(@{$p_data});
			
			foreach my $mapper (@mappers) {
				$mapper->bulkInsert($p_data);
			}
		};
	}
	
	# The storage or the test is done in the callback
	my $p_ENShash = BP::DCCLoader::Parsers::GencodeGTFParser::getGencodeCoordinates($model,$workingDir,$ini,$testmode,\%APPRISannot,$batchProcessor);
	
	# And now, the pathways!
	if(defined($reactome_bundle_local)) {
		my $localReactomeInteractionsFile = File::Spec->catfile($cachingDir,REACTOME_INTERACTIONS_FILE);
		my $localReactomePathwaysFile = File::Spec->catfile($cachingDir,REACTOME_PATHWAYS_FILE);
		my $localReactomeComplexesFile = File::Spec->catfile($cachingDir,REACTOME_COMPLEXES_FILE);
		$LOG->info("Extracting $localReactomeInteractionsFile , $localReactomePathwaysFile and $localReactomeComplexesFile from $reactome_bundle_local");
		if((-f $localReactomeInteractionsFile && -f $localReactomePathwaysFile && -f $localReactomeComplexesFile) || system('tar','-x','-f',$reactome_bundle_local,'-C',$cachingDir,'--transform=s,^.*/,,','--wildcards','*/'.REACTOME_INTERACTIONS_FILE,'*/'.REACTOME_PATHWAYS_FILE,'*/'.REACTOME_COMPLEXES_FILE)==0) {
			$LOG->info("Parsing ".$localReactomePathwaysFile);
			if(open(my $REACTPATH,'<',$localReactomePathwaysFile)) {
				my %pathways = ();
				
				my %config = (
					TabParser::TAG_CONTEXT	=> [$p_ENShash,\%pathways],
					TabParser::TAG_CALLBACK => \&parseReactomePathways,
					TabParser::TAG_FOLLOW	=> 1,	# We had to add it in order to avoid crap from Reactome generators
					TabParser::TAG_FETCH_COLS => [0,1,3],
				);
				$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
				TabParser::parseTab($REACTPATH,%config);
				
				close($REACTPATH);
				
				my @pathwaysArr = values(%pathways);
				postprocessFeatures(@pathwaysArr);
				foreach my $mapper (@mappers) {
					unless($testmode) {
						$mapper->bulkInsert(\@pathwaysArr);
					} else {
						$LOG->info("[TESTMODE] Skipping storage of pathways mappings");
						#foreach my $entry (values(%pathways)) {
						#	print $j->encode($entry),"\n";
						#}
						my $entorp = $mapper->validateAndEnactEntry(\@pathwaysArr);
					}
				}
			} else {
				Carp::croak("ERROR: Unable to open Reactome pathways file ".$localReactomePathwaysFile);
			}
			
			$LOG->info("Parsing ".$localReactomeComplexesFile);
			my %complexDesc;
			if(open(my $COMPLEX,'<',$localReactomeComplexesFile)) {
				my %config = (
					TabParser::TAG_CONTEXT	=> \%complexDesc,
					TabParser::TAG_CALLBACK => \&parseReactomeComplexDesc,
					TabParser::TAG_FOLLOW	=> 1,	# We had to add it in order to avoid crap from Reactome generators
					TabParser::TAG_FETCH_COLS => [0,2],
				);
				$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
				TabParser::parseTab($COMPLEX,%config);
				
				close($COMPLEX);
			} else {
				Carp::croak("ERROR: Unable to open Reactome complexes file ".$localReactomeComplexesFile);
			}
			
			$LOG->info("Parsing ".$localReactomeInteractionsFile);
			if(open(my $REACT,'-|',BP::Loader::Tools::GUNZIP,'-c',$localReactomeInteractionsFile)) {
				# Pathways are also stored in features
				
				#my $reactomeConcept = $model->getConceptDomain('external')->conceptHash->{'reactome'};
				#my $reactomeCorrConcept = BP::Loader::CorrelatableConcept->new($reactomeConcept);
				#foreach my $mapper (@mappers) {
				#	$mapper->setDestination($reactomeCorrConcept);
				#}
				
				my %reactions = ();
				my %foundInter = ();
				
				my %config = (
					TabParser::TAG_COMMENT	=>	'#',
					TabParser::TAG_CONTEXT	=> [$p_ENShash,\%reactions,\%foundInter,\%complexDesc],
					TabParser::TAG_CALLBACK => \&parseReactomeInteractions,
					#TabParser::TAG_ERR_CALLBACK => \&parseReactomeInteractions,
					TabParser::TAG_NUM_COLS	=> 9,	# We had to add it in order to avoid crap from Reactome generators
					TabParser::TAG_FOLLOW	=> 1,	# We had to add it in order to avoid crap from Reactome generators
					TabParser::TAG_FETCH_COLS => [1,4,7,6],
					#TabParser::TAG_POS_FILTER => [[6 => 'reaction']],
				);
				$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
				TabParser::parseTab($REACT,%config);
				
				close($REACT);
				
				undef %foundInter;
				
				my @reactionsArr = values(%reactions);
				postprocessFeatures(@reactionsArr);
				foreach my $mapper (@mappers) {
					unless($testmode) {
						$mapper->bulkInsert(\@reactionsArr);
					} else {
						$LOG->info("[TESTMODE] Skipping storage of reactions mappings");
						#foreach my $entry (@reactionsArr) {
						#	print $j->encode($entry),"\n";
						#}
						my $entorp = $mapper->validateAndEnactEntry(\@reactionsArr);
					}
					$mapper->freeDestination();
				}
			} else {
				Carp::croak("ERROR: Unable to open Reactome interactions file ".$localReactomeInteractionsFile);
			}
		} elsif ($? == -1) {
			Carp::croak("failed to execute: $!");
		} elsif ($? & 127) {
			$LOG->error(sprintf("child died with signal %d, %s coredump\n", ($? & 127),  ($? & 128) ? 'with' : 'without'));
			exit 1;
		} else {
			$LOG->error(sprintf("child exited with value %d\n", $? >> 8));
			exit 1;
		}
	}
	$LOG->info("Program has finished");
} else {
	print STDERR "Usage: $0 [-t] [-skipReactome] iniFile cachingDir\n"
}
