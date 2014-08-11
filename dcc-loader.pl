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
use Tie::IxHash;
use XML::LibXML::Reader;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Relational;

use TabParser;

use constant DCC_LOADER_SECTION => 'dcc-loader';

use constant PUBLIC_INDEX => 'public.results.index';

use constant PUBLIC_INDEX_COLS => [
	# sdata_donor
	'DONOR_ID',	# donor_id
	'DONOR_SEX',	# donor_sex
	'DONOR_REGION_OF_RESIDENCE', # donor_region_of_residence
	'DONOR_ETHNICITY',	# donor_ethnicity
	
	# sdata_specimen
	# As we don't have an specimen_id, it will be derived from 'SAMPLE_ID'
	'TISSUE_TYPE',	# tissue_type
	'SAMPLE_DESC_1',	# tissue_depot
	# specimen_term is derived from TISSUE_TYPE
	# collection_method comes from IHEC metadata (COLLECTION_METHOD)
	'DONOR_AGE',	# donor_min_age_at_specimen_acquisition and donor_max_age_at_specimen_acquisition are extracted from this
	'DONOR_HEALTH_STATUS',	# donor_health_status
	'DISEASE_ONTOLOGY_URI',	# donor_disease, it contains the URI of NCI term or nothing if healthy
	'DISEASE',	# donor_disease_text
	'SPECIMEN_PROCESSING',	# specimen_processing and specimen_processing_other
	'SPECIMEN_STORAGE',	# specimen_storage and specimen_storage_other
	'BIOMATERIAL_PROVIDER',	# specimen_biomaterial_provider
	'SAMPLE_DESC_2',	# specimen_biomaterial_id ???
	# specimen_available is unknown
	# donor_id is already got from 'DONOR_ID'
	
	# sdata_sample
	'SAMPLE_ID',	# sample_id
	'SAMPLE_ONTOLOGY_URI',	# purified_cell_type, it contains the URI of CellOntology term
	# culture_conditions comes from IHEC metadata (CULTURE_CONDITIONS)
	# markers comes from IHEC metadata (MARKERS)
	# analyzed_sample_type must be set to 11 (other) to simplify
	'SAMPLE_DESC_3',	# analyzed_sample_type_other
	# analyzed_sample_interval is unknown
	# specimen_id is already derived from 'SAMPLE_ID'
	
	
];

my %SEXCV = (
	'Male'	=>	'm',
	'Female'	=>	'f',
);

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

sub parseIHECSample($$$$) {
	my($bpDataServer,$metadataPath,$sample_id,$cachingDir) = @_;
	
	my $localIHECSample = cachedGet($bpDataServer,join('/',$metadataPath,substr($sample_id,0,6),$sample_id.'.xml'),$cachingDir);
	
	if(defined($localIHECSample)) {
		my $ihec = XML::LibXML::Reader->new(location=>$localIHECSample);
		my %IHECsample = ();
		
		eval {
			while($ihec->nextElement('SAMPLE_ATTRIBUTE')>0) {
				if($ihec->nextElement('TAG')>0) {
					my $tag = $ihec->value();
					my $value = undef;
					
					$value = $ihec->value()  if($ihec->nextSiblingElement('VALUE')>0);
					
					$IHECsample{$tag} = $value;
				}
			}
		};
		
		$ihec->close();
		
		return \%IHECsample;
	} else {
		Carp::croak("Unable to fetch metadata file about sample $sample_id");
	}
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
	my $metadataPath = undef;
	
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
	
	if($ini->exists(DCC_LOADER_SECTION,'metadata-path')) {
		$metadataPath = $ini->val(DCC_LOADER_SECTION,'metadata-path');
	} else {
		Carp::croak("Configuration file $iniFile must have 'metadata-path'");
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
	
	my $localIndexPath = cachedGet($bpDataServer,$indexPath.'/'.PUBLIC_INDEX,$cachingDir);
	
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
		
		my %storageModels = ();
		
		# Setting up the loader storage model(s)
		Carp::croak('ERROR: undefined destination storage model')  unless($ini->exists($BP::Loader::Mapper::SECTION,'loaders'));
		my $loadModelNames = $ini->val($BP::Loader::Mapper::SECTION,'loaders');
		
		my @loadModels = ();
		foreach my $loadModelName (split(/,/,$loadModelNames)) {
			unless(exists($storageModels{$loadModelName})) {
				$storageModels{$loadModelName} = BP::Loader::Mapper->newInstance($loadModelName,$model,$ini);
				push(@loadModels,$loadModelName);
			}
		}
		# Now, do we need to push the metadata there?
		if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
			foreach my $mapper (@storageModels{@loadModels}) {
				$mapper->storeNativeModel();
			}
		}
		
		my %donors = ();
		tie(%donors,'Tie::IxHash');
		
		my %specimens = ();
		tie(%specimens,'Tie::IxHash');
		
		my %samples = ();
		tie(%samples,'Tie::IxHash');
		
		foreach my $mapper (@storageModels{@loadModels}) {
			# First, let's parse the public.site.index, the backbone
			if(open(my $PSI,'<',$localIndexPath)) {
				my %config = (
					TabParser::TAG_HAS_HEADER	=> 1,
					TabParser::TAG_FETCH_COLS => PUBLIC_INDEX_COLS,
					TabParser::TAG_CALLBACK => sub {
						my(
							$donor_id,
							$donor_sex,
							$donor_region_of_residence,
							$donor_ethnicity,
							
							$tissue_type,
							$tissue_depot,
							$donor_age,
							$donor_health_status,
							$donor_disease_uri,
							$donor_disease_text,
							$specimen_processing,
							$specimen_storage,
							$specimen_biomaterial_provider,
							$specimen_biomaterial_id,
							
							$sample_id,
							$purified_cell_type_uri,
							$analyzed_sample_type_other
						)=@_;
						
						unless(exists($donors{$donor_id})) {
							my $donor_region_of_residence_term = undef;
							$donor_region_of_residence_term = 'ALIAS:EAL'  if($donor_region_of_residence eq "East Anglia");
							
							my %donor = (
								'donor_id'	=>	$donor_id,
								'donor_sex'	=>	exists($SEXCV{$donor_sex})?$SEXCV{$donor_sex}:undef,
								'donor_region_of_residence'	=>	[$donor_region_of_residence_term],
								'donor_ethnicity'	=>	$donor_ethnicity,
							);
							$donors{$donor_id} = \%donor;
						}
						
						my $specimen_id = $sample_id . '_spec';
						my $p_IHECsample = undef;
						unless(exists($specimens{$specimen_id})) {
							my $donor_min_age_at_specimen_acquisition = undef;
							my $donor_max_age_at_specimen_acquisition = undef;
							if($donor_age =~ /(\d+)\s+-\s+(\d+)/) {
								$donor_min_age_at_specimen_acquisition = 'P'.$1.'Y';
								$donor_max_age_at_specimen_acquisition = 'P'.$2.'Y';
							} elsif($donor_age =~ /(\d+)\s+weeks/) {
								$donor_min_age_at_specimen_acquisition = $donor_max_age_at_specimen_acquisition = 'P'.$1.'M';
							}
							
							my $donor_disease = ($donor_disease_text eq 'None')? 'EFO_0000761': undef;
							$donor_disease = $1  if($donor_disease_uri =~ /code=([^= ]+)/);
							
							my $specimen_term = undef;
							if ($tissue_type eq "Peripheral blood"){
								$specimen_term = "UBERON:0013756";
							} elsif ($tissue_type eq "Cord blood"){
								$specimen_term = "UBERON:0012168";
							} elsif ($tissue_type eq "Tonsil"){
								$specimen_term = "UBERON:0002372";
							} elsif ($tissue_type eq "Bone marrow"){
								$specimen_term = "UBERON:0002371";
							}
							
							$p_IHECsample = parseIHECSample($bpDataServer,$metadataPath,$sample_id,$cachingDir);
							
							my %specimen = (
								'specimen_id'	=>	$specimen_id,
								'tissue_type'	=>	$tissue_type,
								'tissue_depot'	=>	$tissue_depot,
								'specimen_term'	=>	$specimen_term,
								'collection_method'	=>	exists($p_IHECsample->{COLLECTION_METHOD})?$p_IHECsample->{COLLECTION_METHOD}:undef,
								'donor_min_age_at_specimen_acquisition'	=>	$donor_min_age_at_specimen_acquisition,
								'donor_max_age_at_specimen_acquisition'	=>	$donor_max_age_at_specimen_acquisition,
								'donor_health_status'	=>	$donor_health_status,
								'donor_disease'	=>	$donor_disease,
								'donor_disease_text'	=>	$donor_disease_text,
								'specimen_processing'	=>	9,
								'specimen_processing_other'	=>	$specimen_processing,
								'specimen_storage'	=>	7,
								'specimen_storage_other'	=>	$specimen_storage,
								'specimen_biomaterial_provider'	=>	$specimen_biomaterial_provider,
								'specimen_biomaterial_id'	=>	$specimen_biomaterial_id,
								'specimen_available'	=>	undef,
								'donor_id'	=>	$donor_id,
							);
							$specimens{$specimen_id} = \%specimen;
						}
						
						unless(exists($samples{$sample_id})) {
							$p_IHECsample = parseIHECSample($bpDataServer,$metadataPath,$sample_id,$cachingDir)  unless(defined($p_IHECsample));
							
							my $purified_cell_type = undef;
							$purified_cell_type = $1  if($purified_cell_type_uri =~ /obo\/([^\/]+)/);
							my %sample = (
								'sample_id'	=>	$sample_id,
								'purified_cell_type'	=>	$purified_cell_type,
								'culture_conditions'	=>	exists($p_IHECsample->{CULTURE_CONDITIONS})?$p_IHECsample->{CULTURE_CONDITIONS}:undef,
								'markers'	=>	exists($p_IHECsample->{MARKERS})?$p_IHECsample->{MARKERS}:undef,
								'analyzed_sample_type'	=>	11,
								'analyzed_sample_type_other'	=>	$analyzed_sample_type_other,
								'analyzed_sample_interval'	=>	undef,
								'specimen_id'	=>	$specimen_id,
							);
							$samples{$sample_id} = \%sample;
						}
						
						
					},
				);
				TabParser::parseTab($PSI,%config);
				close($PSI);
			
				# Several hacks in a row... Yuck!
				if(!defined($modelDomain) || $modelDomain eq 'sdata') {
				}
				
				if(!defined($modelDomain) || $modelDomain eq 'pdna') {
				}
				
				if(!defined($modelDomain) || $modelDomain eq 'rnaseq') {
				}
				
				if(!defined($modelDomain) || $modelDomain eq 'dnase') {
				}
				
				if(!defined($modelDomain) || $modelDomain eq 'meth') {
				}
			}
		}
	} else {
		Carp::croak("FATAL ERROR: Unable to fetch index $indexPath from $host");
	}
	
	$bpDataServer->disconnect()  if($bpDataServer->can('disconnect'));
	$bpDataServer->quit()  if($bpDataServer->can('quit'));
	
} else {
	print STDERR "Usage: $0 iniFile cachingDir [sdata|pdna|rnaseq|dnase|meth]\n"
}
