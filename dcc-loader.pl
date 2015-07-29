#!/usr/bin/perl -w

use strict;

use diagnostics;

use Config::IniFiles;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";
use lib $FindBin::Bin."/libs";

use Carp;
use File::Basename;
use File::Spec;
# Give autoflush to STDOUT and STDERR
use IO::Handle;
use Net::FTP::AutoReconnect;
#use Net::SFTP::Foreign 1.76;
use Tie::IxHash;
use XML::LibXML::Reader;

use BP::Model;
use BP::Loader::Tools;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Autoload::Relational;
use BP::Loader::Mapper::Autoload::Elasticsearch;
use BP::Loader::Mapper::Autoload::MongoDB;

use BP::DCCLoader::Parsers;
use BP::DCCLoader::Parsers::CpGInsertionParser;
use BP::DCCLoader::Parsers::DNASEBedInsertionParser;
use BP::DCCLoader::Parsers::MACSBedInsertionParser;
use BP::DCCLoader::Parsers::MethRegionsBedInsertionParser;
use BP::DCCLoader::Parsers::RNASeqGFFInsertionParser;
use BP::DCCLoader::Parsers::RNASeqStarInsertionParser;
use BP::DCCLoader::Parsers::WigglerInsertionParser;

use BP::DCCLoader::WorkingDir;

use TabParser;

use constant {
	PUBLIC_INDEX_DEFAULT	=>	'public.results.index',
	DATA_INDEX_DEFAULT	=>	'data_files.index',
	EXPERIMENTS2DATASETS_DEFAULT	=>	'experiments2datasets.txt'
};

use constant PUBLIC_INDEX_COLS => [
	# sdata_donor
	'DONOR_ID',	# donor_id
	'CELL_LINE',	# alternate donor_id
	'POOL_ID',	# alternate donor_id
	'POOLED_DONOR_IDS',
	'DONOR_SEX',	# donor_sex
	'DONOR_REGION_OF_RESIDENCE', # donor_region_of_residence
	'DONOR_ETHNICITY',	# donor_ethnicity
	
	# sdata_specimen
	# As we don't have an specimen_id, it will be derived from 'SAMPLE_ID'
	'TISSUE_TYPE',	# tissue_type
	'TISSUE_DEPOT',	# tissue_depot
	# specimen_term is derived from TISSUE_TYPE
	# collection_method comes from IHEC metadata (COLLECTION_METHOD)
	'DONOR_AGE',	# donor_min_age_at_specimen_acquisition and donor_max_age_at_specimen_acquisition are extracted from this
	'DONOR_HEALTH_STATUS',	# donor_health_status
	'DISEASE_ONTOLOGY_URI',	# donor_disease, it contains the URI of NCI term or nothing if healthy
	'DISEASE',	# donor_disease_text
	'TREATMENT',	# The treatment being received by donor on specimen adcquisition
	
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
	
	# lab_*
	'EXPERIMENT_ID',	# experiment_id
	# analyzed_sample_id is already got from 'SAMPLE_ID'
	# experiment_type is got from IHEC metadata (EXPERIMENT_TYPE)
	# library_strategy is got from IHEC metadata (LIBRARY_STRATEGY)
	'CENTER_NAME',	# it helps to get the experimental_group_id
	'LIBRARY_STRATEGY',	# it helps to get the kind of experiment
	'INSTRUMENT_MODEL',	# platform is derived from the contents of this column
	
	# *_m
	'FILE_TYPE',	# We decide the kind of analysis from this, and the postfix to build the analysis_id
	'FILE',	# In some cases this one defines a variant in the analysis_id
	'NSC',	# NSC
	'RSC',	# RSC
	#'WITHDRAWN',	# Has been withdrawn this analysis?
];

my %SEXCV = (
	'male'	=>	'm',
	'female'	=>	'f',
	'mixed'	=>	'b',
);

my %GROUPCV = (
	'CNAG'	=>	['11'],
	'NCMLS'	=>	['1'],
	'NCMLS_CU'	=>	['1','3b'],
	'MPIMG'	=>	['12d'],
	'IDIBAPS'	=>	['10'],
);

my %EXPERIMENTCV = (
	'Bisulfite-Seq'	=>	'wgbs',
	'ChIP-Seq'	=>	'cs',
	'DNase-Hypersensitivity'	=>	'chro',
	'RNA-Seq'	=>	'mrna',
);

my %DOMAIN2EXPANAL = (
	'pdna'	=>	['cs',['pdna']],
	'rnaseq'	=>	['mrna',['exp','jcn']],
	'dnase'	=>	['chro',['rreg']],
	'meth'	=>	['wgbs',['dlat']],
);

my %INSTRUMENT2PLATFORM = (
	'Illumina HiSeq 2000'	=>	60,
	'NextSeq 500'	=>	100,
);

# Constants to access the data on each domain in primary_anal
use constant {
	P_ANALYSIS_ID	=>	0,
	P_CONCEPTNAME	=>	1,
	P_METHOD	=>	2,
	P_FILE	=>	3,
};

#####
# Method prototypes
#####
sub parseIHECsample($$$$);
sub parseIHECexperiment($$$$);

#####
# Method callbacks
#####

sub experiments_to_datasets_callback(\%$$$) {
	my($p_exp2EGA,$experiment_id, $title, $dataset_id) = @_;

	$p_exp2EGA->{$experiment_id} = $dataset_id;
}

sub data_files_callback {
	my(
		$payload,
		
		$donor_id,
		$cell_line,
		$pool_id,
		$pooled_donor_ids,
		$donor_sex,
		$donor_region_of_residence,
		$donor_ethnicity,
		
		$tissue_type,
		$tissue_depot,
		$donor_age,
		$donor_health_status,
		$donor_disease_uri,
		$donor_disease_text,
		$donor_treatment,
		
		$specimen_processing,
		$specimen_storage,
		$specimen_biomaterial_provider,
		$specimen_biomaterial_id,
		
		$sample_id,
		$purified_cell_type_uri,
		$analyzed_sample_type_other,
		
		$experiment_id,
		$center_name,
		$library_strategy,
		$instrument_model,
		
		$file_type,
		$remote_file_path,
		$NSC,
		$RSC,
	)=@_;
	
	
	if(exists($payload->{experiments}{$experiment_id})) {
		# And this is the analysis metadata
		if(exists($payload->{FILETYPE2ANAL}{$file_type})) {
			my $ftype = $payload->{FILETYPE2ANAL}{$file_type};
			my $analDomain = $ftype->[BP::DCCLoader::Parsers::F_DOMAIN];
			
			# Analysis id building
			my $an_postfix = undef;
			
			if(defined($ftype->[BP::DCCLoader::Parsers::F_PATTERN_POSTFIX])) {
				foreach my $p_pat_post (@{$ftype->[BP::DCCLoader::Parsers::F_PATTERN_POSTFIX]}) {
					my($pattern,$postfix)=@{$p_pat_post};
					
					if(index($remote_file_path,$pattern)!=-1) {
						$an_postfix = $postfix;
						last;
					}
				}
			}
			
			$an_postfix = $ftype->[BP::DCCLoader::Parsers::F_POSTFIX]  unless(defined($an_postfix));
			
			# No postfix, no processing!!!!
			if(defined($an_postfix)) {
				my $analysis_id = $experiment_id.'_'.$an_postfix;
				
				unless(exists($payload->{reg_analysis}{$analysis_id})) {
					my $f_metadata = $ftype->[BP::DCCLoader::Parsers::F_METADATA];
					$f_metadata = {}  unless(defined($f_metadata));
					
					my %analysis = (
						'analysis_id'	=>	$analysis_id,
						'experiment_id'	=>	$experiment_id,
						'analysis_group_id'	=>	$ftype->[BP::DCCLoader::Parsers::F_ANALYSIS_GROUP_ID],
						'data_status'	=>	defined($ftype->[BP::DCCLoader::Parsers::F_PRIMARY])?2:0,
						'assembly_version'	=>	$f_metadata->{assembly_version},
						'ensembl_version'	=>	$payload->{ensembl_version},
						'gencode_version'	=>	$payload->{gencode_version},
					);
					$analysis{NSC} = $NSC  if($NSC ne '-');
					$analysis{RSC} = $RSC  if($RSC ne '-');
					@analysis{keys(%{$f_metadata})} = values(%{$f_metadata});
					if(defined($ftype->[BP::DCCLoader::Parsers::F_PARENT_POSTFIX])) {
						$analysis{'base_analysis_id'} = $analysis_id.'_'.$ftype->[BP::DCCLoader::Parsers::F_PARENT_POSTFIX];
					}
					
					# Last, register it!
					$payload->{anal}{$analDomain} = []  unless(exists($payload->{anal}{$analDomain}));
					push(@{$payload->{anal}{$analDomain}},\%analysis);
					$payload->{reg_analysis}{$analysis_id} = undef;
				}
				
				# Preparing the field
				if(defined($ftype->[BP::DCCLoader::Parsers::F_PRIMARY])) {
					$payload->{primary_anal}{$analDomain} = []  unless(exists($payload->{primary_anal}{$analDomain}));
					
					push(@{$payload->{primary_anal}{$analDomain}},[$analysis_id,$ftype->[BP::DCCLoader::Parsers::F_PRIMARY],$ftype->[BP::DCCLoader::Parsers::F_PARSER],$remote_file_path]);
				}
			}
		}
	}
}

sub public_results_callback {
	my(
		$payload,
		
		$donor_id,
		$cell_line,
		$pool_id,
		$pooled_donor_ids,
		$donor_sex,
		$donor_region_of_residence,
		$donor_ethnicity,
		
		$tissue_type,
		$tissue_depot,
		$donor_age,
		$donor_health_status,
		$donor_disease_uri,
		$donor_disease_text,
		$donor_treatment,
		
		$specimen_processing,
		$specimen_storage,
		$specimen_biomaterial_provider,
		$specimen_biomaterial_id,
		
		$sample_id,
		$purified_cell_type_uri,
		$analyzed_sample_type_other,
		
		$experiment_id,
		$center_name,
		$library_strategy,
		$instrument_model,
		
		$file_type,
		$remote_file_path,
		$NSC,
		$RSC,
	)=@_;
	
	my $donor_kind;
	
	if($donor_id eq '-') {
		if($cell_line ne '-') {
			$donor_id = $cell_line;
			$donor_kind = 'c';
		} elsif($pool_id ne '-') {
			$donor_id = $pool_id;
			$donor_kind = 'p';
		} else {
			Carp::croak("Unable to identify the kind of donor");
		}
	} else {
		$donor_kind = 'd';
	}
	
	unless(exists($payload->{donors}{$donor_id})) {
		$donor_ethnicity = undef  if($donor_ethnicity eq 'NA' || $donor_ethnicity eq '-');
		
		my $donor_region_of_residence_term = '001';	# World
		$donor_region_of_residence_term = 'ALIAS:EAL'  if($donor_region_of_residence eq "East Anglia");
		
		# To avoid case matching problems
		$donor_sex = lc($donor_sex);
		
		my %donor = (
			'donor_id'	=>	$donor_id,
			'donor_kind'	=>	$donor_kind,
			'donor_sex'	=>	exists($SEXCV{$donor_sex})?$SEXCV{$donor_sex}:undef,
			'donor_region_of_residence'	=>	[$donor_region_of_residence_term],
			'donor_ethnicity'	=>	$donor_ethnicity,
			'notes'	=> ($cell_line ne '-')?'Cell line':undef,
		);
		# For the case of pooled donor ids
		$donor{'pooled_donor_ids'} = [split(/,/,$pooled_donor_ids)]  if($donor_kind eq 'p' && $pooled_donor_ids ne '-');
		
		$payload->{donors}{$donor_id} = \%donor;
	}
	
	my $specimen_id = $sample_id . '_spec';
	my $p_IHECsample = undef;
	unless(exists($payload->{specimens}{$specimen_id})) {
		$donor_health_status = undef  if($donor_health_status eq 'NA' || $donor_health_status eq '-');
		$specimen_processing = undef  if($specimen_processing eq 'NA' || $specimen_processing eq '-');
		$specimen_storage = undef  if($specimen_storage eq 'NA' || $specimen_storage eq '-');
		
		my $donor_min_age_at_specimen_acquisition = undef;
		my $donor_max_age_at_specimen_acquisition = undef;
		if($donor_age =~ /(\d+)\s*-\s*(\d+)/) {
			$donor_min_age_at_specimen_acquisition = 'P'.$1.'Y';
			$donor_max_age_at_specimen_acquisition = 'P'.$2.'Y';
		} elsif($donor_age =~ /(\d+)\s+weeks/) {
			$donor_min_age_at_specimen_acquisition = $donor_max_age_at_specimen_acquisition = 'P'.$1.'M';
		} elsif($donor_age =~ /(\d+)\s+days/) {
			$donor_min_age_at_specimen_acquisition = $donor_max_age_at_specimen_acquisition = 'P'.$1.'D';
		} else {
			$donor_min_age_at_specimen_acquisition = 'P0Y';
		}
		
		my $donor_disease = ($donor_disease_text eq 'None')? 'http://www.ebi.ac.uk/efo/EFO_0000761': undef;
		
		my $specimen_term = undef;
		
		my @purified_term_uris = split(/;/,$purified_cell_type_uri);
		
		foreach my $term_uri (@purified_term_uris) {
			if($term_uri =~ /obo\/(?:((?:UBERON)|(?:CLO))_([^\/]+))/ || $term_uri =~ /efo\/(EFO)_([^\/]+)/) {
				$specimen_term = $term_uri;
				last;
			}
		}
		
		unless(defined($specimen_term)) {
			if ($tissue_type eq "Peripheral blood" || $tissue_type eq "Venous blood"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0013756";
			} elsif($tissue_type eq "Cord blood"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0012168";
			} elsif($tissue_type eq "Tonsil"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0002372";
			} elsif($tissue_type eq "Bone marrow"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0002371";
			} elsif($tissue_type eq "Thymus") {
				# This is more specific, but as I'm unsure....
				# $specimen_term = "http://purl.obolibrary.org/obo/UBERON_0003483";
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0002370";
			}
		}
		
		# Last resort, look at the cache
		if(defined($specimen_term)) {
			$payload->{cellSpecimenTerm}{$cell_line} = $specimen_term  if($cell_line ne '-' && !exists($payload->{cellSpecimenTerm}{$cell_line}));
		} elsif($cell_line ne '-' && exists($payload->{cellSpecimenTerm}{$cell_line})) {
			$specimen_term = $payload->{cellSpecimenTerm}{$cell_line};
		}
		
		Carp::croak("Undefined specimen term for $specimen_id!!!!")  unless(defined($specimen_term));
		
		$p_IHECsample = parseIHECsample($payload->{bpDataServer},$payload->{metadataPath},$sample_id,$payload->{workingDir});
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
		$specimen{'donor_treatment'} = $donor_treatment  if($donor_treatment ne '-' && $donor_treatment ne 'NA');
		$payload->{specimens}{$specimen_id} = \%specimen;
	}
	
	unless(exists($payload->{samples}{$sample_id})) {
		$p_IHECsample = parseIHECsample($payload->{bpDataServer},$payload->{metadataPath},$sample_id,$payload->{workingDir})  unless(defined($p_IHECsample));
		
		my $purified_cell_type = undef;

		my @purified_term_uris = split(/;/,$purified_cell_type_uri);
		
		foreach my $term_uri (@purified_term_uris) {
			if($term_uri =~ /obo\/(?:((?:CLO)|(?:CL))_([^\/]+))/ || $term_uri =~ /efo\/(EFO)_([^\/]+)/) {
				$purified_cell_type = $term_uri;
				last;
			}
		}
		
		# Last resort, look at the cache
		if(defined($purified_cell_type)) {
			$payload->{cellPurifiedTerm}{$cell_line} = $purified_cell_type  if($cell_line ne '-' && !exists($payload->{cellPurifiedTerm}{$cell_line}));
		} elsif($cell_line ne '-' && exists($payload->{cellPurifiedTerm}{$cell_line})) {
			$purified_cell_type = $payload->{cellPurifiedTerm}{$cell_line};
		} else {
			# When we don't know, it is the most general term, hematopoietic cell
			$purified_cell_type = "http://purl.obolibrary.org/obo/CL_0000988";
		}
		
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
		$payload->{samples}{$sample_id} = \%sample;
	}
	
	if(exists($EXPERIMENTCV{$library_strategy})) {
		# This is the experimental metadata
		unless(exists($payload->{experiments}{$experiment_id})) {
			my $labexp = $EXPERIMENTCV{$library_strategy};
			
			my($p_IHECexperiment,$ihec_library_strategy,$ihec_instrument_model) = parseIHECexperiment($payload->{bpDataServer},$payload->{metadataPath},$experiment_id,$payload->{workingDir});
			
			my %features = map { $_ => { 'feature' => $_ , 'value' => $p_IHECexperiment->{$_}[0], 'units' => $p_IHECexperiment->{$_}[1] } } keys(%{$p_IHECexperiment});
			
			# The common attributes
			my %experiment = (
				'experiment_id'	=>	$experiment_id,
				'analyzed_sample_id'	=>	$sample_id,
				'experiment_type'	=>	exists($p_IHECexperiment->{EXPERIMENT_TYPE})?$p_IHECexperiment->{EXPERIMENT_TYPE}[0]:'',
				'library_strategy'	=>	defined($ihec_library_strategy)?$ihec_library_strategy:$library_strategy,
				'experimental_group_id'	=>	exists($GROUPCV{$center_name})?$GROUPCV{$center_name}:[$center_name],
				'features'	=>	\%features,
				'raw_data_repository'	=>	1,
				'raw_data_accession'	=>	{
									'accession'	=>	exists($payload->{exp2EGA}{$experiment_id})?$payload->{exp2EGA}{$experiment_id}:'',
									'url'	=>	exists($payload->{exp2EGA}{$experiment_id})?('https://www.ebi.ac.uk/ega/datasets/'.$payload->{exp2EGA}{$experiment_id}):'',
								},
				'platform'	=>	exists($INSTRUMENT2PLATFORM{$instrument_model})?$INSTRUMENT2PLATFORM{$instrument_model}:-1,
				'platform_model'	=>	$ihec_instrument_model,
				'seq_coverage'	=>	undef,
				'extraction_protocol'	=>	exists($p_IHECexperiment->{EXTRACTION_PROTOCOL})?$p_IHECexperiment->{EXTRACTION_PROTOCOL}[0]:undef,
			);
			
			# Last, register it!
			$payload->{lab}{$labexp} = []  unless(exists($payload->{lab}{$labexp}));
			push(@{$payload->{lab}{$labexp}},\%experiment);
			$payload->{experiments}{$experiment_id} = undef;
		}
		
		&data_files_callback(@_);
		
	} else {
		Carp::carp("Unknown type of experiment: ".$library_strategy);
	}
}

#####
# Method bodies
#####
sub parseIHECsample($$$$) {
	my($bpDataServer,$metadataPath,$sample_id,$workingDir) = @_;
	
	print "\t* Parsing IHEC sample $sample_id...\n";
	
	my $localIHECsample = $workingDir->cachedGet($bpDataServer,join('/',$metadataPath,'samples',substr($sample_id,0,6),$sample_id.'.xml'));
	
	my %IHECsample = ();
	if(defined($localIHECsample)) {
		my $ihec = XML::LibXML::Reader->new(location=>$localIHECsample);
		
		eval {
			if($ihec->nextElement('SAMPLE_ATTRIBUTES')) {
				while($ihec->nextElement('SAMPLE_ATTRIBUTE')>0) {
					if($ihec->nextElement('TAG')>0) {
						my $tag = $ihec->readInnerXml();
						my $value = undef;
						
						$value = $ihec->readInnerXml()  if($ihec->nextSiblingElement('VALUE')>0);
						
						$IHECsample{$tag} = $value;
					}
				}
			}
		};
		
		$ihec->close();
	} else {
		Carp::carp("Unable to fetch metadata file about sample $sample_id");
	}
	
	return \%IHECsample;
}

sub parseIHECexperiment($$$$) {
	my($bpDataServer,$metadataPath,$experiment_id,$workingDir) = @_;
	
	print "\t* Parsing IHEC experiment $experiment_id...\n";
	
	my $localIHECexperiment = $workingDir->cachedGet($bpDataServer,join('/',$metadataPath,'experiments',substr($experiment_id,0,6),$experiment_id.'.xml'));
	
	my %IHECexperiment = ();
	my $library_strategy = undef;
	my $instrument_model = undef;
	if(defined($localIHECexperiment)) {
		my $pat = XML::LibXML::Pattern->new('//LIBRARY_STRATEGY | //PLATFORM//INSTRUMENT_MODEL | //EXPERIMENT_ATTRIBUTES');
		my $patnext = XML::LibXML::Pattern->new('//PLATFORM//INSTRUMENT_MODEL | //EXPERIMENT_ATTRIBUTES');
		my $ihec = XML::LibXML::Reader->new(location=>$localIHECexperiment);
		
		eval {
			if($ihec->nextPatternMatch($pat)) {
				if($ihec->localName() eq 'LIBRARY_STRATEGY') {
					$library_strategy = $ihec->readInnerXml();
				} elsif($ihec->localName() ne 'EXPERIMENT_ATTRIBUTES') {
					$instrument_model = $ihec->readInnerXml()  if($ihec->localName() eq 'INSTRUMENT_MODEL' || ($ihec->nextPatternMatch($pat) && $ihec->localName() eq 'INSTRUMENT_MODEL'));
				}
				
				if($ihec->localName() eq 'EXPERIMENT_ATTRIBUTES' || $ihec->nextElement('EXPERIMENT_ATTRIBUTES')) {
					while($ihec->nextElement('EXPERIMENT_ATTRIBUTE')>0) {
						if($ihec->nextElement('TAG')>0) {
							my $tag = $ihec->readInnerXml();
							my $value = undef;
							my $units = undef;
							
							$value = $ihec->readInnerXml()  if($ihec->nextSiblingElement('VALUE')>0);
							$units = $ihec->readInnerXml()  if($ihec->nextSiblingElement('UNITS')>0);
							
							$IHECexperiment{$tag} = [$value,$units];
						}
					}
				}
			}
		};
		
		$ihec->close();
	} else {
		Carp::carp("Unable to fetch metadata file about experiment $experiment_id");
	}
	
	return (\%IHECexperiment,$library_strategy,$instrument_model);
}

my $testmode = undef;
if(scalar(@ARGV)>0 && $ARGV[0] eq '-t') {
	$testmode = 1;
	shift(@ARGV);
	print "* [TESTMODE] Enabled test mode (only validating data)\n";
}

if(scalar(@ARGV)>=2) {
	STDOUT->autoflush(1);
	STDERR->autoflush(1);
	my $iniFile = shift(@ARGV);
	
	# Defined outside
	my $cachingDir = shift(@ARGV);
	my $modelDomain = shift(@ARGV);
	
	Carp::croak('ERROR: Unknown knowledge domain '.$modelDomain)  if(defined($modelDomain) && $modelDomain ne 'sdata' && !exists($DOMAIN2EXPANAL{$modelDomain}));
	




	######
	# These were global variables, but as we need them in the different callbacks
	######
	my %donors = ();
	tie(%donors,'Tie::IxHash');
	
	my %specimens = ();
	tie(%specimens,'Tie::IxHash');
	
	my %samples = ();
	tie(%samples,'Tie::IxHash');
	
	my %experiments = ();
	
	# Laboratory experiments and analysis metadata
	my %lab = ();
	my %anal = ();
	
	my $p_FILETYPE2ANAL;
	
	my %reg_analysis = ();
	
	my %primary_anal = ();
	
	my %expfiles = ();
	
	# Correspondence between experiment ids and EGA ids
	my %exp2EGA = ();
	
	# Cache for cell lines
	my %cellSpecimenTerm = ();
	my %cellPurifiedTerm = ();
	
	# Parameters needed
	# Either a Net::FTP or Net::FTP::AutoReconnect instance
	my $bpDataServer = undef;
	my $metadataPath = undef;
	# A BP::DCCLoader::WorkingDir instance
	my $workingDir = undef;
	
	my $ensembl_version = undef;
	my $gencode_version = undef;











	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
	
	# Check the needed parameters for the SSH connection
	my $protocol = undef;
	my $host = undef;
	my $user = undef;
	my $pass = undef;
	my $indexPath = undef;
	# Defined outside
	$metadataPath = undef;
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'protocol')) {
		$protocol = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'protocol');
	} else {
		Carp::croak("Configuration file $iniFile must have 'protocol'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'host')) {
		$host = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'host');
	} else {
		Carp::croak("Configuration file $iniFile must have 'host'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'user')) {
		$user = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'user');
	} else {
		$user = 'ftp'  if($protocol eq 'ftp');
		Carp::croak("Configuration file $iniFile must have 'user'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'pass')) {
		$pass = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'pass');
	} else {
		$pass = ($user eq 'ftp')?'guest@':''  if($protocol eq 'ftp');
		Carp::croak("Configuration file $iniFile must have 'pass'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'index-path')) {
		$indexPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'index-path');
	} else {
		Carp::croak("Configuration file $iniFile must have 'index-path'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-path')) {
		$metadataPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-path');
	} else {
		Carp::croak("Configuration file $iniFile must have 'metadata-path'");
	}
	
	my $publicIndex;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'public-index-file')) {
		$publicIndex = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'public-index-file');
	} else {
		$publicIndex = PUBLIC_INDEX_DEFAULT;
	}
	
	my $dataIndex;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'data-index-file')) {
		$dataIndex = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'data-index-file');
	} else {
		$dataIndex = DATA_INDEX_DEFAULT;
	}
	
	my $exp2datasets;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'exp2datasets-file')) {
		$exp2datasets = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'exp2datasets-file');
	} else {
		$exp2datasets = EXPERIMENTS2DATASETS_DEFAULT;
	}
	
	
	# First, explicitly create the caching directory
	$workingDir = BP::DCCLoader::WorkingDir->new($cachingDir);
	
	print "Connecting to $host...\n";
	# Defined outside
	$bpDataServer = undef;
	if($protocol eq 'ftp') {
		$bpDataServer = Net::FTP::AutoReconnect->new($host,Debug=>0) || Carp::croak("FTP connection to server $host failed: ".$@);
		$bpDataServer->login($user,$pass) || Carp::croak("FTP login to server $host failed: ".$bpDataServer->message());
		$bpDataServer->binary();
		
	} elsif($protocol eq 'sftp') {
		Carp::croak("Unfinished protocol $protocol. Ask the developers to finish it");
		
		$bpDataServer = Net::SFTP::Foreign->new('host' => $host,'user' => $user,'password' => $pass,'fs_encoding' => 'utf8');
		$bpDataServer->die_on_error("SSH connection to server $host failed");
	} else {
		Carp::croak("Unknown protocol $protocol");
	}
	
	my $localIndexPath = $workingDir->cachedGet($bpDataServer,$indexPath.'/'.$publicIndex);
	my $localDataFilesIndexPath = $workingDir->cachedGet($bpDataServer,$indexPath.'/'.$dataIndex);
	my $localExp2Datasets = $workingDir->cachedGet($bpDataServer,$indexPath.'/'.$exp2datasets);
	
	if(defined($localIndexPath)) {
		Carp::carp("WARNING: Unable to fetch experiments to datasets correspondence from $indexPath (host $host)")  unless(defined($localExp2Datasets));

		# Try getting a connection to 
		
		# Let's parse the model
		my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
		# Setting up the right path on relative cases
		$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

		print "Parsing model $modelFile...\n";
		my $model = undef;
		eval {
			$model = BP::Model->new($modelFile);
			$ensembl_version = exists($model->annotations->hash->{EnsemblVer})?$model->annotations->hash->{EnsemblVer}:'';
			$gencode_version = exists($model->annotations->hash->{GENCODEVer})?$model->annotations->hash->{GENCODEVer}:'';
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
		
		# Initializations of parsable file types, needed by next parsing tasks
		$p_FILETYPE2ANAL = BP::DCCLoader::Parsers->getParsableFiletypes({'ini'=>$ini,'model'=>$model,'workingDir'=>$workingDir});

		# First, these correspondences experiment <=> EGA needed by next parse
		if(defined($localExp2Datasets)) {
			print "Parsing ",$exp2datasets,"...\n";
			if(open(my $E2D,'<:encoding(UTF-8)',$localExp2Datasets)) {
				my %e2dConfig = (
					TabParser::TAG_HAS_HEADER	=> 1,
					TabParser::TAG_CONTEXT	=> \%exp2EGA,
					TabParser::TAG_CALLBACK => \&experiments_to_datasets_callback,
				);
				$e2dConfig{TabParser::TAG_VERBOSE} = 1  if($testmode);
				TabParser::parseTab($E2D,%e2dConfig);
				close($E2D);
			} else {
				Carp::croak("Unable to parse $localExp2Datasets, needed to get the EGA dataset identifiers");
			}
		} else {
			Carp::carp("WARNING: $exp2datasets unavailable. raw_data_accession subfields will be empty");
		}
		
		my $publicIndexPayload = {
			'donors'	=>	\%donors,
			'specimens'	=>	\%specimens,
			'samples'	=>	\%samples,
			'experiments'	=>	\%experiments,
			'lab'	=>	\%lab,
			'anal'	=>	\%anal,
			'FILETYPE2ANAL'	=>	$p_FILETYPE2ANAL,
			'reg_analysis'	=>	\%reg_analysis,
			'primary_anal'	=>	\%primary_anal,
			'expfiles'	=>	\%expfiles,
			'exp2EGA'	=>	\%exp2EGA,
			'cellSpecimenTerm'	=>	\%cellSpecimenTerm,
			'cellPurifiedTerm'	=>	\%cellPurifiedTerm,
			'workingDir'	=>	$workingDir,
			'bpDataServer'	=>	$bpDataServer,
			'metadataPath'	=>	$metadataPath,
			'ensembl_version'	=>	$ensembl_version,
			'gencode_version'	=>	$gencode_version,
		};
		
		print "Parsing ",$publicIndex,"...\n";
		# Now, let's parse the public.site.index, the backbone
		if(open(my $PSI,'<:encoding(UTF-8)',$localIndexPath)) {
			my %indexConfig = (
				TabParser::TAG_HAS_HEADER	=> 1,
				TabParser::TAG_CONTEXT	=> $publicIndexPayload,
				TabParser::TAG_FETCH_COLS => PUBLIC_INDEX_COLS,
				TabParser::TAG_CALLBACK => \&public_results_callback,
			);
			$indexConfig{TabParser::TAG_VERBOSE} = 1  if($testmode);
			TabParser::parseTab($PSI,%indexConfig);
			close($PSI);
		} else {
			Carp::croak("Unable to parse $localIndexPath, the main metadata holder");
		}

		print "Parsing ",$dataIndex,"...\n";
		# Now, let's parse the data_files.index, the backbone
		if(open(my $DFI,'<:encoding(UTF-8)',$localDataFilesIndexPath)) {
			my %indexConfig = (
				TabParser::TAG_HAS_HEADER	=> 1,
				TabParser::TAG_CONTEXT	=> $publicIndexPayload,
				TabParser::TAG_FETCH_COLS => PUBLIC_INDEX_COLS,
				TabParser::TAG_POS_FILTER	=> [['FILE_TYPE' => 'BS_METH_TABLE_CYTOSINES_CNAG']],
				TabParser::TAG_CALLBACK => \&data_files_callback,
			);
			$indexConfig{TabParser::TAG_VERBOSE} = 1  if($testmode);
			TabParser::parseTab($DFI,%indexConfig);
			close($DFI);
		} else {
			Carp::croak("Unable to parse $localDataFilesIndexPath, the accessory metadata holder");
		}

		
		# For each data model
		foreach my $loadModelName (@loadModels) {
			print "Storing data using $loadModelName mapper\n";
			my $mapper = $storageModels{$loadModelName};
			
			# Now, do we need to push the metadata there?
			if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
				if($testmode) {
					print "\t [TESTMODE]Skipping storage of metadata model\n";
				} else {
					print "\t* Storing native model\n";
					$mapper->storeNativeModel();
				}
			}
			
			# Several hacks in a row... Yuck!
			if(!defined($modelDomain) || $modelDomain eq 'sdata') {
				my $conceptDomain = $model->getConceptDomain('sdata');
				print "Storing ",$conceptDomain->fullname,"\n";

				my %corrConcepts = map { $_ => BP::Loader::CorrelatableConcept->new($conceptDomain->conceptHash->{$_}) } keys(%{$conceptDomain->conceptHash});
				
				my $destination = undef;
				my $bulkData = undef;
				my $entorp = undef;
				my @bulkArray = ();
				
				# donor
				$mapper->setDestination($corrConcepts{'donor'});
				
				@bulkArray = values(%donors);
				$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
				unless($testmode) {
					$destination = $mapper->getInternalDestination();
					$bulkData = $mapper->_bulkPrepare($entorp);
					$mapper->_bulkInsert($destination,$bulkData);
				} else {
					print "[TESTMODE] Skipping storage of donors\n";
				}
				
				$destination = undef;
				$mapper->freeDestination();
				@bulkArray = ();
				$bulkData = undef;
				$entorp = undef;
				
				# specimen
				$mapper->setDestination($corrConcepts{'specimen'});
				
				@bulkArray = values(%specimens);
				$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
				unless($testmode) {
					$destination = $mapper->getInternalDestination();
					$bulkData = $mapper->_bulkPrepare($entorp);
					$mapper->_bulkInsert($destination,$bulkData);
				} else {
					print "[TESTMODE] Skipping storage of specimens\n";
				}
				
				$destination = undef;
				$mapper->freeDestination();
				@bulkArray = ();
				$bulkData = undef;
				$entorp = undef;
				
				# sample
				$mapper->setDestination($corrConcepts{'sample'});
				
				@bulkArray = values(%samples);
				$entorp = $mapper->validateAndEnactEntry(\@bulkArray);
				unless($testmode) {
					$destination = $mapper->getInternalDestination();
					$bulkData = $mapper->_bulkPrepare($entorp);
					$mapper->_bulkInsert($destination,$bulkData);
				} else {
					print "[TESTMODE] Skipping storage of samples\n";
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
				
				my @modelDomains = defined($modelDomain)?($DOMAIN2EXPANAL{$modelDomain}) : values(%DOMAIN2EXPANAL);
				
				foreach my $p_modelDomain (@modelDomains) {
					my($expDomain,$p_analDomains) = @{$p_modelDomain};
					
					if(exists($lab{$expDomain})) {
						my $destination = undef;
						my $bulkData = undef;
						my $entorp = undef;

						print "Storing $labFullname\n\t* ",$labConceptDomain->conceptHash->{$expDomain}->fullname,"...\n";
						$mapper->setDestination(BP::Loader::CorrelatableConcept->new($labConceptDomain->conceptHash->{$expDomain}));
						$entorp = $mapper->validateAndEnactEntry($lab{$expDomain});
						unless($testmode) {
							$destination = $mapper->getInternalDestination();
							$bulkData = $mapper->_bulkPrepare($entorp);
							$mapper->_bulkInsert($destination,$bulkData);
						} else {
							print "\t[TESTMODE] Skipping storage of IHEC experiment data $labFullname\n";
						}
						
						$destination = undef;
						$mapper->freeDestination();
						$bulkData = undef;
						$entorp = undef;
						
						foreach my $analDomain (@{$p_analDomains}) {
							if(exists($anal{$analDomain})) {
								my $conceptDomain = $model->getConceptDomain($analDomain);
								my %corrConcepts = map { $_ => BP::Loader::CorrelatableConcept->new($conceptDomain->conceptHash->{$_}) } keys(%{$conceptDomain->conceptHash});
								
								print "Storing ",$conceptDomain->fullname,"\n";
								
								if(exists($corrConcepts{'m'})) {
									print "\t* ",$corrConcepts{'m'}->concept->fullname,"...\n";
									$mapper->setDestination($corrConcepts{'m'});
									$entorp = $mapper->validateAndEnactEntry($anal{$analDomain});
									unless($testmode) {
										$destination = $mapper->getInternalDestination();
										$bulkData = $mapper->_bulkPrepare($entorp);
										$mapper->_bulkInsert($destination,$bulkData);
									} else {
										print "\t[TESTMODE] Skipping storage of analysis metadata ".$corrConcepts{'m'}->concept->fullname."\n";
									}
									
									$destination = undef;
									$mapper->freeDestination();
									$bulkData = undef;
									$entorp = undef;
									
									# And here the different bulk load
									if(exists($primary_anal{$analDomain})) {
										foreach my $p_primary (@{$primary_anal{$analDomain}}) {
											my($analysis_id,$conceptName,$instance,$remote_file) = @{$p_primary};
											print "\t* ",$corrConcepts{$conceptName}->concept->fullname," ($remote_file)...\n";
											
											$mapper->setDestination($corrConcepts{$conceptName});
											
											my $p_remote_files = (ref($remote_file) eq 'ARRAY')?$remote_file:[$remote_file];
											
											foreach my $r_file (@{$p_remote_files}) {
												my $local_file = $workingDir->cachedGet($bpDataServer,$r_file);
												
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
														unless($testmode) {
															eval {
																$instance->insert($F,$analysis_id,$mapper);
															};
															
															if($@) {
																Carp::carp("Errors while processing $remote_file: ".$@);
															}
														} else {
															print "\t[TESTMODE] Skipping storage of ".$corrConcepts{$conceptName}->concept->fullname." ($local_file)\n";
														}
														close($F);
													} else {
														Carp::carp("File $local_file (fetched from $remote_file) not processed. Reason: ".$!);
													}
													
													# At the end, free space of the huge downloaded file
													unlink($local_file);
												} else {
													Carp::carp("File $remote_file not processed (unable to fetch it). Reason: ".$bpDataServer->message);
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
	} elsif(!defined($localIndexPath)) {
		Carp::croak("FATAL ERROR: Unable to fetch index from $indexPath (host $host)");
	}
	
	$bpDataServer->disconnect()  if($bpDataServer->can('disconnect'));
	$bpDataServer->quit()  if($bpDataServer->can('quit'));
	
} else {
	print STDERR "Usage: $0 [-t] iniFile cachingDir [",join('|','sdata',sort(keys(%DOMAIN2EXPANAL))),"]\n"
}
