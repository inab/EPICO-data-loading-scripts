#!/usr/bin/perl

use warnings;
use strict;

use Carp;
use BP::DCCLoader::Parsers;

package BP::DCCLoader::Parsers::BlueprintDCCMetadataParser;

use File::Basename;
use File::Spec;
# Give autoflush to STDOUT and STDERR
use IO::Handle;
use Log::Log4perl;
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

##
# Needed specific parsers must be registered BEFORE the instance is created
##
#use BP::DCCLoader::Parsers::AbstractInsertionParser;
#use BP::DCCLoader::Parsers::CpGInsertionParser;
#use BP::DCCLoader::Parsers::DNASEBedInsertionParser;
#use BP::DCCLoader::Parsers::MACSBedInsertionParser;
#use BP::DCCLoader::Parsers::MethRegionsBedInsertionParser;
#use BP::DCCLoader::Parsers::RNASeqGFFInsertionParser;
#use BP::DCCLoader::Parsers::RNASeqStarInsertionParser;
#use BP::DCCLoader::Parsers::WigglerInsertionParser;

use BP::DCCLoader::Parsers::AbstractInsertionParser;
use BP::DCCLoader::WorkingDir;

use TabParser;

##
# Constants and static structures
##
use constant {
	PUBLIC_INDEX_DEFAULT	=>	'public.results.index',
	DATA_INDEX_DEFAULT	=>	'data_files.index',
	DEFAULT_FTP_REL	=>	'/',
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
	['SAMPLE_BARCODE','SAMPLE_DESC_2'],	# specimen_biomaterial_id ???
	# specimen_available is unknown
	# donor_id is already got from 'DONOR_ID'
	
	# sdata_sample
	'SAMPLE_ID',	# sample_id
	'SAMPLE_NAME',	# sample_name
	'SAMPLE_ONTOLOGY_URI',	# purified_cell_type, it contains the URI of CellOntology term
	# culture_conditions comes from IHEC metadata (CULTURE_CONDITIONS)
	# markers comes from IHEC metadata (MARKERS)
	# analyzed_sample_type must be set to 11 (other) to simplify
	['SAMPLE_DESCRIPTION','SAMPLE_DESC_3'],	# analyzed_sample_type_other
	# analyzed_sample_interval is unknown
	# specimen_id is already derived from 'SAMPLE_ID'
	
	# lab_*
	'STUDY_ID',	# study_id
	'STUDY_NAME',	# study_name
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
	'FILE_SIZE',	# The expected file size
	'FILE_MD5',	# The expected MD5 sum
	'NSC',	# NSC
	'RSC',	# RSC
	#'WITHDRAWN',	# Has been withdrawn this analysis?
];

use constant NORMAL_DISEASE_URI	=>	'http://purl.obolibrary.org/obo/PATO_0000461';

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

use constant {
	UNSPECIFIED_PLATFORM	=>	-1,
	UNSPECIFIED_PLATFORM_MODEL	=>	'unspecified',
};
my %INSTRUMENT2PLATFORM = (
	UNSPECIFIED_PLATFORM_MODEL()	=>	UNSPECIFIED_PLATFORM,
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
sub parseIHECsample($$$$$);
sub parseIHECexperiment($$$$$);

#####
# Method callbacks
#####

sub experiments_to_datasets_callback(\%$$$) {
	my($p_exp2EGA,$experiment_id, $title, $dataset_id) = @_;
	
	# The file could contain spureous contentse
	$p_exp2EGA->{$experiment_id} = $dataset_id  if(index($title,'GRCh38')!=-1);
}

sub datasets_parse_callback(\%$) {
	my($payload,$experiment_id) = @_;

	$payload->{exp2EGA}->{$experiment_id} = $payload->{dataset_id};
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
		$sample_name,
		$purified_cell_type_uri,
		$analyzed_sample_type_other,
		
		$study_id,
		$study_name,
		
		$experiment_id,
		$center_name,
		$library_strategy,
		$instrument_model,
		
		$file_type,
		$remote_file_path,
		$file_size,
		$file_md5,
		$NSC,
		$RSC,
	)=@_;
	
	
	if(exists($payload->{experiments}{$experiment_id})) {
		# And this is the analysis metadata
		if(exists($payload->{FILETYPE2ANAL}{$file_type})) {
			my $ftypeArr = $payload->{FILETYPE2ANAL}{$file_type};
			
			# Iterate over all the registered filetypes, until one matches
			foreach my $ftype (@{$ftypeArr}) {
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
				
				# No postfix, skip! No processing!!!!
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
						
						my $remote_file = {
							'r_file'	=>	$remote_file_path,
							'expectedSize'	=>	$file_size,
							'expectedMD5'	=>	$file_md5
						};
						push(@{$payload->{primary_anal}{$analDomain}},[$analysis_id,$ftype->[BP::DCCLoader::Parsers::F_PRIMARY],$ftype->[BP::DCCLoader::Parsers::F_PARSER],$remote_file]);
					}
					last;
				}
			}
		}
	}
}

# Remember to put in lowercase the keys
my %TISSUE_MAP = (
	"peripheral blood"	=>	"http://purl.obolibrary.org/obo/UBERON_0013756",
	"venous blood"		=>	"http://purl.obolibrary.org/obo/UBERON_0013756",
	"cord blood"		=>	"http://purl.obolibrary.org/obo/UBERON_0012168",
	"tonsil"		=>	"http://purl.obolibrary.org/obo/UBERON_0002372",
	"bone marrow"		=>	"http://purl.obolibrary.org/obo/UBERON_0002371",
	"thymus"		=>	"http://purl.obolibrary.org/obo/UBERON_0002370",
	"thymus lymphoid tissue"		=>	"http://purl.obolibrary.org/obo/UBERON_0003483",
	"liver"			=>	"http://purl.obolibrary.org/obo/UBERON_0002107",
);

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
		$sample_name,
		$purified_cell_type_uri,
		$analyzed_sample_type_other,
		
		$study_id,
		$study_name,
		
		$experiment_id,
		$center_name,
		$library_strategy,
		$instrument_model,
		
		$file_type,
		$remote_file_path,
		$file_size,
		$file_md5,
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
			$payload->{LOG}->logdie("Unable to identify the kind of donor");
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
		
		# Tweaks
		
		# First pass
		$donor_disease_uri =~ s/&amp;/&/g;
		# Second pass
		$donor_disease_uri =~ s/&amp;/&/g;
		$donor_disease_uri =~ s/%2520/%20/g;
		
		if($donor_disease_uri eq '-' || $donor_disease_uri eq 'NA') {
			$donor_disease_uri = NORMAL_DISEASE_URI;
		} elsif(index($donor_disease_uri,'http')!=0) {
			print STDERR "FIXME: $specimen_id $donor_disease_uri\n";
			$donor_disease_uri = undef;
		}
		
		unless(defined($donor_disease_uri)) {
			if($donor_disease_text eq 'None') {
				$donor_disease_uri = NORMAL_DISEASE_URI;
			}
			
			print STDERR "TODO: empty donor disease: $specimen_id $donor_disease_text\n"  unless(defined($donor_disease_uri));
		}
		
		my $specimen_term = undef;
		
		my @purified_term_uris = split(/;/,$purified_cell_type_uri);
		
		foreach my $term_uri (@purified_term_uris) {
			if($term_uri =~ /obo\/(?:((?:UBERON)|(?:CLO))_([^\/]+))/ || $term_uri =~ /efo\/(EFO)_([^\/]+)/) {
				$specimen_term = $term_uri;
				last;
			}
		}
		
		unless(defined($specimen_term)) {
			my $lc_tissue_type = lc($tissue_type);
			$specimen_term = $TISSUE_MAP{$lc_tissue_type}  if(exists($TISSUE_MAP{$lc_tissue_type}));
		}
		
		# Last resort, look at the cache
		if(defined($specimen_term)) {
			$payload->{cellSpecimenTerm}{$cell_line} = $specimen_term  if($cell_line ne '-' && !exists($payload->{cellSpecimenTerm}{$cell_line}));
		} elsif($cell_line ne '-' && exists($payload->{cellSpecimenTerm}{$cell_line})) {
			$specimen_term = $payload->{cellSpecimenTerm}{$cell_line};
		}
		
		$payload->{LOG}->logdie("Undefined specimen term for $specimen_id!!!!")  unless(defined($specimen_term));
		
		$p_IHECsample = parseIHECsample($payload->{bpMetadataServer},$payload->{metadataPath},$sample_id,$payload->{workingDir},$payload->{LOG});
		my %specimen = (
			'specimen_id'	=>	$specimen_id,
			'tissue_type'	=>	$tissue_type,
			'tissue_depot'	=>	$tissue_depot,
			'specimen_term'	=>	$specimen_term,
			'collection_method'	=>	exists($p_IHECsample->{COLLECTION_METHOD})?$p_IHECsample->{COLLECTION_METHOD}:undef,
			'donor_min_age_at_specimen_acquisition'	=>	$donor_min_age_at_specimen_acquisition,
			'donor_max_age_at_specimen_acquisition'	=>	$donor_max_age_at_specimen_acquisition,
			'donor_health_status'	=>	$donor_health_status,
			'donor_disease'	=>	$donor_disease_uri,
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
		$p_IHECsample = parseIHECsample($payload->{bpMetadataServer},$payload->{metadataPath},$sample_id,$payload->{workingDir},$payload->{LOG})  unless(defined($p_IHECsample));
		
		my %sampleFeatures = map { $_ => { 'feature' => $_ , 'value' => $p_IHECsample->{$_} } } keys(%{$p_IHECsample});
		
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
			'sample_name'	=>	$sample_name,
			'purified_cell_type'	=>	$purified_cell_type,
			'culture_conditions'	=>	exists($p_IHECsample->{CULTURE_CONDITIONS})?$p_IHECsample->{CULTURE_CONDITIONS}:undef,
			'markers'	=>	exists($p_IHECsample->{MARKERS})?$p_IHECsample->{MARKERS}:undef,
			'analyzed_sample_type'	=>	11,
			'analyzed_sample_type_other'	=>	$analyzed_sample_type_other,
			'analyzed_sample_interval'	=>	undef,
			'specimen_id'	=>	$specimen_id,
			'features'	=>	\%sampleFeatures,
		);
		$payload->{samples}{$sample_id} = \%sample;
	}
	
	if(exists($EXPERIMENTCV{$library_strategy})) {
		unless(exists($payload->{studies}{$study_id})) {
			my %study = (
				'study_id'	=>	$study_id,
				'study_name'	=>	$study_name
			);
			
			$payload->{studies}{$study_id} = \%study;
		}
		
		# This is the experimental metadata
		unless(exists($payload->{experiments}{$experiment_id})) {
			my $labexp = $EXPERIMENTCV{$library_strategy};
			
			my($p_IHECexperiment,$ihec_library_strategy,$ihec_instrument_model) = parseIHECexperiment($payload->{bpMetadataServer},$payload->{metadataPath},$experiment_id,$payload->{workingDir},$payload->{LOG});
			
			my %features = map { my $val = { 'feature' => $_ , 'value' => $p_IHECexperiment->{$_}[0] }; $val->{'units'} = $p_IHECexperiment->{$_}[1]  if(defined($p_IHECexperiment->{$_}[1])); $_ => $val } keys(%{$p_IHECexperiment});
			
			# Curating instrument_model
			$ihec_instrument_model = $instrument_model  if(!defined($ihec_instrument_model) && length($instrument_model) > 0);
			unless(defined($ihec_instrument_model)) {
				$ihec_instrument_model = UNSPECIFIED_PLATFORM_MODEL;
				$payload->{LOG}->logwarn("FIXME: On experiment $experiment_id, missing instrument model");
			}
			my $ihec_instrument = exists($INSTRUMENT2PLATFORM{$ihec_instrument_model})?$INSTRUMENT2PLATFORM{$ihec_instrument_model}:UNSPECIFIED_PLATFORM;
			
			# The common attributes
			my %experiment = (
				'experiment_id'	=>	$experiment_id,
				'analyzed_sample_id'	=>	$sample_id,
				'study_id'	=>	$study_id,
				'experiment_type'	=>	exists($p_IHECexperiment->{EXPERIMENT_TYPE})?$p_IHECexperiment->{EXPERIMENT_TYPE}[0]:'',
				'library_strategy'	=>	defined($ihec_library_strategy)?$ihec_library_strategy:$library_strategy,
				'experimental_group_id'	=>	exists($GROUPCV{$center_name})?$GROUPCV{$center_name}:[$center_name],
				'features'	=>	\%features,
				'raw_data_repository'	=>	1,
				'raw_data_accession'	=>	{
									'accession'	=>	exists($payload->{exp2EGA}{$experiment_id})?$payload->{exp2EGA}{$experiment_id}:'',
									'url'	=>	exists($payload->{exp2EGA}{$experiment_id})?('https://www.ebi.ac.uk/ega/datasets/'.$payload->{exp2EGA}{$experiment_id}):'',
								},
				'platform'	=>	$ihec_instrument,
				'platform_model'	=>	$ihec_instrument_model,
				'seq_coverage'	=>	undef,
				'extraction_protocol'	=>	exists($p_IHECexperiment->{EXTRACTION_PROTOCOL})?$p_IHECexperiment->{EXTRACTION_PROTOCOL}[0]:($payload->{testmode}?'':undef),
			);
			
			$payload->{LOG}->logwarn("FIXME: On experiment $experiment_id, missing extraction protocol")  unless(exists($p_IHECexperiment->{EXTRACTION_PROTOCOL}));
			
			# Last, register it!
			$payload->{lab}{$labexp} = []  unless(exists($payload->{lab}{$labexp}));
			push(@{$payload->{lab}{$labexp}},\%experiment);
			$payload->{experiments}{$experiment_id} = undef;
		}
		
		&data_files_callback(@_);
		
	} else {
		$payload->{LOG}->logwarn("Unknown type of experiment: ".$library_strategy);
	}
}

#####
# Method bodies
#####
sub parseIHECsample($$$$$) {
	my($bpDataServer,$metadataPath,$sample_id,$workingDir,$LOG) = @_;
	
	$LOG->info("\t* Parsing IHEC sample $sample_id...");
	
	my($localIHECsample, $reason) = $workingDir->cachedGet($bpDataServer,join('/',$metadataPath,'samples',substr($sample_id,0,6),$sample_id.'.xml'));
	
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
		$LOG->logwarn("Unable to fetch metadata file about sample $sample_id . Reason: ".$reason);
	}
	
	return \%IHECsample;
}

sub parseIHECexperiment($$$$$) {
	my($bpDataServer,$metadataPath,$experiment_id,$workingDir,$LOG) = @_;
	
	$LOG->info("\t\t* Parsing IHEC experiment $experiment_id...");
	
	my($localIHECexperiment, $reason) = $workingDir->cachedGet($bpDataServer,join('/',$metadataPath,'experiments',substr($experiment_id,0,6),$experiment_id.'.xml'));
	
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
		$LOG->logwarn("Unable to fetch metadata file about experiment $experiment_id . Reason: ".$reason);
	}
	
	return (\%IHECexperiment,$library_strategy,$instrument_model);
}

##############
# Constructor
############
sub new($$) {
	my $self = shift;
	my $class = ref($self) || $self;
	
	my($iniFile, $cachingDir) = @_;
	
	my $LOG = Log::Log4perl->get_logger(__PACKAGE__);
	
	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
	
	# Check the needed parameters for the SSH connection
	my $protocol = undef;
	my $host = undef;
	my $user = undef;
	my $pass = undef;
	my $blueprintFTPRel = undef;
	my $indexPath = undef;

	my $metadataProtocol = undef;
	my $metadataHost = undef;
	my $metadataUser = undef;
	my $metadataPass = undef;
	my $blueprintMetadataFTPRel = undef;
	my $metadataIndexPath = undef;
	
	# Defined outside
	my $metadataPath = undef;
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'protocol')) {
		$protocol = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'protocol');
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'protocol'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'host')) {
		$host = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'host');
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'host'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'user')) {
		$user = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'user');
	} else {
		$user = 'ftp'  if($protocol eq 'ftp');
		$LOG->logdie("Configuration file $iniFile must have 'user'")  unless(defined($user));
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'pass')) {
		$pass = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'pass');
	} else {
		$pass = ($user eq 'ftp')?'guest@':''  if($protocol eq 'ftp');
		$LOG->logdie("Configuration file $iniFile must have 'pass'")  unless(defined($pass));
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'blueprint-rel')) {
		$blueprintFTPRel = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'blueprint-rel');
	} else {
		$blueprintFTPRel = DEFAULT_FTP_REL;
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'index-path')) {
		$indexPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'index-path');
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'index-path'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-protocol')) {
		$metadataProtocol = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-protocol');
	} else {
		$metadataProtocol = $protocol;
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-host')) {
		$metadataHost = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-host');
	} elsif($metadataProtocol eq $protocol) {
		$metadataHost = $host;
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'metadata-host'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-user')) {
		$metadataUser = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-user');
	} elsif($metadataProtocol eq $protocol && $metadataHost eq $host) {
		$metadataUser = $user;
	} else {
		$metadataUser = 'ftp'  if($metadataProtocol eq 'ftp');
		$LOG->logdie("Configuration file $iniFile must have 'metadata-user'")  unless(defined($metadataUser));
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-pass')) {
		$metadataPass = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-pass');
	} elsif($metadataProtocol eq $protocol && $metadataHost eq $host && $metadataUser eq $user) {
		$metadataPass = $pass;
	} else {
		$metadataPass = ($metadataUser eq 'ftp')?'guest@':''  if($metadataProtocol eq 'ftp');
		$LOG->logdie("Configuration file $iniFile must have 'metadata-pass'")  unless(defined($metadataPass));
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'blueprint-metadata-rel')) {
		$blueprintMetadataFTPRel = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'blueprint-metadata-rel');
	} elsif($metadataProtocol eq $protocol && $metadataHost eq $host && $metadataUser eq $user) {
		$blueprintMetadataFTPRel = $blueprintFTPRel;
	} else {
		$blueprintMetadataFTPRel = DEFAULT_FTP_REL;
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-index-path')) {
		$metadataIndexPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-index-path');
	} elsif($metadataProtocol eq $protocol && $metadataHost eq $host && $metadataUser eq $user) {
		$metadataIndexPath = $indexPath;
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'metadata-index-path'");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-path')) {
		$metadataPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'metadata-path');
	} else {
		$LOG->logdie("Configuration file $iniFile must have 'metadata-path'");
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
	}
	
	my $egaDataSetsPath;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'egadatasets-path')) {
		$egaDataSetsPath = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,'egadatasets-path');
	}
	
	# Use the internal file as reference in case of the EGA datasets path not defined
	unless(defined($egaDataSetsPath) || defined($exp2datasets)) {
		$exp2datasets = EXPERIMENTS2DATASETS_DEFAULT;
	}
	
	$self = {
		'iniFile'	=>	$iniFile,
		'ini'	=>	$ini,
		'cachingDir'	=>	$cachingDir,
		'LOG'	=>	$LOG,
		
		# The needed parameters for the connections
		'protocol'	=>	$protocol,
		'host'	=>	$host,
		'user'	=>	$user,
		'pass'	=>	$pass,
		'blueprintFTPRel'	=>	$blueprintFTPRel,
		'indexPath'	=>	$indexPath,

		'metadataProtocol'	=>	$metadataProtocol,
		'metadataHost'	=>	$metadataHost,
		'metadataUser'	=>	$metadataUser,
		'metadataPass'	=>	$metadataPass,
		'blueprintMetadataFTPRel'	=>	$blueprintMetadataFTPRel,
		'metadataIndexPath'	=>	$metadataIndexPath,
		
		# Defined outside
		'metadataPath'	=>	$metadataPath,
		
		'publicIndex'	=>	$publicIndex,
		'dataIndex'	=>	$dataIndex,
		'exp2datasets'	=>	$exp2datasets,
		'egaDataSetsPath'	=>	$egaDataSetsPath,
	};
	
	$self = bless($self,$class);
	
	return $self;
}

sub getConfigIni() {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	return $self->{ini};
}

sub getWorkingDir() {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	unless(exists($self->{workingDir})) {
		# First, explicitly create the caching directory
		my $workingDir = BP::DCCLoader::WorkingDir->new($self->{cachingDir});

		$self->{workingDir} = $workingDir;
	}
	
	return $self->{workingDir};
}

sub getModel() {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	unless(exists($self->{model})) {
		my $ini = $self->{ini};
		my $iniFile = $self->{iniFile};
		
		# Let's parse the model
		my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
		# Setting up the right path on relative cases
		$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

		$LOG->info("Parsing model $modelFile...");
		my $model = undef;
		my $ensembl_version = undef;
		my $gencode_version = undef;
		eval {
			$model = BP::Model->new($modelFile);
			$ensembl_version = exists($model->annotations->hash->{EnsemblVer})?$model->annotations->hash->{EnsemblVer}:'';
			$gencode_version = exists($model->annotations->hash->{GENCODEVer})?$model->annotations->hash->{GENCODEVer}:'';
		};
		
		if($@) {
			$LOG->logdie('ERROR: Model parsing and validation failed. Reason: '.$@);
		}
		$LOG->info("\tDONE!");
		
		# Saving the model and some additional parameters
		$self->{model} = $model;
		$self->{ensembl_version} = $ensembl_version;
		$self->{gencode_version} = $gencode_version;
	}
	
	return $self->{model};
}

sub getLoadAndStorageModels() {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	unless(exists($self->{p_loadModels})) {
		# Let's parse the model
		my $ini = $self->{ini};
		my $model = $self->getModel();
		
		my %storageModels = ();
		
		# Setting up the loader storage model(s)
		$LOG->logdie('ERROR: undefined destination storage model')  unless($ini->exists($BP::Loader::Mapper::SECTION,'loaders'));
		my $loadModelNames = $ini->val($BP::Loader::Mapper::SECTION,'loaders');
		
		my @loadModels = ();
		foreach my $loadModelName (split(/,/,$loadModelNames)) {
			unless(exists($storageModels{$loadModelName})) {
				$storageModels{$loadModelName} = BP::Loader::Mapper->newInstance($loadModelName,$model,$ini);
				push(@loadModels,$loadModelName);
			}
		}
		
		$self->{p_loadModels} = \@loadModels;
		$self->{p_storageModels} = \%storageModels;
	}
	
	return ($self->{p_loadModels},$self->{p_storageModels});
}

sub getPublicIndexPayload($) {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	unless(exists($self->{p_publicIndexPayload})) {
		my($testmode) = @_;
		
		my $ini = $self->{ini};
		# A BP::DCCLoader::WorkingDir instance
		my $workingDir = $self->getWorkingDir();
		
		
		
		######
		# These were global variables, but as we need them in the different callbacks
		######
		my %donors = ();
		tie(%donors,'Tie::IxHash');
		
		my %specimens = ();
		tie(%specimens,'Tie::IxHash');
		
		my %samples = ();
		tie(%samples,'Tie::IxHash');
		
		my %studies = ();
		tie(%studies,'Tie::IxHash');
		
		my %experiments = ();
		
		# Laboratory experiments and analysis metadata
		my %lab = ();
		my %anal = ();
		
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
		my $bpMetadataServer = undef;
		
		my $ensembl_version = undef;
		my $gencode_version = undef;



		# Check the needed parameters for the SSH connection
		my $protocol = $self->{'protocol'};
		my $host = $self->{'host'};
		my $user = $self->{'user'};
		my $pass = $self->{'pass'};
		my $blueprintFTPRel = $self->{'blueprintFTPRel'};
		my $indexPath = $self->{'indexPath'};

		my $metadataProtocol = $self->{'metadataProtocol'};
		my $metadataHost = $self->{'metadataHost'};
		my $metadataUser = $self->{'metadataUser'};
		my $metadataPass = $self->{'metadataPass'};
		my $blueprintMetadataFTPRel = $self->{'blueprintMetadataFTPRel'};
		my $metadataIndexPath = $self->{'metadataIndexPath'};
		
		# Defined outside
		my $metadataPath = $self->{'metadataPath'};

		my $publicIndex = $self->{'publicIndex'};
		my $dataIndex = $self->{'dataIndex'};
		my $exp2datasets = $self->{'exp2datasets'};
		my $egaDataSetsPath = $self->{'egaDataSetsPath'};






		
		$LOG->info("Connecting to host $protocol://$host...");
		# Defined outside
		$bpDataServer = $workingDir->doBPConnect($protocol,$host,$user,$pass);
		
		$bpMetadataServer = undef;
		if($metadataProtocol eq $protocol && $metadataHost eq $host && $metadataUser eq $user && $metadataPass eq $pass) {
			# Connection sharing
			$bpMetadataServer = $bpDataServer;
		} else {
			$LOG->info("Connecting to metadata host $metadataProtocol://$metadataHost...");
			$bpMetadataServer = $workingDir->doBPConnect($metadataProtocol,$metadataHost,$metadataUser,$metadataPass);
		}
		
		my($localIndexPath, $indexPathReason) = $workingDir->cachedGet($bpDataServer,join('/',$blueprintFTPRel , $indexPath , $publicIndex));
		my($localDataFilesIndexPath, $dataFilePathReason) = $workingDir->cachedGet($bpMetadataServer,join('/',$blueprintMetadataFTPRel , $metadataIndexPath , $dataIndex));
		
		if(defined($localIndexPath)) {

			# Let's parse the model
			my $model = $self->getModel();
			
			# Initializations of parsable file types, needed by next parsing tasks
			my $p_FILETYPE2ANAL = BP::DCCLoader::Parsers->getParsableFiletypes({
				BP::DCCLoader::Parsers::AbstractInsertionParser::K_INI		=>	$ini,
				BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL	=>	$model,
				BP::DCCLoader::Parsers::AbstractInsertionParser::K_WORKINGDIR	=>	$workingDir,
				BP::DCCLoader::Parsers::AbstractInsertionParser::K_TESTMODE	=>	$testmode,
			});

			# First, these correspondences experiment <=> EGA needed by next parse
			if(defined($egaDataSetsPath)) {
				$LOG->info("Parsing datasets from $egaDataSetsPath...");
				my $egaDataSetsRemotePath = join('/',$blueprintFTPRel,$egaDataSetsPath);
				my $p_egaDataSetsFiles = $workingDir->listing($bpDataServer,$egaDataSetsRemotePath);
				
				if(defined($p_egaDataSetsFiles)) {
					# Let's fetch all the dataset definitions
					foreach my $dataSetPath (@{$p_egaDataSetsFiles}) {
						if($dataSetPath =~ /(EGA[DS][0-9]*)_[^\/]+_analysis_files.tsv$/) {
							my $dataset_id = $1;
							
							my($localDataSetFile, $dataSetFileReason) = $workingDir->cachedGet($bpDataServer,$dataSetPath);
							if(defined($localDataSetFile)) {
								$LOG->info("* Parsing $dataSetPath...");
								
								if(open(my $E2DF,'<:encoding(UTF-8)',$localDataSetFile)) {
									my %e2dfConfig = (
										TabParser::TAG_HAS_HEADER	=> 1,
										TabParser::TAG_CONTEXT	=> {
											'exp2EGA'	=>	\%exp2EGA,
											'dataset_id'	=>	$dataset_id
										},
										TabParser::TAG_FETCH_COLS => ['EXPERIMENT_ID'],
										TabParser::TAG_CALLBACK => \&datasets_parse_callback,
									);
									$e2dfConfig{TabParser::TAG_VERBOSE} = 1  if($testmode);
									TabParser::parseTab($E2DF,%e2dfConfig);
									close($E2DF);
								} else {
									$LOG->warn("Unable to parse $localDataSetFile, needed to get the correspondence to EGA dataset identifiers. Reason: ".$!);
								}
							} else {
								$LOG->warn("Unable to fetch experiments to datasets correspondence from $dataSetPath (host $host). Reason: $dataSetFileReason");
							}
						}
					}
				} else {
					$LOG->warn("$egaDataSetsPath unavailable. raw_data_accession subfields will be empty");
				}
			} else {
				my($localExp2Datasets, $exp2datasetsReason) = $workingDir->cachedGet($bpMetadataServer,join('/',$blueprintMetadataFTPRel , $metadataIndexPath , $exp2datasets));
				
				if(defined($localExp2Datasets)) {
					$LOG->info("Parsing $exp2datasets...");
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
						$LOG->logdie("Unable to parse $localExp2Datasets, needed to get the EGA dataset identifiers. Reason: ".$!);
					}
				} else {
					$LOG->warn("Unable to fetch experiments to datasets correspondence from $blueprintMetadataFTPRel $metadataIndexPath (host $metadataHost). Reason: $exp2datasetsReason");
					$LOG->warn("$exp2datasets unavailable. raw_data_accession subfields will be empty");
				}
			}
			
			my $publicIndexPayload = {
				'donors'	=>	\%donors,
				'specimens'	=>	\%specimens,
				'samples'	=>	\%samples,
				'studies'	=>	\%studies,
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
				'bpMetadataServer'	=>	$bpMetadataServer,
				'blueprintFTPRel'	=>	$blueprintFTPRel,
				'blueprintMetadataFTPRel'	=>	$blueprintMetadataFTPRel,
				'metadataPath'	=>	join('/',$blueprintMetadataFTPRel,$metadataPath),
				'ensembl_version'	=>	$ensembl_version,
				'gencode_version'	=>	$gencode_version,
				'testmode'	=>	$testmode,
				'LOG'	=>	$LOG,
			};
			
			$LOG->info("Parsing $publicIndex...");
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
				$LOG->logdie("Unable to parse $localIndexPath, the main metadata holder");
			}
			
			if(defined($localDataFilesIndexPath)) {
				$LOG->info("Parsing $dataIndex...");
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
					$LOG->logdie("Unable to parse $localDataFilesIndexPath, the accessory metadata holder");
				}
			} else {
				$LOG->warn("Unable to fetch the accessory metadata holder from $blueprintMetadataFTPRel $metadataIndexPath $dataIndex (host $metadataHost). Reason: ".$dataFilePathReason)  unless(defined($dataIndex));

			}
			
			$self->{p_publicIndexPayload} = $publicIndexPayload;
		} elsif(!defined($localIndexPath)) {
			$LOG->logdie("FATAL ERROR: Unable to fetch index from $indexPath (host $host). Reason: $indexPathReason");
		}
	}
	
	return $self->{p_publicIndexPayload};
	
}

sub dataServerGet($$$) {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my($remote_file_path,$expectedSize,$expectedMD5) = @_;
	
	my $publicIndexPayload = $self->getPublicIndexPayload();
	
	return $publicIndexPayload->{workingDir}->cachedGet($publicIndexPayload->{bpDataServer},join('/',$publicIndexPayload->{blueprintFTPRel},$remote_file_path),$expectedSize,$expectedMD5);
}

sub disconnect() {
	my $self = shift;
	
	my $LOG = $self->{LOG};
	
	$LOG->logdie((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	# Disconnecting the FTP servers
	if(exists($self->{p_publicIndexPayload})) {
		my $publicIndexPayload = $self->{p_publicIndexPayload};
		
		if(exists($publicIndexPayload->{bpMetadataServer})) {
			my $bpDataServer = $publicIndexPayload->{bpDataServer};
			my $bpMetadataServer = $publicIndexPayload->{bpMetadataServer};
			
			if($bpMetadataServer != $bpDataServer) {
				$bpMetadataServer->disconnect()  if($bpMetadataServer->can('disconnect'));
				$bpMetadataServer->quit()  if($bpMetadataServer->can('quit'));
			}
			
			$bpDataServer->disconnect()  if($bpDataServer->can('disconnect'));
			$bpDataServer->quit()  if($bpDataServer->can('quit'));
			
			delete($publicIndexPayload->{bpDataServer});
			delete($publicIndexPayload->{bpMetadataServer});
		}
	}
}

1;
