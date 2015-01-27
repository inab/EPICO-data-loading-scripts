#!/usr/bin/perl -w

use strict;

use diagnostics;

use Config::IniFiles;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";

use Carp;
use File::Basename;
use File::Path;
use File::Spec;
use IO::Handle;
use Net::FTP::AutoReconnect;
#use Net::SFTP::Foreign 1.76;
use Tie::IxHash;
use XML::LibXML::Reader;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Autoload::Relational;
use BP::Loader::Mapper::Autoload::Elasticsearch;
use BP::Loader::Mapper::Autoload::MongoDB;

use TabParser;

use constant DCC_LOADER_SECTION => 'dcc-loader';

use constant {
	PUBLIC_INDEX	=>	'public.results.index',
	DATA_FILES_INDEX	=>	'data_files.index',
	EXPERIMENTS2DATASETS	=>	'experiments2datasets.txt'
};

use constant PUBLIC_INDEX_COLS => [
	# sdata_donor
	'DONOR_ID',	# donor_id
	'CELL_LINE',	# alternate donor_id
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
	'Male'	=>	'm',
	'Female'	=>	'f',
);

my %GROUPCV = (
	'CNAG'	=>	['11'],
	'NCMLS'	=>	['1'],
	'NCMLS_CU'	=>	['1','3b'],
	'MPIMG'	=>	['12d'],
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
);

#####
# Parser methods
# --------------
# Each method must take these parameters
#	F: A filehandler with the content
#	analysis_id: The analysis_id for each entry
#	mapper: A BP::Loader::Mapper instance
#####

sub macsBedParser($$$);

sub rnaGFFQuantParser($$$);

sub dsHotspotsBedParser($$$);

sub dlatBedHypoMParser($$$);
sub dlatBedHyperMParser($$$);
sub __dlatBedParser($$$$);

sub dlatTxtCpGParser($$$);



use constant UNK_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'',
		'url'	=>	'',
	},
	'other_analysis_algorithm'	=>	[
	],
};

use constant WIGGLER_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'',
		'url'	=>	'',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'WIGGLER',
			'url'	=>	'https://code.google.com/p/align2rawsignal/'
		}
	],
};

use constant CS_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
		{
			'program'	=>	'BWA',
			'version'	=>	'0.5.9',
		},
		{
			'program'	=>	'samtools',
			'version'	=>	'0.1.18',
		},
		{
			'program'	=>	'phantompeakqualtools',
			'version'	=>	'1.1',
		},
		{
			'program'	=>	'spp',
			'version'	=>	'1.11',
		},
		{
			'program'	=>	'macs2',
			'version'	=>	'2.0.10.20120913',
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'BWA',
		'url'	=>	'http://bio-bwa.sourceforge.net/',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'Samtools',
			'url'	=>	'http://samtools.sourceforge.net',
		},
		{
			'name'	=>	'PhantomPeakQualTools',
			'url'	=>	'http://code.google.com/p/phantompeakqualtools/',
		},
		{
			'name'	=>	'Macs2',
			'url'	=>	'https://pypi.python.org/pypi/MACS2',
		}
	],
};

use constant CBR_METADATA => {
	'assembly_version'	=>	8,
	'program_versions'	=>	[
		{
			'program'	=>	'Bowtie',
			'version'	=>	'0.12.8'
		},
		{
			'program'	=>	'MMSEQ',
			'version'	=>	'1.0.5'
		}
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'Bowtie',
		'url'	=>	'http://bowtie-bio.sourceforge.net/index.shtml',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'MMSEQ',
			'url'	=>	'https://github.com/eturro/mmseq',
		}
	],
};

use constant CRG_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEMTools library',
		'url'	=>	'http://github.com/gemtools',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'Samtools',
			'url'	=>	'http://samtools.sourceforge.net',
		},
		{
			'name'	=>	'Cufflinks',
			'url'	=>	'http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_input',
		},
		{
			'name'	=>	'Flux Capacitor',
			'url'	=>	'http://sammeth.net/confluence/display/FLUX/Home',
		}
	],
};

use constant DS_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
		{
			'program'	=>	'BWA',
			'version'	=>	'0.5.9'
		},
		{
			'program'	=>	'samtools',
			'version'	=>	'0.1.18'
		},
		{
			'program'	=>	'phantompeakqualtools',
			'version'	=>	'1.1'
		},
		{
			'program'	=>	'spp',
			'version'	=>	'1.11'
		},
		{
			'program'	=>	'macs2',
			'version'	=>	'2.0.10.20120913'
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'BWA',
		'url'	=>	'http://bio-bwa.sourceforge.net/',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'samtools',
			'url'	=>	'http://samtools.sourceforge.net/',
		},
		{
			'name'	=>	'Hotspot',
			'url'	=>	'http://www.uwencode.org/proj/hotspot-ptih/',
		},
		{
			'name'	=>	'Bedops',
			'url'	=>	'http://code.google.com/p/bedops/',
		},
	],
};

use constant METH_CPG_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
		{
			'program'	=>	'GEM',
			'version'	=>	'1.242',
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEM',
		'url'	=>	'http://big.crg.cat/services/gem_genome_multi_tool_library',
	},
	'other_analysis_algorithm'	=>	[
	],
	'mr_type'	=>	'cpg',
};

use constant METH_HYPER_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
		{
			'program'	=>	'GEM',
			'version'	=>	'1.242',
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEM',
		'url'	=>	'http://big.crg.cat/services/gem_genome_multi_tool_library',
	},
	'other_analysis_algorithm'	=>	[
	],
	'mr_type'	=>	'hyper',
};

use constant METH_HYPO_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
		{
			'program'	=>	'GEM',
			'version'	=>	'1.242',
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEM',
		'url'	=>	'http://big.crg.cat/services/gem_genome_multi_tool_library',
	},
	'other_analysis_algorithm'	=>	[
	],
	'mr_type'	=>	'hypo',
};




use constant {
	F_DOMAIN	=>	0,
	F_ANALYSIS_GROUP_ID	=>	1,
	F_POSTFIX	=>	2,
	F_PATTERN_POSTFIX	=>	3,
	F_PRIMARY	=>	4,
	F_PARSER	=>	5,
	F_METADATA	=>	6,
	F_PARENT_POSTFIX	=>	7,
};

use constant {
	CNAG_CPGS_POSTFIX	=>	'bs_cpg',
	CNAG_CYTOSINES_POSTFIX	=>	'bs_c'
};

# 0. The concept domain
# 1. The group
# 2. The postfix to the experiment_id in order to build the analysis_id
#
my %FILETYPE2ANAL = (
	'CS_BROAD_MACS2'	=>	[
		'pdna',
		['1'],
		undef,
		[['bed.gz' => 'cs_broad_peaks']],
		'p',
		\&macsBedParser,
		CS_METADATA,
		undef
	],
	
	'CS_MACS2'	=>	[
		'pdna',
		['1'],
		undef,
		[['bed.gz' => 'cs_peaks']],
		'p',
		\&macsBedParser,
		CS_METADATA,
		undef
	],
	
	'CS_WIGGLER'	=>	[
		'pdna',
		['1'],
		'cs_wiggler',
		undef,
		undef,
		undef,
		WIGGLER_METADATA,
		undef
	],
	
	'RNA_GENE_QUANT_CBR'	=>	[
		'exp',
		['3b'],
		undef,
		[['.gff.gz' => 'gq_cbr']],
		'g',
		\&rnaGFFQuantParser,
		CBR_METADATA,
		undef
	],
	
	'RNA_GENE_QUANT_CRG'	=>	[
		'exp',
		['15b'],
		undef,
		[['.gff' => 'gq_crg']],
		'g',
		\&rnaGFFQuantParser,
		CRG_METADATA,
		undef
	],
	
	'RNA_TRANSCRIPT_QUANT_CBR'	=>	[
		'exp',
		['3b'],
		undef,
		[['.gff.gz' => 'tq_cbr']],
		't',
		\&rnaGFFQuantParser,
		CBR_METADATA,
		undef
	],
	
	'RNA_TRANSCRIPT_QUANT_CRG'	=>	[
		'exp',
		['15b'],
		undef,
		[['.gtf' => 'tq_crg']],
		't',
		\&rnaGFFQuantParser,
		CRG_METADATA,
		undef
	],
	
	'RNA_JUNCTIONS_CRG'	=>	[
		'jcn',
		['15b'],
		'junctions',
		undef,
		undef,
		undef,
		CRG_METADATA,
		undef
	],
	
	'DS_HOTSPOT'	=>	[
		'rreg',
		['8'],
		'ds_hotspots',
		[['peaks' => 'ds_hotspots_peaks']],
		'p',
		\&dsHotspotsBedParser,
		DS_METADATA,
		undef
	],
	
	'DS_WIGGLER'	=>	[
		'rreg',
		['8'],
		'ds_wiggler',
		undef,
		undef,
		undef,
		WIGGLER_METADATA,
		undef
	],
	
	'BS_HYPER_METH_BED_CNAG'	=>	[
		'dlat',
		['11'],
		'bs_hyper',
		undef,
		'mr',
		\&dlatBedHyperMParser,
		METH_HYPER_METADATA,
		CNAG_CPGS_POSTFIX
	],
	
	'BS_HYPO_METH_BED_CNAG'	=>	[
		'dlat',
		['11'],
		'bs_hypo',
		undef,
		'mr',
		\&dlatBedHypoMParser,
		METH_HYPO_METADATA,
		undef	#CNAG_CPGS_POSTFIX
	],
	
	'BS_METH_TABLE_CYTOSINES_CNAG'	=>	[
		'dlat',
		['11'],
		undef,
		[['cpgs.bs_call' => +CNAG_CPGS_POSTFIX]],
		'cpg',
		\&dlatTxtCpGParser,
		METH_CPG_METADATA,
		undef	#CNAG_CYTOSINES_POSTFIX
	],
);

######
# Shouldn't be global variables, but we need them in the different callbacks
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

my %reg_analysis = ();

# Constants to access the data on each domain in primary_anal
use constant {
	P_ANALYSIS_ID	=>	0,
	P_CONCEPTNAME	=>	1,
	P_METHOD	=>	2,
	P_FILE	=>	3,
};

my %primary_anal = ();

my %expfiles = ();

# Correspondence between experiment ids and EGA ids
my %exp2EGA = ();

# Cache for cell lines
my %cellSpecimenTerm = ();
my %cellPurifiedTerm = ();

# Parameters needed
my $bpDataServer = undef;
my $metadataPath = undef;
my $cachingDir = undef;

#####
# Method prototypes
#####
sub cachedGet($$$);
sub parseIHECsample($$$$);
sub parseIHECexperiment($$$$);

#####
# Method callbacks
#####

sub experiments_to_datasets_callback(\%$$$) {
	my($p_exp2EGA,$experiment_id, $title, $dataset_id) = @_;

	$p_exp2EGA->{$experiment_id} = $dataset_id;
}

my $ensembl_version = undef;
my $gencode_version = undef;

sub data_files_callback {
	my(
		$donor_id,
		$cell_line,
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
	
	if(exists($experiments{$experiment_id})) {
		# And this is the analysis metadata
		if(exists($FILETYPE2ANAL{$file_type})) {
			my $ftype = $FILETYPE2ANAL{$file_type};
			my $analDomain = $ftype->[F_DOMAIN];
			
			# Analysis id building
			my $an_postfix = undef;
			
			if(defined($ftype->[F_PATTERN_POSTFIX])) {
				foreach my $p_pat_post (@{$ftype->[F_PATTERN_POSTFIX]}) {
					my($pattern,$postfix)=@{$p_pat_post};
					
					if(index($remote_file_path,$pattern)!=-1) {
						$an_postfix = $postfix;
						last;
					}
				}
			}
			
			$an_postfix = $ftype->[F_POSTFIX]  unless(defined($an_postfix));
			
			# No postfix, no processing!!!!
			if(defined($an_postfix)) {
				my $analysis_id = $experiment_id.'_'.$an_postfix;
				
				unless(exists($reg_analysis{$analysis_id})) {
					my $f_metadata = $ftype->[F_METADATA];
					$f_metadata = {}  unless(defined($f_metadata));
					
					my %analysis = (
						'analysis_id'	=>	$analysis_id,
						'experiment_id'	=>	$experiment_id,
						'analysis_group_id'	=>	$ftype->[F_ANALYSIS_GROUP_ID],
						'data_status'	=>	defined($ftype->[F_PRIMARY])?2:0,
						'assembly_version'	=>	$f_metadata->{assembly_version},
						'ensembl_version'	=>	$ensembl_version,
						'gencode_version'	=>	$gencode_version,
					);
					$analysis{NSC} = $NSC  if($NSC ne '-');
					$analysis{RSC} = $RSC  if($RSC ne '-');
					@analysis{keys(%{$f_metadata})} = values(%{$f_metadata});
					if(defined($ftype->[F_PARENT_POSTFIX])) {
						$analysis{'base_analysis_id'} = $analysis_id.'_'.$ftype->[F_PARENT_POSTFIX];
					}
					
					# Last, register it!
					$anal{$analDomain} = []  unless(exists($anal{$analDomain}));
					push(@{$anal{$analDomain}},\%analysis);
					$reg_analysis{$analysis_id} = undef;
				}
				
				# Preparing the field
				if(defined($ftype->[F_PRIMARY])) {
					$primary_anal{$analDomain} = []  unless(exists($primary_anal{$analDomain}));
					
					push(@{$primary_anal{$analDomain}},[$analysis_id,$ftype->[F_PRIMARY],$ftype->[F_PARSER],$remote_file_path]);
				}
			}
		}
	}
}

sub public_results_callback {
	my(
		$donor_id,
		$cell_line,
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
	
	$donor_id = $cell_line  if($donor_id eq '-');
	
	unless(exists($donors{$donor_id})) {
		$donor_ethnicity = undef  if($donor_ethnicity eq 'NA' || $donor_ethnicity eq '-');
		
		my $donor_region_of_residence_term = '001';	# World
		$donor_region_of_residence_term = 'ALIAS:EAL'  if($donor_region_of_residence eq "East Anglia");
		
		my %donor = (
			'donor_id'	=>	$donor_id,
			'donor_sex'	=>	exists($SEXCV{$donor_sex})?$SEXCV{$donor_sex}:undef,
			'donor_region_of_residence'	=>	[$donor_region_of_residence_term],
			'donor_ethnicity'	=>	$donor_ethnicity,
			'notes'	=> ($cell_line ne '-')?'Cell line':undef,
		);
		$donors{$donor_id} = \%donor;
	}
	
	my $specimen_id = $sample_id . '_spec';
	my $p_IHECsample = undef;
	unless(exists($specimens{$specimen_id})) {
		$donor_health_status = undef  if($donor_health_status eq 'NA' || $donor_health_status eq '-');
		$specimen_processing = undef  if($specimen_processing eq 'NA' || $specimen_processing eq '-');
		$specimen_storage = undef  if($specimen_storage eq 'NA' || $specimen_storage eq '-');
		
		my $donor_min_age_at_specimen_acquisition = undef;
		my $donor_max_age_at_specimen_acquisition = undef;
		if($donor_age =~ /(\d+)\s+-\s+(\d+)/) {
			$donor_min_age_at_specimen_acquisition = 'P'.$1.'Y';
			$donor_max_age_at_specimen_acquisition = 'P'.$2.'Y';
		} elsif($donor_age =~ /(\d+)\s+weeks/) {
			$donor_min_age_at_specimen_acquisition = $donor_max_age_at_specimen_acquisition = 'P'.$1.'M';
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
			if ($tissue_type eq "Peripheral blood"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0013756";
			} elsif($tissue_type eq "Cord blood"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0012168";
			} elsif($tissue_type eq "Tonsil"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0002372";
			} elsif($tissue_type eq "Bone marrow"){
				$specimen_term = "http://purl.obolibrary.org/obo/UBERON_0002371";
			}
		}
		
		# Last resort, look at the cache
		if(defined($specimen_term)) {
			$cellSpecimenTerm{$cell_line} = $specimen_term  if($cell_line ne '-' && !exists($cellSpecimenTerm{$cell_line}));
		} elsif($cell_line ne '-' && exists($cellSpecimenTerm{$cell_line})) {
			$specimen_term = $cellSpecimenTerm{$cell_line};
		}
		
		Carp::croak("Undefined specimen term for $specimen_id!!!!")  unless(defined($specimen_term));
		
		$p_IHECsample = parseIHECsample($bpDataServer,$metadataPath,$sample_id,$cachingDir);
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
		$p_IHECsample = parseIHECsample($bpDataServer,$metadataPath,$sample_id,$cachingDir)  unless(defined($p_IHECsample));
		
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
			$cellPurifiedTerm{$cell_line} = $purified_cell_type  if($cell_line ne '-' && !exists($cellPurifiedTerm{$cell_line}));
		} elsif($cell_line ne '-' && exists($cellPurifiedTerm{$cell_line})) {
			$purified_cell_type = $cellPurifiedTerm{$cell_line};
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
		$samples{$sample_id} = \%sample;
	}
	
	if(exists($EXPERIMENTCV{$library_strategy})) {
		# This is the experimental metadata
		unless(exists($experiments{$experiment_id})) {
			my $labexp = $EXPERIMENTCV{$library_strategy};
			
			my($p_IHECexperiment,$ihec_library_strategy,$ihec_instrument_model) = parseIHECexperiment($bpDataServer,$metadataPath,$experiment_id,$cachingDir);
			
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
									'accession'	=>	exists($exp2EGA{$experiment_id})?$exp2EGA{$experiment_id}:'',
									'url'	=>	exists($exp2EGA{$experiment_id})?('https://www.ebi.ac.uk/ega/datasets/'.$exp2EGA{$experiment_id}):'',
								},
				'platform'	=>	exists($INSTRUMENT2PLATFORM{$instrument_model})?$INSTRUMENT2PLATFORM{$instrument_model}:-1,
				'platform_model'	=>	$ihec_instrument_model,
				'seq_coverage'	=>	undef,
				'extraction_protocol'	=>	exists($p_IHECexperiment->{EXTRACTION_PROTOCOL})?$p_IHECexperiment->{EXTRACTION_PROTOCOL}[0]:undef,
			);
			
			# Last, register it!
			$lab{$labexp} = []  unless(exists($lab{$labexp}));
			push(@{$lab{$labexp}},\%experiment);
			$experiments{$experiment_id} = undef;
		}
		
		&data_files_callback(@_);
		
	} else {
		Carp::carp("Unknown type of experiment: ".$library_strategy);
	}
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	F: A filehandler with the content
#	analysis_id: The analysis_id for each entry
#	mapper: A BP::Loader::Mapper instance
#####

sub macsBedParser($$$) {
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->{'batch-size'};
	
	my $numBatch = 0;
	my @batch = ();
	
	my %macsBedParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				$chromosome_start,
				$chromosome_end,
				$protein_dna_interaction_id,
				undef,
				undef,
				$fold_enrichment,	# fold_enrichment
				$log10_pvalue, # -log10(pvalue)
				$log10_qvalue, # -log10(qvalue)
			) = @_;
			
			my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
			
			$chromosome = 'MT'  if($chromosome eq 'M');
			
			my $protein_stable_id = ($protein_dna_interaction_id =~ /^[^.]+\.([^.]+)/)?$1:'';
			
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'protein_dna_interaction_id'	=>	$protein_dna_interaction_id,
				'chromosome'	=>	$chromosome,
				'chromosome_start'	=>	$chromosome_start+1,	# Bed holds the data 0-based
				'chromosome_end'	=>	$chromosome_end,	# Bed holds the end coordinate as exclusive, so it does not change
				'rank'	=>	[
					{
						'rank'	=>	'fold_enrichment',
						'value'	=>	$fold_enrichment
					}
				],
				'protein_stable_id'	=>	$protein_stable_id,
				'log10_pvalue'	=>	$log10_pvalue,
				'log10_qvalue'	=>	$log10_qvalue,
			);
			
			push(@batch,\%entry);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				$mapper->bulkInsert(\@batch);
				
				@batch = ();
				$numBatch = 0;
			}
		},
	);
	TabParser::parseTab($F,%macsBedParserConfig);
	
	# Last step
	if($numBatch > 0) {
		$mapper->bulkInsert(\@batch);
		
		@batch = ();
	}
}

sub rnaGFFQuantParser($$$) {
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->{'batch-size'};
	
	my $numBatch = 0;
	my @batch = ();
	
	my %rnaGFFQuantParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				undef, # source
				$feature, # feature
				$chromosome_start,
				$chromosome_end,
				undef,
				$chromosome_strand,
				undef, # frame
				$attributes_str, # attributes following .ace format
			) = @_;
			
			my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
			
			$chromosome = 'MT'  if($chromosome eq 'M');
			
			my %attributes = ();
			
			my @tokens = split(/\s*;\s*/,$attributes_str);
			foreach my $token (@tokens) {
				my($key,$value) = split(/\s+/,$token,2);
				
				# Removing double quotes
				$value =~ tr/"//d;
				
				$attributes{$key} = $value;
			}
			
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'chromosome'	=>	$chromosome,
				'chromosome_start'	=>	$chromosome_start,
				'chromosome_end'	=>	$chromosome_end,
				'chromosome_strand'	=>	(($chromosome_strand eq '-')?-1:1),
				'normalized_read_count'	=>	$attributes{'RPKM'},
				'raw_read_count'	=>	int($attributes{'reads'} + 0.5),
				'is_annotated'	=>	1,
			);
			$entry{'gene_stable_id'} = $attributes{'gene_id'}  if(exists($attributes{'gene_id'}));
			$entry{'transcript_stable_id'} = $attributes{'transcript_id'}  if($feature eq 'transcript' && exists($attributes{'transcript_id'}));
			
			push(@batch,\%entry);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				$mapper->bulkInsert(\@batch);
				
				@batch = ();
				$numBatch = 0;
			}
		},
	);
	TabParser::parseTab($F,%rnaGFFQuantParserConfig);
	
	# Last step
	if($numBatch > 0) {
		$mapper->bulkInsert(\@batch);
		
		@batch = ();
	}
}

sub dsHotspotsBedParser($$$) {
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->{'batch-size'};
	
	my $numBatch = 0;
	my @batch = ();
	
	my %dsHotspotsBedParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				$chromosome_start,
				$chromosome_end,
				$zscore,
				$zscore_peak,
			) = @_;
			
			my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
			
			$chromosome = 'MT'  if($chromosome eq 'M');
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'chromosome'	=>	$chromosome,
				'chromosome_start'	=>	$chromosome_start+1,	# Bed holds the data 0-based
				'chromosome_end'	=>	$chromosome_end,	# Bed holds the end coordinate as exclusive, so it does not change
				'z_score'	=>	defined($zscore_peak)?$zscore_peak:$zscore,
			);
			
			push(@batch,\%entry);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				$mapper->bulkInsert(\@batch);
				
				@batch = ();
				$numBatch = 0;
			}
		},
	);
	TabParser::parseTab($F,%dsHotspotsBedParserConfig);
	
	# Last step
	if($numBatch > 0) {
		$mapper->bulkInsert(\@batch);
		
		@batch = ();
	}
}

sub dlatBedHyperMParser($$$) {
	return __dlatBedParser($_[0],$_[1],$_[2],'hyper');
}

sub dlatBedHypoMParser($$$) {
	return __dlatBedParser($_[0],$_[1],$_[2],'hypo');
}

sub __dlatBedParser($$$$) {
	my($F,$analysis_id,$mapper,$hyperhypo) = @_;
	
	# UGLY
	my $BMAX = $mapper->{'batch-size'};
	
	my $numBatch = 0;
	my @batch = ();
	
	my %dlatBedParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				$chromosome_start,
				$chromosome_end,
				undef,	# Size of region in base pairs
				$avg_meth_level,	# Average methylation level in region
				undef,	# Number of CpGs in region
				$d_lated_reads,	# Median number of non-converted reads at CpGs in region
				$converted_reads,	# Median number of converted reads at CpGs in region
				$total_reads,	# Median number of total reads at CpGs in region
				undef,	# Island/Shelf/Shore (union of CpG Island annotations for all CpGs in region)
				undef,	# refGene annotation (union of refGene  annotations for all CpGs in region)
			) = @_;
			
			my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
			
			$chromosome = 'MT'  if($chromosome eq 'M');
			
			$chromosome_start = $chromosome_start+1;	# Bed holds the data 0-based
			my $d_lated_fragment_id = $hyperhypo.'|'.$chro.'_'.$chromosome_start.'_'.$chromosome_end;
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'd_lated_fragment_id'	=>	$d_lated_fragment_id,
				'chromosome'	=>	$chromosome,
				'chromosome_start'	=>	$chromosome_start,
				'chromosome_end'	=>	$chromosome_end,	# Bed holds the end coordinate as exclusive, so it does not change
				'total_reads'	=>	$total_reads,
				'c_total_reads'	=>	($d_lated_reads + $converted_reads),
				'd_lated_reads'	=>	$d_lated_reads,
				'meth_level'	=>	$avg_meth_level
				
			);
			
			push(@batch,\%entry);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				$mapper->bulkInsert(\@batch);
				
				@batch = ();
				$numBatch = 0;
			}
		},
	);
	TabParser::parseTab($F,%dlatBedParserConfig);
	
	# Last step
	if($numBatch > 0) {
		$mapper->bulkInsert(\@batch);
		
		@batch = ();
	}
}

sub dlatTxtCpGParser($$$) {
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->{'batch-size'};
	
	my $numBatch = 0;
	my @batch = ();
	
	my %dlatTxtCpGParserConfig = (
		TabParser::TAG_CALLBACK => sub {

			my(
				$chro,			# Chromosome
				$chromosome_start,	# Position of CC (offset 1)
				$probability,		# Phred scaled probability of genotype *not* being CC/GG
				$avg_meth_level,	# Methylation probability (combined estimate from the weighted average of the MLEs at the two positions)
				undef,			# Standard deviation of methylation probability (from weighted average)
				$d_lated_reads,		# No. of non-converted C reads (sum of counts at both positions)
				$converted_reads,	# No. of converted C reads (idem)
				undef,			# Total reads supporting genotype call (idem)
				$total_reads,		# Total reads (idem)
			) = @_;
			
			my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
			
			$chromosome = 'MT'  if($chromosome eq 'M');
			
			my $chromosome_end = $chromosome_start+1;	# As it is a CpG, it is one more
			my $d_lated_fragment_id = 'cpg|'.$chro.'_'.$chromosome_start.'_'.$chromosome_end;
			
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				#'d_lated_fragment_id'	=>	$d_lated_fragment_id,
				'chromosome'	=>	$chromosome,
				'chromosome_start'	=>	$chromosome_start,	# This txt had the coordinates 1-based
				#'chromosome_end'	=>	$chromosome_end,
				'total_reads'	=>	$total_reads,
				'c_total_reads'	=>	($d_lated_reads + $converted_reads),
				'd_lated_reads'	=>	$d_lated_reads,
				'meth_level'	=>	$avg_meth_level,
				'probability'	=>	$probability
			);
			
			push(@batch,\%entry);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				$mapper->bulkInsert(\@batch);
				
				@batch = ();
				$numBatch = 0;
			}
		},
	);
	TabParser::parseTab($F,%dlatTxtCpGParserConfig);
	
	# Last step
	if($numBatch > 0) {
		$mapper->bulkInsert(\@batch);
		
		@batch = ();
	}
}

#####
# Method bodies
#####
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
		$remotePath = '/'.$remotePath  unless(substr($remotePath,0,1) eq '/');
		File::Path::make_path($localDir);
		#print STDERR join(" -=- ",$remotePath,$cachingDir,$localPath,$localBasePath,$localRelDir,$localDir),"\n";
		my $targetLocalPath = $localPath;
		$localPath = $bpDataServer->get($remotePath,$localPath);
		print STDERR "DEBUGFTP: ($remotePath -> $targetLocalPath) ".$bpDataServer->message."\n"  unless(defined($localPath));
		utime($filedate,$filedate,$localPath)  if(defined($localPath));
	}
	
	return $localPath;
}

sub parseIHECsample($$$$) {
	my($bpDataServer,$metadataPath,$sample_id,$cachingDir) = @_;
	
	print "\t* Parsing IHEC sample $sample_id...\n";
	
	my $localIHECsample = cachedGet($bpDataServer,join('/',$metadataPath,'samples',substr($sample_id,0,6),$sample_id.'.xml'),$cachingDir);
	
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
	my($bpDataServer,$metadataPath,$experiment_id,$cachingDir) = @_;
	
	print "\t* Parsing IHEC experiment $experiment_id...\n";
	
	my $localIHECexperiment = cachedGet($bpDataServer,join('/',$metadataPath,'experiments',substr($experiment_id,0,6),$experiment_id.'.xml'),$cachingDir);
	
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
	$cachingDir = shift(@ARGV);
	my $modelDomain = shift(@ARGV);
	
	Carp::croak('ERROR: Unknown knowledge domain '.$modelDomain)  if(defined($modelDomain) && $modelDomain ne 'sdata' && !exists($DOMAIN2EXPANAL{$modelDomain}));
	
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
	
	my $localIndexPath = cachedGet($bpDataServer,$indexPath.'/'.PUBLIC_INDEX,$cachingDir);
	my $localDataFilesIndexPath = cachedGet($bpDataServer,$indexPath.'/'.DATA_FILES_INDEX,$cachingDir);
	my $localExp2Datasets = cachedGet($bpDataServer,$indexPath.'/'.EXPERIMENTS2DATASETS,$cachingDir);
	
	if(defined($localIndexPath) && defined($localExp2Datasets)) {
		# Try getting a connection to 
		
		# Let's parse the model
		my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
		# Setting up the right path on relative cases
		$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

		print "Parsing model $modelFile...\n";
		my $model = undef;
		eval {
			$model = BP::Model->new($modelFile);
			$ensembl_version = exists($model->annotations->hash->{GENCODEVer})?$model->annotations->hash->{GENCODEVer}:'';
			$gencode_version = exists($model->annotations->hash->{GENCODEVer})?$model->annotations->hash->{EnsemblVer}:'';
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
		
		# First, these correspondences experiment <=> EGA needed by next parse
		print "Parsing ",EXPERIMENTS2DATASETS,"...\n";
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
		
		print "Parsing ",PUBLIC_INDEX,"...\n";
		# Now, let's parse the public.site.index, the backbone
		if(open(my $PSI,'<:encoding(UTF-8)',$localIndexPath)) {
			my %indexConfig = (
				TabParser::TAG_HAS_HEADER	=> 1,
				TabParser::TAG_FETCH_COLS => PUBLIC_INDEX_COLS,
				TabParser::TAG_CALLBACK => \&public_results_callback,
			);
			$indexConfig{TabParser::TAG_VERBOSE} = 1  if($testmode);
			TabParser::parseTab($PSI,%indexConfig);
			close($PSI);
		} else {
			Carp::croak("Unable to parse $localIndexPath, the main metadata holder");
		}

		print "Parsing ",DATA_FILES_INDEX,"...\n";
		# Now, let's parse the public.site.index, the backbone
		if(open(my $DFI,'<:encoding(UTF-8)',$localDataFilesIndexPath)) {
			my %indexConfig = (
				TabParser::TAG_HAS_HEADER	=> 1,
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
											my($analysis_id,$conceptName,$method,$remote_file) = @{$p_primary};
											print "\t* ",$corrConcepts{$conceptName}->concept->fullname," ($remote_file)...\n";
											
											$mapper->setDestination($corrConcepts{$conceptName});
											
											my $p_remote_files = (ref($remote_file) eq 'ARRAY')?$remote_file:[$remote_file];
											
											foreach my $r_file (@{$p_remote_files}) {
												my $local_file = cachedGet($bpDataServer,$r_file,$cachingDir);
												
												if(defined($local_file)) {
													my $f_mode = undef;
													my @f_params = ();
													
													# Compressed file detection
													if($local_file =~ /\.gz$/) {
														$f_mode = '-|';
														push(@f_params,'gunzip','-c',$local_file);
													} else {
														$f_mode = '<';
														push(@f_params,$local_file);
													}
													
													if(open(my $F,$f_mode,@f_params)) {
														unless($testmode) {
															eval {
																$method->($F,$analysis_id,$mapper);
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
	} else {
		Carp::croak("FATAL ERROR: Unable to fetch experiments to datasets correspondence from $indexPath (host $host)");
	}
	
	$bpDataServer->disconnect()  if($bpDataServer->can('disconnect'));
	$bpDataServer->quit()  if($bpDataServer->can('quit'));
	
} else {
	print STDERR "Usage: $0 [-t] iniFile cachingDir [",join('|','sdata',keys(%DOMAIN2EXPANAL)),"]\n"
}
