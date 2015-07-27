#!/usr/bin/perl

use warnings;
use strict;

use TabParser;

package BP::DCCLoader::Parsers::MACSBedInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

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

# This is the empty constructor
sub new(;$$) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new()  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	return {
		'CHIP_MACS2_BROAD_BED'	=>	[
			'pdna',
			['1'],
			undef,
			[['bed.gz' => 'cs_broad_peaks']],
			'p',
			$self,
			CS_METADATA,
			undef
		],
		
		'CHIP_MACS2_BED'	=>	[
			'pdna',
			['1'],
			undef,
			[['bed.gz' => 'cs_peaks']],
			'p',
			$self,
			CS_METADATA,
			undef
		],
	};
	
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	F: A filehandler with the content
#	analysis_id: The analysis_id for each entry
#	mapper: A BP::Loader::Mapper instance
#####
sub insert($$$) {
	my($self)=shift;
	
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->bulkBatchSize();
	
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

1;
