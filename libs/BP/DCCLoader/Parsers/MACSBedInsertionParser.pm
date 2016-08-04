#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use TabParser;

package BP::DCCLoader::Parsers::MACSBedInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant CS_METADATA => {
	'assembly_version'	=>	38,
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
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	my $MACS_BROAD_BED = [
		'pdna',
		['1'],
		undef,
		[['bed.gz' => 'cs_broad_peaks']],
		'p',
		$self,
		CS_METADATA,
		undef
	];
		
	my $MACS_BED = [
		'pdna',
		['1'],
		undef,
		[['bed.gz' => 'cs_peaks']],
		'p',
		$self,
		CS_METADATA,
		undef
	];

	return {
		# The two firsts are the old ones
		'CS_BROAD_MACS2'	=>	[ $MACS_BROAD_BED ],
		'CS_MACS2'	=>	[ $MACS_BED ],
		
		'CHIP_MACS2_BROAD_BED'	=>	[ $MACS_BROAD_BED ],
		'CHIP_MACS2_BED'	=>	[ $MACS_BED ],
	};
	
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$) {
	my($self)=shift;
	
	my($analysis_id,$p_insertMethod) = @_;
	
	my $model = $self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL};
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	my $termChro = '';
	my $term = undef;
	my $chromosome = undef;
	
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
			
			if($termChro ne $chro) {
				$termChro = $chro;
				$term = $chroCV->getTerm($chro);
				if($term) {
					$chromosome = $term->key();
				} else {
					print STDERR "\tdiscarding entries from unknown chromosome $chro\n";
					$chromosome = undef;
				}
			}
			
			if(defined($chromosome)) {
				my $protein_stable_id = ($protein_dna_interaction_id =~ /^[^.]+\.([^.]+)/)?$1:'';
				
				
				my %entry = (
					'analysis_id'	=>	$analysis_id,
					'protein_dna_interaction_id'	=>	$protein_dna_interaction_id,
					'chromosome'	=>	$chromosome,
					'chromosome_start'	=>	$chromosome_start+1,	# Bed holds the data 0-based
					'chromosome_end'	=>	$chromosome_end+0,	# Bed holds the end coordinate as exclusive, so it does not change
					'rank'	=>	[
						{
							'rank'	=>	'fold_enrichment',
							'value'	=>	$fold_enrichment+0.0
						}
					],
					'protein_stable_id'	=>	$protein_stable_id,
					'log10_pvalue'	=>	$log10_pvalue+0.0,
					'log10_qvalue'	=>	$log10_qvalue+0.0,
				);
				
				$p_insertMethod->(\%entry);
			}
		},
	);
	return \%macsBedParserConfig;
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
