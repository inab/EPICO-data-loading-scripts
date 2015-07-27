#!/usr/bin/perl

use warnings;
use strict;

use TabParser;

package BP::DCCLoader::Parsers::CpGInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

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

use constant CNAG_CPGS_POSTFIX	=>	'bs_cpg';
#use constant CNAG_CYTOSINES_POSTFIX	=>	'bs_c';

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
		'BS_METH_TABLE_CYTOSINES_CNAG'	=>	[
			'dlat',
			['11'],
			undef,
			[['cpgs.bs_call' => +CNAG_CPGS_POSTFIX]],
			'cpg',
			undef,	# Disabling the insertion parser (at least, register the CpGs)
			METH_CPG_METADATA,
			undef	#CNAG_CYTOSINES_POSTFIX
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

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
