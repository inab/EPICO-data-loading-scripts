#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use TabParser;

# For the reused constant
use BP::DCCLoader::Parsers::CpGInsertionParser;


# This package contains the common methods
package BP::DCCLoader::Parsers::MethRegionsBedInsertionParser;

sub __dlatBedParser($$$$);

sub __dlatBedParser($$$$) {
	my($analysis_id,$p_insertMethod,$chroCV,$hyperhypo) = @_;
	
	my $termChro = '';
	my $term = undef;
	my $chromosome = undef;
	
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
				$chromosome_start = $chromosome_start+1;	# Bed holds the data 0-based
				my $d_lated_fragment_id = $hyperhypo.'|'.$chro.'_'.$chromosome_start.'_'.$chromosome_end;
				
				my %entry = (
					'analysis_id'	=>	$analysis_id,
					'd_lated_fragment_id'	=>	$d_lated_fragment_id,
					'chromosome'	=>	$chromosome,
					'chromosome_start'	=>	$chromosome_start,
					'chromosome_end'	=>	$chromosome_end+0,	# Bed holds the end coordinate as exclusive, so it does not change
					'total_reads'	=>	$total_reads,
					'c_total_reads'	=>	($d_lated_reads + $converted_reads),
					'd_lated_reads'	=>	$d_lated_reads,
					'meth_level'	=>	$avg_meth_level+0e0,
					
				);
				
				$p_insertMethod->(\%entry);
			}
		},
	);
	return \%dlatBedParserConfig;
}


# Specific metadata for hypermethylation
package BP::DCCLoader::Parsers::MethRegionsBedInsertionParser::Hyper;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant METH_HYPER_METADATA => {
	'assembly_version'	=>	38,
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

# This is the empty constructor
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	return {
		'BS_HYPER_METH_BED_CNAG'	=>	[
			'dlat',
			['11'],
			'bs_hyper',
			undef,
			'mr',
			$self,
			METH_HYPER_METADATA,
			BP::DCCLoader::Parsers::CpGInsertionParser::CNAG_CPGS_POSTFIX
		],
	};
	
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$$) {
	my($self)=shift;
	
	my $model = $self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL};
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	
	return BP::DCCLoader::Parsers::MethRegionsBedInsertionParser::__dlatBedParser($_[0],$_[1],$chroCV,'hyper');
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);


package BP::DCCLoader::Parsers::MethRegionsBedInsertionParser::Hypo;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant METH_HYPO_METADATA => {
	'assembly_version'	=>	38,
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

# This is the empty constructor
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	return {
		'BS_HYPO_METH_BED_CNAG'	=>	[
			'dlat',
			['11'],
			'bs_hypo',
			undef,
			'mr',
			$self,
			METH_HYPO_METADATA,
			BP::DCCLoader::Parsers::CpGInsertionParser::CNAG_CPGS_POSTFIX
		],
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
	
	my $model = $self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL};
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	
	return BP::DCCLoader::Parsers::MethRegionsBedInsertionParser::__dlatBedParser($_[0],$_[1],$chroCV,'hypo');
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
