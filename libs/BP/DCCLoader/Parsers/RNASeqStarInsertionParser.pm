#!/usr/bin/perl

use warnings;
use strict;

use TabParser;

package BP::DCCLoader::Parsers::RNASeqStarInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant CRG_STAR_METADATA => {
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
		'RNA_GENE_QUANT_STAR_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.results' => 'gq_crg']],
			'g',
			$self,
			CRG_STAR_METADATA,
			undef
		],
		
		'RNA_TRANSCRIPT_QUANT_STAR_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.results' => 'tq_crg']],
			't',
			$self,
			CRG_STAR_METADATA,
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
	
	my %rnaGFFQuantParserConfig = (
		TabParser::TAG_HAS_HEADER	=>	1,
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

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
