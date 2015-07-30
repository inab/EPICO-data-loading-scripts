#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use TabParser;

package BP::DCCLoader::Parsers::RNASeqGFFInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

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
		'RNA_GENE_QUANT_CBR'	=>	[
			'exp',
			['3b'],
			undef,
			[['.gff.gz' => 'gq_cbr']],
			'g',
			$self,
			CBR_METADATA,
			undef
		],
		
		'RNA_GENE_QUANT_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.gff' => 'gq_crg']],
			'g',
			$self,
			CRG_METADATA,
			undef
		],
		
		'RNA_TRANSCRIPT_QUANT_CBR'	=>	[
			'exp',
			['3b'],
			undef,
			[['.gff.gz' => 'tq_cbr']],
			't',
			$self,
			CBR_METADATA,
			undef
		],
		
		'RNA_TRANSCRIPT_QUANT_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.gtf' => 'tq_crg']],
			't',
			$self,
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
			
			$p_insertMethod->(\%entry);
		},
	);
	return \%rnaGFFQuantParserConfig;
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
