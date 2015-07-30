#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use BP::DCCLoader::WorkingDir;
use TabParser;

package BP::DCCLoader::Parsers::RNASeqStarInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use BP::DCCLoader::Parsers::EnsemblGTParser;

use constant CRG_STAR_METADATA => {
	'assembly_version'	=>	38,
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
			'name'	=>	'Star/RSEM',
			'url'	=>	'http://deweylab.biostat.wisc.edu/rsem/',
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

sub getEnsemblCoordinates() {
	my $self = shift;
	
	unless(exists($self->{ENShash})) {
		$self->{ENShash} = BP::DCCLoader::Parsers::EnsemblGTParser::getEnsemblCoordinates($self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL},$self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_WORKINGDIR},$self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_INI});
	}
	
	return $self->{ENShash};
}

my @CommonKeysDelete = (
	'gene_id',
	'transcript_id(s)',
	'transcript_id',
	'length',
	'effective_length',
	'expected_count',
	'TPM',
	'FPKM',
);

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub commonInsert($$\@;$) {
	my $self = shift;
	
	my($analysis_id,$p_insertMethod,$p_colnames,$isTranscript) = @_;
	
	# Get the Ensembl hash
	my $p_ensHash = $self->getEnsemblCoordinates();
	
	my %rnaRSemStarGeneParserConfig = (
		TabParser::TAG_HAS_HEADER	=>	1,
		TabParser::TAG_FETCH_COLS	=>	$p_colnames,
		TabParser::TAG_CALLBACK => sub {
			my %metrics = ();
			@metrics{@{$p_colnames}} = @_;
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'gene_stable_id'	=>	$metrics{'gene_id'},
				'length'	=>	$metrics{'length'},
				'effective_length'	=>	$metrics{'effective_length'},
				'expected_count'	=>	$metrics{'expected_count'},
				'TPM'	=>	$metrics{'TPM'},
				'FPKM'	=>	$metrics{'FPKM'},
				'is_annotated'	=>	1,
				'metrics'	=>	\%metrics
			);
			
			# Choosing the ensId
			my $ensId;
			if($isTranscript) {
				$ensId = $entry{'transcript_stable_id'} = $metrics{'transcript_id'};
			} else {
				$ensId = $metrics{'gene_id'}
			}
			
			# Remove the revision from the ensId
			$ensId = substr($ensId,0,rindex($ensId,'.'));
			
			# Now, let's get the chromosomical coordinates
			if(exists($p_ensHash->{$ensId})) {
				my $p_data = $p_ensHash->{$ensId};
				$entry{'chromosome'} = $p_data->{'chromosome'};
				$entry{'chromosome_start'} = $p_data->{'chromosome_start'};
				$entry{'chromosome_end'} = $p_data->{'chromosome_end'};
				
				# Removing redundant data
				delete @metrics{@CommonKeysDelete};
				
				$p_insertMethod->(\%entry);
			} else {
				Carp::carp("ERROR: Next entry was rejected as its Ensembl Id does not match with the ones defined\n".join("\t",@{$p_colnames})."\n".join("\t",@_));
			}
		},
	);
	return \%rnaRSemStarGeneParserConfig;
}


package BP::DCCLoader::Parsers::RNASeqStarInsertionParser::Genes;

use base qw(BP::DCCLoader::Parsers::RNASeqStarInsertionParser);

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
		'RNA_GENE_QUANT_STAR_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.results' => 'gq_crg']],
			'g',
			$self,
			BP::DCCLoader::Parsers::RNASeqStarInsertionParser::CRG_STAR_METADATA,
			undef
		],
	};
	
}

my $p_colnamesGenes = [
	'gene_id',
	'transcript_id(s)',
	'length',
	'effective_length',
	'expected_count',
	'TPM',
	'FPKM',
	'posterior_mean_count',
	'posterior_standard_deviation_of_count',
	'pme_TPM',
	'pme_FPKM',
	'TPM_ci_lower_bound',
	'TPM_ci_upper_bound',
	'FPKM_ci_lower_bound',
	'FPKM_ci_upper_bound',
];

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$) {
	my $self = shift;
	
	return $self->commonInsert($_[0],$_[1],$p_colnamesGenes);
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);


package BP::DCCLoader::Parsers::RNASeqStarInsertionParser::Transcripts;

use base qw(BP::DCCLoader::Parsers::RNASeqStarInsertionParser);

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
		'RNA_TRANSCRIPT_QUANT_STAR_CRG'	=>	[
			'exp',
			['15b'],
			undef,
			[['.results' => 'tq_crg']],
			't',
			$self,
			BP::DCCLoader::Parsers::RNASeqStarInsertionParser::CRG_STAR_METADATA,
			undef
		],
	};
	
}

my $p_colnamesTranscripts = [
	'transcript_id',
	'gene_id',
	'length',
	'effective_length',
	'expected_count',
	'TPM',
	'FPKM',
	'IsoPct',
	'posterior_mean_count',
	'posterior_standard_deviation_of_count',
	'pme_TPM',
	'pme_FPKM',
	'IsoPct_from_pme_TPM',
	'TPM_ci_lower_bound',
	'TPM_ci_upper_bound',
	'FPKM_ci_lower_bound',
	'FPKM_ci_upper_bound',
];

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$) {
	my $self = shift;
	
	return $self->commonInsert($_[0],$_[1],$p_colnamesTranscripts,1);
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
