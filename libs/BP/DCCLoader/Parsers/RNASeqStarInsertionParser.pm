#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use BP::DCCLoader::WorkingDir;
use TabParser;

package BP::DCCLoader::Parsers::RNASeqStarInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use BP::DCCLoader::Parsers::GencodeGTFParser;

use constant CRG_STAR_METADATA => {
	'assembly_version'	=>	38,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEMTools library',
		'url'	=>	'https://github.com/gemtools/gemtools',
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

{
my $p_GThash = undef;	

sub getGencodeCoordinates() {
	my $self = shift;
	
	unless($p_GThash) {
		$p_GThash = BP::DCCLoader::Parsers::GencodeGTFParser::getGencodeCoordinates($self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL},$self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_WORKINGDIR},$self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_INI});
	}
	
	return $p_GThash;
}

}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub commonInsert($$\@\@\@;$) {
	my $self = shift;
	
	my($analysis_id,$p_insertMethod,$p_colnames,$p_baseIndexes,$p_metricsIndexes,$isTranscript) = @_;
	
	# Get the Ensembl hash
	my $p_ensHash = $self->getGencodeCoordinates();
	
	my %rnaRSemStarGeneParserConfig = (
		TabParser::TAG_HAS_HEADER	=>	1,
		TabParser::TAG_FETCH_COLS	=>	$p_colnames,
		TabParser::TAG_CALLBACK => sub {
			my %metrics = ();
			
			# Normalizing the content
			@metrics{@{$p_colnames}[@{$p_metricsIndexes}]} = map { $_+0e0 } @_[@{$p_metricsIndexes}];
			
			my %entry = (
				'analysis_id'	=>	$analysis_id,
				'is_annotated'	=>	1,
				'metrics'	=>	\%metrics
			);
			
			my @values = @_[@{$p_baseIndexes}];
			foreach my $val (@values[1..$#values]) {
				# Normalizing the content
				$val += 0e0;
			}
			
			@entry{'gene_stable_id','length','effective_length','expected_count','TPM','FPKM'} = @values;
			
			# Choosing the ensId
			my $ensId;
			if(defined($isTranscript)) {
				$ensId = $entry{'transcript_stable_id'} = $_[$isTranscript];
			} else {
				$ensId = $entry{'gene_stable_id'};
			}
			$ensId = substr($ensId,0,rindex($ensId,'.'));
			
			# Now, let's get the chromosomical coordinates
			if(exists($p_ensHash->{$ensId})) {
				my $p_data = $p_ensHash->{$ensId};
				my $p_coordinates = $p_data->{'coordinates'}[0];
				$entry{'chromosome'} = $p_coordinates->{'chromosome'};
				$entry{'chromosome_start'} = $p_coordinates->{'chromosome_start'};
				$entry{'chromosome_end'} = $p_coordinates->{'chromosome_end'};
				
				$p_insertMethod->(\%entry);
			} else {
				print STDERR "DEBUGRNA: ",$ensId,"\n";
				Carp::carp("ERROR: Next entry was rejected as its Ensembl Id does not match with the ones defined\n".join("\t",@{$p_colnames})."\n".join("\t",@_));
			}
			
			1;
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
			[
				'exp',
				['15b'],
				undef,
				[['.results' => 'gq_crg']],
				'g',
				$self,
				BP::DCCLoader::Parsers::RNASeqStarInsertionParser::CRG_STAR_METADATA,
				undef
			]
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
	
	return $self->commonInsert($_[0],$_[1],$p_colnamesGenes,[0,2,3,4,5,6],[7..14]);
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
			[
				'exp',
				['15b'],
				undef,
				[['.results' => 'tq_crg']],
				't',
				$self,
				BP::DCCLoader::Parsers::RNASeqStarInsertionParser::CRG_STAR_METADATA,
				undef
			]
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
	
	return $self->commonInsert($_[0],$_[1],$p_colnamesTranscripts,[1,2,3,4,5,6],[7..16],0);
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
