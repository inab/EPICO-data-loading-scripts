#!/usr/bin/perl

use warnings;
use strict;

use boolean qw();

use BP::DCCLoader::Parsers;
use BP::DCCLoader::WorkingDir;
use TabParser;

package BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant CRG_COSI_METADATA => {
	'assembly_version'	=>	38,
	'program_versions'	=>	[
		{
			'program'	=>	'IPSA',
			'version'	=>	'v3.3'
		}
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'GEMTools library',
		'url'	=>	'https://github.com/gemtools/gemtools',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'IPSA',
			'url'	=>	'https://github.com/pervouchine/ipsa-full',
		},
	],
};

# This is the empty constructor
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub commonInsert($$$) {
	my $self = shift;
	
	my($analysis_id,$p_insertMethod,$p_prepareMethod) = @_;
	
	my %rnaCoSIParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				undef, # source
				$feature, # feature
				$chromosome_start,
				$chromosome_end,
				$score,
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
				'chromosome_start'	=>	$chromosome_start+0,
				'chromosome_end'	=>	$chromosome_end+0,
				'chromosome_strand'	=>	(($chromosome_strand eq '-')?-1:1),
			);
			
			$p_insertMethod->($p_prepareMethod->($feature,\%entry,\%attributes,$score));
			
			1;
		},
	);
	return \%rnaCoSIParserConfig;
}


package BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser::SplicingRatios;

use base qw(BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser);

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
		'RNA_COSI_STAR_CRG'	=>	[
			[
				'cosi',
				['15b'],
				undef,
				[['.splicing_ratios.' => 'cosi']],
				'spr',
				$self,
				BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser::CRG_COSI_METADATA,
				undef
			]
		],
	};
	
}

sub _prepareSplicingRatiosInternal($\%\%$) {
	my($feature,$p_entry,$p_attributes,$score) = @_;
	
	my $cosi;
	
	if(exists($p_attributes->{'cosi'})) {
		$cosi = $p_attributes->{'cosi'} + 0e0  unless($p_attributes->{'cosi'} eq 'NA');
	} else {
		my $cosi3 = ($p_attributes->{'cosi3'} ne 'NA') ? $p_attributes->{'cosi3'} + 0e0 : undef;
		my $cosi5 = ($p_attributes->{'cosi5'} ne 'NA') ? $p_attributes->{'cosi5'} + 0e0 : undef;
		
		my $doCalc = defined($cosi3) || defined($cosi5);
		$cosi5 = 0e0  unless(defined($cosi5));
		$cosi3 = 0e0  unless(defined($cosi3));
		
		$cosi = ($cosi3 + $cosi5) / 2.0  if($doCalc);
	}
	
	my $psi;
	
	if(exists($p_attributes->{'psi'})) {
		$psi = $p_attributes->{'psi'} + 0e0  unless($p_attributes->{'psi'} eq 'NA');
	} else {
		my $psi3 = ($p_attributes->{'psi3'} ne 'NA') ? $p_attributes->{'psi3'} + 0e0 : undef;
		my $psi5 = ($p_attributes->{'psi5'} ne 'NA') ? $p_attributes->{'psi5'} + 0e0 : undef;
		
		my $doCalc = defined($psi3) || defined($psi5);
		$psi5 = 0e0  unless(defined($psi5));
		$psi3 = 0e0  unless(defined($psi3));
		
		$psi = ($psi3 + $psi5) / 2.0  if($doCalc);
	}
	
	foreach my $metric_key (keys(%{$p_attributes})) {
		if($p_attributes->{$metric_key} eq 'NA') {
			# Not a Number
			$p_attributes->{$metric_key} = 'NaN' + 0;
		} else {
			$p_attributes->{$metric_key} += 0e0;
		}
	}
	
	$p_entry->{'cosi'} = $cosi;
	$p_entry->{'psi'} = $psi;
	$p_entry->{'is_exon'} = ($feature eq 'exon') ? boolean::true : boolean::false ;
	$p_entry->{'metrics'} = $p_attributes;
	
	return $p_entry;
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$) {
	my $self = shift;
	
	return $self->commonInsert($_[0],$_[1],\&_prepareSplicingRatiosInternal);
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);


package BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser::IPSAJunctions;

use base qw(BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser);

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
		'RNA_COSI_STAR_CRG'	=>	[
			[
				'cosi',
				['15b'],
				undef,
				[['.ipsa_junctions.' => 'cosi']],
				'jnc',
				$self,
				BP::DCCLoader::Parsers::RNASeqCoSIInsertionParser::CRG_COSI_METADATA,
				undef
			]
		],
	};
	
}

sub _prepareJunctionsInternal($\%\%$) {
	my($feature,$p_entry,$p_attributes,$score) = @_;
	
	$p_entry->{'idr'} = $p_attributes->{'IDR'} + 0e0;
	$p_entry->{'total_reads'} = $p_attributes->{'count'} + 0;
	$p_entry->{'stagg_reads'} = $p_attributes->{'stagg'} + 0;
	$p_entry->{'entropy'} = $p_attributes->{'entr'} + 0e0;
	$p_entry->{'annot_status'} = $p_attributes->{'annot'} + 0;
	$p_entry->{'splice_site_nucleotides'} = $p_attributes->{'nucl'};
	$p_entry->{'score'} = ($score + 0e0)/ 100.0;
	
	return $p_entry;
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#####
sub _insertInternal($$) {
	my $self = shift;
	
	return $self->commonInsert($_[0],$_[1],\&_prepareJunctionsInternal);
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
