#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers;
use TabParser;

package BP::DCCLoader::Parsers::DNASEBedInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant DS_METADATA => {
	'assembly_version'	=>	38,
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

# This is the empty constructor
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	my $DS_HOTSPOT_BED = [
		'rreg',
		['8'],
		'ds_hotspots',
		[['peaks' => 'ds_hotspots_peaks']],
		'p',
		$self,
		DS_METADATA,
		undef
	];
	
	return {
		# The first is the old one
		'DS_HOTSPOT'	=>	[ $DS_HOTSPOT_BED ],
		
		'DNASE_HOTSPOT_BED'	=>	[ $DS_HOTSPOT_BED ],
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
				'chromosome_end'	=>	$chromosome_end+0,	# Bed holds the end coordinate as exclusive, so it does not change
				'z_score'	=>	(defined($zscore_peak)?$zscore_peak:$zscore)+0.0,
			);
			
			$p_insertMethod->(\%entry);
		},
	);
	return \%dsHotspotsBedParserConfig;
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
