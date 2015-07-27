#!/usr/bin/perl

use warnings;
use strict;

package BP::DCCLoader::Parsers::WigglerInsertionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant WIGGLER_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'',
		'url'	=>	'',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'WIGGLER',
			'url'	=>	'https://code.google.com/p/align2rawsignal/'
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
		'CHIP_WIGGLER'	=>	[
			'pdna',
			['1'],
			'cs_wiggler',
			undef,
			undef,
			undef,
			WIGGLER_METADATA,
			undef
		],
		'DNASE_WIGGLER'	=>	[
			'rreg',
			['8'],
			'ds_wiggler',
			undef,
			undef,
			undef,
			WIGGLER_METADATA,
			undef
		],
	};
}

1;
