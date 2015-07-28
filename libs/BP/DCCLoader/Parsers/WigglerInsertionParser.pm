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
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

sub getParsingFeatures() {
	my($self)=shift;
	
	my $CHIP_WIGGLER = [
		'pdna',
		['1'],
		'cs_wiggler',
		undef,
		undef,
		undef,
		WIGGLER_METADATA,
		undef
	];
	
	my $DNASE_WIGGLER = [
		'rreg',
		['8'],
		'ds_wiggler',
		undef,
		undef,
		undef,
		WIGGLER_METADATA,
		undef
	];
	
	return {
		# The two firsts are the old ones
		'CS_WIGGLER'	=>	$CHIP_WIGGLER,
		'DS_WIGGLER'	=>	$DNASE_WIGGLER,
		
		'CHIP_WIGGLER'	=>	$CHIP_WIGGLER,
		'DNASE_WIGGLER'	=>	$DNASE_WIGGLER,
	};
}

# This call registers the parser
BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
