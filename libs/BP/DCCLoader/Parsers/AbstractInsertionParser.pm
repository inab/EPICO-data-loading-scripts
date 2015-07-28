#!/usr/bin/perl

use warnings;
use strict;

use Carp;

package BP::DCCLoader::Parsers::AbstractInsertionParser;

sub new() {
	# Very special case for multiple inheritance handling
	# This is the seed
	my($facet)=shift;
	my($class)=ref($facet) || $facet;
	
	my $ref = shift;
	my %href = ();
	
	%href = %{$ref}  if(ref($ref) eq 'HASH');
	
	return bless(\%href,$class);
}

# A template of metadata
use constant UNK_METADATA => {
	'assembly_version'	=>	1,
	'program_versions'	=>	[
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'',
		'url'	=>	'',
	},
	'other_analysis_algorithm'	=>	[
	],
};

# It returns a hash reference with the parsing features to be inserted in %FILETYPE2ANAL
sub getParsingFeatures() {
	my($self)=shift;
	
	Carp::carp("Empty method! You should implement it (as well the parser)!");
	
	return {};
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
	Carp::croak("Abstract method! You should implement it!");
}

# Uncomment this on the different parsers in order to self-register the singletons
# BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
