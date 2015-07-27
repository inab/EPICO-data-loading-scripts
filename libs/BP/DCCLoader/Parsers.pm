#!/usr/bin/perl

use warnings;
use strict;

use Carp;

package BP::DCCLoader::Parsers;

# These constants define the elements of the arrays used as values for
# the hash returned by getParsingFeatures
use constant {
	F_DOMAIN	=>	0,
	F_ANALYSIS_GROUP_ID	=>	1,
	F_POSTFIX	=>	2,
	F_PATTERN_POSTFIX	=>	3,
	F_PRIMARY	=>	4,
	F_PARSER	=>	5,
	F_METADATA	=>	6,
	F_PARENT_POSTFIX	=>	7,
};

{
my %FILETYPE2ANAL = ();

sub _registerParsableFiletypes(\%) {
	if(ref($_[0]) eq 'HASH') {
		my $p_parsable = $_[0];
		
		foreach my $newtype (keys(%{$p_parsable})) {
			if(exists($FILETYPE2ANAL{$newtype})) {
				Carp::carp("WARNING: Trying to register twice the filetype $newtype. Ignoring...");
			} else {
				$FILETYPE2ANAL{$newtype} = $p_parsable->{$newtype};
			}
		}
	}
}

# It returns the registered parsable filetypes
sub getParsableFiletypes() {
	return \%FILETYPE2ANAL;
}

}

1;
