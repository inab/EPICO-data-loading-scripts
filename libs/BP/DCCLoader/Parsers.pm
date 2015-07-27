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

1;
