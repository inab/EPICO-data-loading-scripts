#!/usr/bin/perl

use warnings;
use strict;

use Carp;

package BP::DCCLoader::Parsers;

use constant DCC_LOADER_SECTION	=> 'dcc-loader';
use constant EPICO_LOADER_SECTION	=> 'epico-loader';
use constant ENSEMBL_FTP_BASE_TAG	=> 'ensembl-ftp-base-uri';

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
my @PARSERS = ();

sub _registerParsableFiletypes(@) {
	push(@PARSERS,@_);
}

my %FILETYPE2ANAL = ();

# It returns the registered parsable filetypes
sub getParsableFiletypes(@) {
	my $class = shift;
	
	if(scalar(@PARSERS) > 0) {
		foreach my $clazz (@PARSERS) {
			my $instance = $clazz->new(@_);
			my $p_parsable = $instance->getParsingFeatures();
			
			foreach my $newtype (keys(%{$p_parsable})) {
				if(exists($FILETYPE2ANAL{$newtype})) {
					Carp::carp("WARNING: Trying to register twice the filetype $newtype while loading $clazz. Ignoring...");
				} else {
					$FILETYPE2ANAL{$newtype} = $p_parsable->{$newtype};
				}
			}
		}
		
		@PARSERS = ();
	}
	
	return \%FILETYPE2ANAL;
}

}

1;
