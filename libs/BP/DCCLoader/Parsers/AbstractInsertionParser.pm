#!/usr/bin/perl

use warnings;
use strict;

use Carp;
use BP::DCCLoader::Parsers;

package BP::DCCLoader::Parsers::AbstractInsertionParser;

use TabParser;

# Next keywords are recognized:
use constant {
	K_INI	=>	'ini',
	K_MODEL	=>	'model',
	K_WORKINGDIR	=>	'workingDir',
	K_TESTMODE	=>	'testmode',
};

sub new(;\%) {
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
	'assembly_version'	=>	38,
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
	my($self)=shift;
	
	my($F,$analysis_id,$mapper) = @_;
	
	# UGLY
	my $BMAX = $mapper->bulkBatchSize();
	
	my $numBatch = 0;
	my @batch = ();
	
	my $p_insertMethod;
	if(exists($self->{K_TESTMODE()}) && $self->{K_TESTMODE()}) {
		$p_insertMethod = sub {
			eval {
				$mapper->validateAndEnactEntry($_[0]);
			};
			
			if($@) {
				print STDERR "DEBUG: WTF! $@\n";
			}
			1;
		};
	} else {
		$p_insertMethod = sub {
			push(@batch,$_[0]);
			$numBatch++;
			
			if($numBatch >= $BMAX) {
				eval {
					$mapper->bulkInsert(\@batch);
				};
				
				if($@) {
					print STDERR "DEBUG: WTF! $@\n";
				}
				
				@batch = ();
				$numBatch = 0;
			}
			1;
		};
	}
	
	my $p_tabConfig = $self->_insertInternal($analysis_id,$p_insertMethod);
	TabParser::parseTab($F,%{$p_tabConfig});
	
	# Last step
	if($numBatch > 0) {
		eval {
			if(exists($self->{K_TESTMODE()}) && $self->{K_TESTMODE()}) {
				$mapper->validateAndEnactEntry(\@batch);
			} else {
				$mapper->bulkInsert(\@batch);
			}
		};
		
		if($@) {
			print STDERR "DEBUG: WTF! $@\n";
		}
		@batch = ();
	}
}

#####
# Parser method bodies
# --------------
# Each method must take these parameters
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
# It returns a configuration hash usable with TabParser::parseTab
#####
sub _insertInternal($$) {
	Carp::croak("Abstract method! You should implement it!");
}

# Uncomment this on the different parsers in order to self-register the singletons
# BP::DCCLoader::Parsers::_registerParsableFiletypes(__PACKAGE__);

1;
