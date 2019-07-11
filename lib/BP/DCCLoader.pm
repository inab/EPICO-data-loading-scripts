#!/usr/bin/perl

use v5.12;
use strict;
use warnings 'all';

package BP::DCCLoader;

use version;

our $VERSION = version->declare('v1.0.0');

1;
__END__

=encoding utf8

=head1 NAME

EPICO data loading scripts

Then, you have to setup the project profile, either changing of copying L<epico-setup.ini|share/epico-setup.ini> file.

=over

=item *

Supported database paradigms (through _B<loaders>_ configuration parameter in B<mapper> section) are I<relational> (through _B<sql-dialect>_ configuration
  parameter from B<relational> section, subtypes I<postgresql>, I<mysql> or I<sqlite3>), I<mongodb> and I<elasticsearch> (these last ones have their own dedicated sections).
  On the dedicated section, uncomment and fill in only the parameters you need.

=item *

You also have to setup several connection parameters, like _B<host>_, _B<index-path>_, _B<user>_ and _B<pass>_ on B<dcc-loader> section in case you are using the BLUEPRINT data loading specific scripts.

=item *

The usage of these scripts is:

=over 4

C<perl epico-feature-loader.pl [-t] [-skipReactome] {iniFile} {cachindDir}>
	
C<perl blueprint-dcc-loader.pl [-t|-tt] [-s|-ss] {iniFile} {cachingDir} [sdata|dnase|meth|pdna|rnaseq]>

=back

=over

=item *
C<{inifile}> is the customized copy of L<epico-setup.ini|share/epico-setup.ini> you have previously prepared.
  
=item *
C<{cachingDir}> is a directory where the files to be processes are downloaded and cached. For primary data, those files are removed once they are processed.
  
=item *
Optional flags C<-t> and C<-tt> are used to put the scripts in test mode, so all the tasks are done but the insertion in the database.
  
=item *
Optional flags C<-s> and C<-ss> are to skip the insertion of primary data and its related specific analysis metadata.
  
=item *
Optional flag C<-skipReactome> is usually used in combination with C<-t> flag, in order to skip Reactome fetch and processing.
  
=back

=back

=head1 SEE ALSO

=cut