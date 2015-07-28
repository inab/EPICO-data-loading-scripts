#!/usr/bin/perl

use warnings;
use strict;

use BP::Loader::Tools;

use Carp;
use SQL::Parser ();

package BP::DCCLoader::Parsers::MySQLSchemaParser;

# constructor parameters:
#	file: The path to a gzip compressed MySQL SQL schema
# This method parses a mysql schema, returning the columns for each declared table
sub new($) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	my %tables = ();
	$self = bless(\%tables,$class);
	
	my $file = shift;

	my $content;
	if(open(my $SQL,'-|',BP::Loader::Tools::GUNZIP,'-c',$file)) {
		local $/;
		$content = <$SQL>;
		close($SQL);
	}

	my $sqlP = SQL::Parser->new('ANSI');


	my @sentences = split(/;[\n ]*/,$content);
	foreach my $sentence (@sentences) {
		# This removes the tail after the definition
		$sentence =~ s/\)[^\)]+$/)/;
		
		$sentence =~ s/`(module|value|absolute|sequence|timestamp|position|translation|end)`/"${1}"/mig;
		$sentence =~ s/`([^`_,\(\n]+)`/$1/mg;
		$sentence =~ s/`([^`,\(\n]+)`/"${1}"/mg;
		$sentence =~ s/ int(\([0-9]*\))/ integer${1}/mig;
		$sentence =~ s/tinyint\(/smallint(/mig;
		$sentence =~ s/smallinteger\(/smallint(/mig;
		$sentence =~ s/ longtext([ ,])/ varchar(4294967296)${1}/mig;
		$sentence =~ s/ text([ ,])/ varchar(2147483648)${1}/mig;
		$sentence =~ s/ mediumtext([ ,])/ varchar(65536)${1}/mig;
		$sentence =~ s/ tinytext([ ,])/ varchar(4096)${1}/mig;
		$sentence =~ s/ mediumblob([ ,])/ varchar(2147483648)${1}/mig;
		$sentence =~ s/ double([ ,])/ double precision${1}/mig;
		$sentence =~ s/ datetime([ ,])/ timestamp${1}/mig;
		$sentence =~ s/(?:enum|set)\([^)]+\)/char(256)/mig;
		$sentence =~ s/^( +)key /${1}unique key /mig;
		$sentence =~ s/ unique key [^(]+\(/ unique (/mig;
		
		# Removing stuff which hurts the parser
		$sentence =~ s/ unsigned / /mig;
		$sentence =~ s/ auto_increment */ /mig;
		$sentence =~ s/ default [^,\n]*//mig;
		$sentence =~ s/ CHARACTER SET [^ ]+ COLLATE [^ ,]+//mig;
		
		# And removing unique keys, which hurts it A LOT
		$sentence =~ s/ unique .*$//mig;
		
		$sqlP->parse($sentence);
		
		my $structure = $sqlP->structure;
		Carp::croak("INTERNAL ERROR: could not correctly parse next table:\n\n$sentence\n\nReason: ".$structure->{'errstr'})  if(exists($structure->{'errstr'}));
		#print STDERR Dumper($structure),"\n";
		my $tablename = $structure->{'org_table_names'}[0];
		$tablename =~ tr/"//d;
		
		my %colnames = ();
		my $colpos = 0;
		foreach my $colname (@{$structure->{'org_col_names'}}) {
			$colname =~ tr/"//d;
			$colnames{$colname} = $colpos;
			$colpos++;
		}
		
		$tables{$tablename} = \%colnames;
	}
	
	return $self;
}

# This method takes a table name and a list of columns,
# and it returns the column positions (needed to parse MySQL dumps)
sub mapTableColumns($@) {
	my $self = shift;
	
	Carp::croak((caller(0))[3].' is an instance method!')  unless(ref($self));
	
	my($tableName,@columnNames)=@_;
	
	my @colpos = ();
	
	if(exists($self->{$tableName})) {
		my $p_table = $self->{$tableName};
		foreach my $colname (@columnNames) {
			if(exists($p_table->{$colname})) {
				push(@colpos,$p_table->{$colname});
			} else {
				Carp::croak("ERROR: could not find mapping for column $colname in table $tableName");
			}
		}
	} else {
		Carp::croak("ERROR: could not find table $tableName in SQL definition");
	}
	
	return \@colpos;
}

1;
