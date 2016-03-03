#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers::MySQLSchemaParser;
use BP::DCCLoader::WorkingDir;

use BP::Loader::Tools;

use Carp;

use Net::FTP::AutoReconnect;

use TabParser;

use URI;

package BP::DCCLoader::Parsers::EnsemblGTParser;

use Log::Log4perl;

my $LOG;
sub BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::WARN, layout => '%-5p - %m%n' } );

	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
}

use constant {
	ENSEMBL_SQL_FILE	=> 'homo_sapiens_core_{EnsemblVer}_{GRChVer}.sql.gz',
	ENSEMBL_SEQ_REGION_FILE	=> 'seq_region.txt.gz',
	ENSEMBL_GENE_FILE	=> 'gene.txt.gz',
	ENSEMBL_TRANSCRIPT_FILE	=> 'transcript.txt.gz',
	ENSEMBL_EXON_TRANSCRIPT_FILE	=> 'exon_transcript.txt.gz',
	ENSEMBL_EXON_FILE	=> 'exon.txt.gz',
	ENSEMBL_XREF_FILE	=> 'xref.txt.gz',
	ENSEMBL_EXTERNAL_DB_FILE	=> 'external_db.txt.gz',
};

#####
# Method bodies
#####
sub parseSeqRegions($$$) {
	my($payload,$seq_region_id,$name) = @_;
	
	my($p_regionId,$chroCV) = @{$payload};
	
	my $term = $chroCV->getTerm($name);
	$p_regionId->{$seq_region_id} = $term->key()  if($term);
	
	1;
}

use constant ENSEMBL_NS => 'Ensembl';
use constant DESCRIPTION_NS => 'description';

sub parseENS($$$$$$$$) {
	my($payload,$seq_region_id,$chromosome_start,$chromosome_end,$chromosome_strand,$stable_id,$version,$internal_id,$display_xref_id,$description,$internal_gene_id)=@_;
	
	my($p_regionId,$p_ENShash,$p_ENSintHash,$p_geneMap,$chroCV,$p_exonsMap) = @{$payload};
	
	if(exists($p_regionId->{$seq_region_id})) {
		my $fullStableId = $stable_id.'.'.$version;
		my @symbols = (
			{
				'domain'	=>	ENSEMBL_NS,
				'value'	=>	[$fullStableId,$stable_id],
			}
		);
		$description = undef  if(defined($description) && $description eq "\\N");
		
		my $parsedData = undef;
		if(exists($p_ENShash->{$stable_id})) {
			$parsedData = $p_ENShash->{$stable_id};
			push(@{$parsedData->{'symbol'}},@symbols);
		} else {
			my $p_feature_cluster_ids;
			if(defined($internal_gene_id)) {
				$p_feature_cluster_ids = [ @{$p_geneMap->{$internal_gene_id}->{'feature_cluster_id'}}, $fullStableId ];
			} elsif(defined($display_xref_id)) {
				$p_feature_cluster_ids = [ $fullStableId ];
			} else {
				if(exists($p_exonsMap->{$internal_id})) {
					$p_feature_cluster_ids = $p_exonsMap->{$internal_id};
					push(@{$p_feature_cluster_ids},$fullStableId);
				} else {
					# This case should not happen
					$p_feature_cluster_ids = [ $fullStableId ];
				}
			}
			
			$parsedData = {
				#$fullStableId,
				'feature_cluster_id'	=> $p_feature_cluster_ids,
				'feature_ns'	=> ENSEMBL_NS,
				'feature_id'	=> $fullStableId,
				'feature'	=> (defined($internal_gene_id)?'transcript':(defined($display_xref_id) ? 'gene' : 'exon')),
				'coordinates'	=> [
					{
						'feature_ns'	=> ENSEMBL_NS,
						'feature_id'	=> $fullStableId,
						'chromosome'	=> $p_regionId->{$seq_region_id},
						'chromosome_start'	=> ($chromosome_start+0),
						'chromosome_end'	=> ($chromosome_end+0),
						'chromosome_strand'	=> ($chromosome_strand+0),
					}
				],
				'symbol'	=> \@symbols,
			};
			$p_ENShash->{$stable_id} = $parsedData;
			$p_geneMap->{$internal_id} = $parsedData  unless(defined($internal_gene_id) || !defined($display_xref_id));
			$p_exonsMap->{$internal_id} = $parsedData  if(defined($internal_gene_id));
		}
		if(defined($display_xref_id)) {
			push(@{$parsedData->{'symbol'}},{'domain' => DESCRIPTION_NS,'value' => [$description]})  if(defined($description));
			
			$p_ENSintHash->{$display_xref_id} = $parsedData;
		}
	}
	
	1;
}

sub parseET($$$) {
	my($payload,$internal_exon_id,$internal_transcript_id)=@_;
	
	my($p_transcriptsMap,$p_exonsMap) = @{$payload};
	
	if(exists($p_transcriptsMap->{$internal_transcript_id})) {
		if(exists($p_exonsMap->{$internal_exon_id})) {
			foreach my $elem (@{$p_transcriptsMap->{$internal_transcript_id}{'feature_cluster_id'}}) {
				$p_exonsMap->{$internal_exon_id}{$elem} = undef;
			}
		} else {
			my %set = map { $_ => undef } @{$p_transcriptsMap->{$internal_transcript_id}{'feature_cluster_id'}};
			$p_exonsMap->{$internal_exon_id} = \%set;
		}
	#} else {
	#	# This case happens with exons from LRG, so skip them
	#	print $internal_exon_id, ' JARL ', $internal_transcript_id,"\n";
	}
	
	1;
}

sub parseExternalDB($$$) {
	my($p_externalDB,$external_db_id,$db_name)=@_;
	
	$p_externalDB->{$external_db_id} = $db_name;
	
	1;
}

sub parseXREF($$$$$$) {
	my($payload,$display_xref_id,$external_db_id,$dbprimary_acc,$display_label,$description)=@_;
	
	my($p_ENSintHash,$p_externalDB) = @{$payload};
	
	if(exists($p_ENSintHash->{$display_xref_id})) {
		my $ns = exists($p_externalDB->{$external_db_id}) ? $p_externalDB->{$external_db_id} : $external_db_id;
		
		my $p_desc;
		if(exists($p_ENSintHash->{$display_xref_id}->{'attribute'})) {
			$p_desc = $p_ENSintHash->{$display_xref_id}->{'attribute'};
		} else {
			$p_desc = [];
			$p_ENSintHash->{$display_xref_id}->{'attribute'} = $p_desc;
		}
		
		my @values = ($dbprimary_acc);
		push(@values,$display_label)  if($dbprimary_acc ne $display_label);
		push(@values,$description)  unless($description eq "\\N");
		push(@{$p_desc},{
			'domain' => $ns,
			'value' => \@values,
		});
	}
	
	1;
}

# parseEnsemblGenesAndTranscripts parameters:
#	chroCV: A BP::Model::CV instance. These are the known chromosomes
#	testmode:
#	localSQLFile: path to the fetched Ensembl MySQL schema
#	localSeqRegion: path to the Ensembl MySQL dump of the region types
#	localGenes: path to the Ensembl MySQL dump of the known genes
#	localTranscripts: path to the Ensembl MySQL dump of the known transcripts
#	localExonTranscripts: path to the Ensembl MySQL dump of the correspondence between transcripts and exons
#	localExons: path to the Ensembl MySQL dump of the known exons
#	localXref: path to the Ensembl MySQL dump of the known cross-references
#	localExternalDB: path to the Ensembl MySQL dump of the external databases used for cross-references
# The returned hash reference holds the data to be dumped later
# The values will be array instances with the structure of
#	ensembl id (with version number)
#	chromosome (BLUEPRINT style)
#	chromosome_start (1-based)
#	chromosome_end (1-based)
#	symbols (an array)
#	feature
sub parseEnsemblGenesAndTranscripts($$$$$$$$$$) {
	my($chroCV,$testmode,$localSQLFile,$localSeqRegion,$localGenes,$localTranscripts,$localExonTranscripts,$localExons,$localXref,$localExternalDB) = @_;
	
	# These are the Ensembl columns for its tables
	my $p_ensTables = BP::DCCLoader::Parsers::MySQLSchemaParser->new($localSQLFile);
	
	# These hashes will hold the data to be dumped later
	# The values will be array instances with the structure of
	#	ensembl id (with version number)
	#	chromosome (BLUEPRINT style)
	#	chromosome_start (1-based)
	#	chromosome_end (1-based)
	#	symbols (an array)
	#	feature
	my %ENShash = ();
	my %ENSintHash = ();
	
	# This hash holds the seq_region_id -> name correspondence
	my %regionId = ();
	
	$LOG->info("Parsing ".$localSeqRegion);
	if(open(my $SEQREG,'-|',BP::Loader::Tools::GUNZIP,'-c',$localSeqRegion)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,$chroCV],
			TabParser::TAG_CALLBACK => \&parseSeqRegions,
			# The columns are 
			#	seq_region_id
			#	name
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('seq_region',
				'seq_region_id',
				'name'
			), # [0,1],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($SEQREG,%config);
		
		close($SEQREG);
	} else {
		Carp::croak("ERROR: Unable to parse Sequence Regions file ".$localSeqRegion);
	}
	
	my %geneMap;
	
	$LOG->info("Parsing ".$localGenes);
	if(open(my $ENSG,'-|',BP::Loader::Tools::GUNZIP,'-c',$localGenes)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,\%ENShash,\%ENSintHash,\%geneMap,$chroCV],
			TabParser::TAG_CALLBACK => \&parseENS,
			# The columns are 
			#	seq_region_id
			#	chromosome_start
			#	chromosome_end
			#	chromosome_strand
			#	display_xref_id
			#	ensembl gene id
			#	ensembl gene id version
			#	description
			#	internal gene id
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('gene',
				'seq_region_id',
				'seq_region_start',
				'seq_region_end',
				'seq_region_strand',
				'stable_id',
				'version',
				'gene_id',
				'display_xref_id',
				'description',
			), # [3,4,5,7,14,15,10,0],
#			TabParser::TAG_NEG_FILTER => [[1 => 'LRG_gene']],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($ENSG,%config);
		
		close($ENSG);
	} else {
		Carp::croak("ERROR: Unable to parse EnsEMBL Genes file ".$localGenes);
	}
	
	$LOG->info("Parsing ".$localTranscripts);
	my %transcriptsMap;
	if(open(my $ENST,'-|',BP::Loader::Tools::GUNZIP,'-c',$localTranscripts)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,\%ENShash,\%ENSintHash,\%geneMap,$chroCV,\%transcriptsMap],
			TabParser::TAG_CALLBACK => \&parseENS,
			# The columns are 
			#	seq_region_id
			#	chromosome_start
			#	chromosome_end
			#	chromosome_strand
			#	display_xref_id
			#	ensembl transcript id
			#	ensembl transcript id version
			#	description
			#	internal transcript id
			#	internal gene id
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('transcript',
				'seq_region_id',
				'seq_region_start',
				'seq_region_end',
				'seq_region_strand',
				'stable_id',
				'version',
				'transcript_id',
				'display_xref_id',
				'description',
				'gene_id'
			), # [3,4,5,7,13,14,10,0,1],
#			TabParser::TAG_NEG_FILTER => [[8 => 'LRG_gene']],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($ENST,%config);
		
		close($ENST);
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL Transcripts file ".$localTranscripts);
	}
	

	# Mapping from exons to transcripts
	my %exonsMap;
	
	$LOG->info("Parsing ".$localExonTranscripts);
	if(open(my $ET,'-|',BP::Loader::Tools::GUNZIP,'-c',$localExonTranscripts)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%transcriptsMap,\%exonsMap],
			TabParser::TAG_CALLBACK => \&parseET,
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('exon_transcript',
				'exon_id',
				'transcript_id',
			),
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($ET,%config);
		
		close($ET);
		
		# Post-processing
		foreach my $set (values(%exonsMap)) {
			$set = [ keys(%{$set}) ];
		}
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL Exons file ".$localExonTranscripts);
	}
	
	%transcriptsMap = ();
	
	$LOG->info("Parsing ".$localExons);
	if(open(my $ENSE,'-|',BP::Loader::Tools::GUNZIP,'-c',$localExons)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,\%ENShash,\%ENSintHash,\%geneMap,$chroCV,\%exonsMap],
			TabParser::TAG_CALLBACK => \&parseENS,
			# The columns are 
			#	seq_region_id
			#	chromosome_start
			#	chromosome_end
			#	chromosome_strand
			#	display_xref_id
			#	ensembl transcript id
			#	ensembl transcript id version
			#	description
			#	internal transcript id
			#	internal gene id
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('exon',
				'seq_region_id',
				'seq_region_start',
				'seq_region_end',
				'seq_region_strand',
				'stable_id',
				'version',
				'exon_id',
			), # [3,4,5,7,13,14,10,0,1],
#			TabParser::TAG_NEG_FILTER => [[8 => 'LRG_gene']],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($ENSE,%config);
		
		close($ENSE);
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL Exons file ".$localExons);
	}
	
	# Freeing unused memory
	undef %regionId;
	undef %geneMap;
	undef %exonsMap;
	
	$LOG->info("Parsing ".$localExternalDB);
	my %externalDB = ();
	if(open(my $ExtDB,'-|',BP::Loader::Tools::GUNZIP,'-c',$localExternalDB)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> \%externalDB,
			TabParser::TAG_CALLBACK => \&parseExternalDB,
			TabParser::TAG_MULTILINE_SEP => "\\",
			# The columns are 
			#	display_xref_id
			#	dbprimary_acc
			#	display_label
			#	description
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('external_db',
				'external_db_id',
				'db_name',
			), # [0,2,3,5],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($ExtDB,%config);
		
		close($ExtDB);
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL external DB file ".$localExternalDB);
	}
	
	$LOG->info("Parsing ".$localXref);
	if(open(my $XREF,'-|',BP::Loader::Tools::GUNZIP,'-c',$localXref)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%ENSintHash,\%externalDB],
			TabParser::TAG_CALLBACK => \&parseXREF,
			TabParser::TAG_MULTILINE_SEP => "\\",
			# The columns are 
			#	display_xref_id
			#	dbprimary_acc
			#	display_label
			#	description
			TabParser::TAG_FETCH_COLS => $p_ensTables->mapTableColumns('xref',
				'xref_id',
				'external_db_id',
				'dbprimary_acc',
				'display_label',
				'description'
			), # [0,2,3,5],
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($XREF,%config);
		
		close($XREF);
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL XREF file ".$localXref);
	}
	# Freeing unused memory
	undef %ENSintHash;
	
	return \%ENShash;
}

sub getEnsemblCoordinates($$$;$) {
	my($model,$workingDir,$ini,$testmode) = @_;
	
	my $ensembl_ftp_base;
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,BP::DCCLoader::Parsers::ENSEMBL_FTP_BASE_TAG)) {
		$ensembl_ftp_base = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,BP::DCCLoader::Parsers::ENSEMBL_FTP_BASE_TAG);
	} else {
		Carp::croak("Configuration file must have '".BP::DCCLoader::Parsers::ENSEMBL_FTP_BASE_TAG."' in '".BP::DCCLoader::Parsers::DCC_LOADER_SECTION."' section");
	}
	
	# Now, let's patch the properies of the different remote resources, using the properties inside the model
	my $ensembl_sql_file = BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_SQL_FILE;
	eval {
		$model->annotations->applyAnnotations(\($ensembl_ftp_base,$ensembl_sql_file));
	};
	
	if($@) {
		Carp::croak("$@ (does not exist in model)");
	}
	
	# And translate these to URI objects
	$ensembl_ftp_base = URI->new($ensembl_ftp_base);
	
	# Defined outside
	my $ftpServer = undef;
	
	# Fetching FTP resources
	$LOG->info("Connecting to $ensembl_ftp_base...");
	
	my $ensemblHost = $ensembl_ftp_base->host();
	$ftpServer = Net::FTP::AutoReconnect->new($ensemblHost,Debug=>0) || Carp::croak("FTP connection to server ".$ensemblHost." failed: ".$@);
	$ftpServer->login(BP::DCCLoader::WorkingDir::ANONYMOUS_USER,BP::DCCLoader::WorkingDir::ANONYMOUS_PASS) || Carp::croak("FTP login to server $ensemblHost failed: ".$ftpServer->message());
	$ftpServer->binary();
	
	my $ensemblPath = $ensembl_ftp_base->path;
	
	my $localSQLFile = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.$ensembl_sql_file);
	my $localSeqRegion = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_SEQ_REGION_FILE);
	my $localGenes = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_GENE_FILE);
	my $localTranscripts = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_TRANSCRIPT_FILE);
	my $localExonTranscripts = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_EXON_TRANSCRIPT_FILE);
	my $localExons = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_EXON_FILE);
	my $localXref = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_XREF_FILE);
	my $localExternalDB = $workingDir->cachedGet($ftpServer,$ensemblPath.'/'.BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_EXTERNAL_DB_FILE);
	
	Carp::croak("FATAL ERROR: Unable to fetch files from $ensemblPath (host $ensemblHost)")  unless(defined($localSQLFile) && defined($localSeqRegion) && defined($localGenes) && defined($localTranscripts) && defined($localXref));
	
	$ftpServer->disconnect()  if($ftpServer->can('disconnect'));
	$ftpServer->quit()  if($ftpServer->can('quit'));
	$ftpServer = undef;
	
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
		
	return parseEnsemblGenesAndTranscripts($chroCV,$testmode,$localSQLFile,$localSeqRegion,$localGenes,$localTranscripts,$localExonTranscripts,$localExons,$localXref,$localExternalDB);
}

1;
