#!/usr/bin/perl -w

use strict;

#use diagnostics;
use FindBin;
use lib $FindBin::Bin."/model/schema+tools/lib";

use Carp;
use File::Basename;
use File::Path;
use File::Spec;
use IO::Handle;
use Net::FTP::AutoReconnect;
use URI;

use BP::Model;
use BP::Loader::CorrelatableConcept;
use BP::Loader::Mapper;
use BP::Loader::Mapper::Relational;
use BP::Loader::Mapper::Elasticsearch;
use BP::Loader::Mapper::MongoDB;

use TabParser;

use constant DCC_LOADER_SECTION => 'dcc-loader';

use constant {
	ENSEMBL_FTP_BASE_TAG	=> 'ensembl-ftp-base-uri',
	GENCODE_FTP_BASE_TAG	=> 'gencode-ftp-base-uri',
	GENCODE_GTF_FILE_TAG	=> 'gencode-gtf',
	ENSEMBL_SEQ_REGION_FILE	=> 'seq_region.txt.gz',
	ENSEMBL_GENE_FILE	=> 'gene.txt.gz',
	ENSEMBL_TRANSCRIPT_FILE	=> 'transcript.txt.gz',
	ENSEMBL_XREF_FILE	=> 'xref.txt.gz',
	ANONYMOUS_USER	=> 'ftp',
	ANONYMOUS_PASS	=> 'guest@',
};

sub cachedGet($$$);

sub parseSeqRegions($$$);
sub parseENS($$$$$$$$);
sub parseXREF($$$$$);

#####
# Method bodies
#####
sub cachedGet($$$) {
	my($ftpServer,$remotePath,$cachingDir)=@_;
	
	my $filedate = $ftpServer->mdtm($remotePath);
	my $filesize = $ftpServer->size($remotePath);
	
	my $localPath = File::Spec->catfile($cachingDir,$remotePath);
	my $localBasePath = File::Basename::basename($remotePath);
	my $localRelDir = File::Basename::dirname($remotePath);
	my $localDir = File::Spec->catdir($cachingDir,$localRelDir);
	
	my $mirrored = undef;
	if(-f $localPath) {
		my($localsize,$localdate) = ( stat($localPath) )[7,9];
		$mirrored = $filedate == $localdate && $filesize == $localsize;
	}
	
	unless($mirrored) {
		$remotePath = '/'.$remotePath  unless(substr($remotePath,0,1) eq '/');
		File::Path::make_path($localDir);
		#print STDERR join(" -=- ",$remotePath,$cachingDir,$localPath,$localBasePath,$localRelDir,$localDir),"\n";
		$localPath = $ftpServer->get($remotePath,$localPath);
		print STDERR "DEBUGFTP: ".$ftpServer->message."\n"  unless(defined($localPath));
		utime($filedate,$filedate,$localPath)  if(defined($localPath));
	}
	
	return $localPath;
}

sub parseSeqRegions($$$) {
	my($p_regionId,$seq_region_id,$name) = @_;
	
	$p_regionId->{$seq_region_id} = $name;
	
	1;
}

sub parseENS($$$$$$$$) {
	my($payload,$seq_region_id,$chromosome_start,$chromosome_end,$display_xref_id,$stable_id,$version,$description)=@_;
	
	my($p_regionId,$p_ENShash,$p_ENSintHash,$feature) = @{$payload};
	
	if(exists($p_regionId->{$seq_region_id})) {
		my $fullStableId = $stable_id.'.'.$version;
		my @symbols = ($stable_id,$fullStableId);
		push(@symbols,$description)  unless($description eq "\\N");
		
		my $parsedData = undef;
		if(exists($p_ENShash->{$stable_id})) {
			$parsedData = $p_ENShash->{$stable_id};
			push(@{$parsedData->{'symbol'}},@symbols);
		} else {
			$parsedData = {
				#$fullStableId,
				'chromosome'	=> $p_regionId->{$seq_region_id},
				'chromosome_start'	=> $chromosome_start,
				'chromosome_end'	=> $chromosome_end,
				'symbol'	=> \@symbols,
				'feature'	=> $feature
			};
			$p_ENShash->{$stable_id} = $parsedData;
		}
		
		$p_ENSintHash->{$display_xref_id} = $parsedData;
	}
	
	1;
}

sub parseXREF($$$$$) {
	my($p_ENSintHash,$display_xref_id,$dbprimary_acc,$display_label,$description)=@_;
	
	if(exists($p_ENSintHash->{$display_xref_id})) {
		my $p_desc = $p_ENSintHash->{$display_xref_id}->{'symbol'};
		push(@{$p_desc},$dbprimary_acc);
		
		push(@{$p_desc},$display_label)  if($dbprimary_acc ne $display_label);
		push(@{$p_desc},$description)  unless($description eq "\\N");
	}
	
	1;
}

sub parseGTF(@) {
	my(
		$payload,
		$chro,
		undef, # source
		$feature, # feature
		$chromosome_start,
		$chromosome_end,
		undef,
		$chromosome_strand,
		undef, # frame
		$attributes_str, # attributes following .ace format
	) = @_;
	
	my($p_ENShash,$p_GencodeObjects) = @{$payload};
	
	my %attributes = ();
	
	#print STDERR "DEBUG: ",join(" ",@_),"\n"  unless(defined($attributes_str));
	my @tokens = split(/\s*;\s*/,$attributes_str);
	foreach my $token (@tokens) {
		my($key,$value) = split(/\s+/,$token,2);
		
		# Removing double quotes
		$value =~ tr/"//d;
		
		$attributes{$key} = $value;
	}
	
	my $p_regionData = undef;
	
	my $ensIdKey = undef;
	my $p_ensFeatures = undef;
	if($feature eq 'gene') {
		$ensIdKey = 'gene_id';
		$p_ensFeatures = ['gene_name','havana_gene'];
	} else {
		$ensIdKey = 'transcript_id';
		$p_ensFeatures = ['transcript_name','havana_transcript','ccdsid'];
	}
	
	my $fullEnsemblId = $attributes{$ensIdKey};
	my $ensemblId = substr($fullEnsemblId,0,index($fullEnsemblId,'.'));
		
	unless($feature eq 'gene' || $feature eq 'transcript') {
		$p_regionData = $p_ENShash->{$ensemblId}  if(exists($p_ENShash->{$ensemblId}));
	}
	
	unless(defined($p_regionData)) {
		my $chromosome = (index($chro,'chr')==0)?substr($chro,3):$chro;
		
		$chromosome = 'MT'  if($chromosome eq 'M');
		
		$p_regionData = {
			#$fullEnsemblId,
			'chromosome'	=> $chromosome,
			'chromosome_start'	=> $chromosome_start,
			'chromosome_end'	=> $chromosome_end,
			'symbol'	=> [$fullEnsemblId,$ensemblId],
			'feature'	=> $feature
		};
		push(@{$p_GencodeObjects},$p_regionData);
	}
	
	foreach my $ensFeature (@{$p_ensFeatures}) {
		push(@{$p_regionData->{symbol}},$attributes{$ensFeature})  if(exists($attributes{$ensFeature}));
	}
	
	1;
}

	

#####
# main
#####
if(scalar(@ARGV)>=2) {
	STDOUT->autoflush(1);
	STDERR->autoflush(1);
	my $iniFile = shift(@ARGV);
	# Defined outside
	my $cachingDir = shift(@ARGV);
	
	# First, let's read the configuration
	my $ini = Config::IniFiles->new(-file => $iniFile, -default => $BP::Loader::Mapper::DEFAULTSECTION);
	
	my $ensembl_ftp_base = undef;
	my $gencode_ftp_base = undef;
	my $gencode_gtf_file = undef;
	my $ensembl_gene_file = undef;
	my $ensembl_transcript_file = undef;
	my $ensembl_xref_file = undef;
	
	# Check the needed parameters for the construction
	if($ini->exists(DCC_LOADER_SECTION,ENSEMBL_FTP_BASE_TAG)) {
		$ensembl_ftp_base = URI->new($ini->val(DCC_LOADER_SECTION,ENSEMBL_FTP_BASE_TAG));
	} else {
		Carp::croak("Configuration file $iniFile must have '".ENSEMBL_FTP_BASE_TAG."'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,GENCODE_FTP_BASE_TAG)) {
		$gencode_ftp_base = URI->new($ini->val(DCC_LOADER_SECTION,GENCODE_FTP_BASE_TAG));
	} else {
		Carp::croak("Configuration file $iniFile must have '".GENCODE_FTP_BASE_TAG."'");
	}
	
	if($ini->exists(DCC_LOADER_SECTION,GENCODE_GTF_FILE_TAG)) {
		$gencode_gtf_file = $ini->val(DCC_LOADER_SECTION,GENCODE_GTF_FILE_TAG);
	} else {
		Carp::croak("Configuration file $iniFile must have '".GENCODE_GTF_FILE_TAG."'");
	}
	
	# First, explicitly create the caching directory
	File::Path::make_path($cachingDir);
	
	print "Connecting to $ensembl_ftp_base...\n";
	# Defined outside
	my $ftpServer = undef;
	
	my $ensemblHost = $ensembl_ftp_base->host();
	$ftpServer = Net::FTP::AutoReconnect->new($ensemblHost,Debug=>0) || Carp::croak("FTP connection to server ".$ensemblHost." failed: ".$@);
	$ftpServer->login(ANONYMOUS_USER,ANONYMOUS_PASS) || Carp::croak("FTP login to server $ensemblHost failed: ".$ftpServer->message());
	$ftpServer->binary();
	
	my $ensemblPath = $ensembl_ftp_base->path;
	
	my $localSeqRegion = cachedGet($ftpServer,$ensemblPath.'/'.ENSEMBL_SEQ_REGION_FILE,$cachingDir);
	my $localGenes = cachedGet($ftpServer,$ensemblPath.'/'.ENSEMBL_GENE_FILE,$cachingDir);
	my $localTranscripts = cachedGet($ftpServer,$ensemblPath.'/'.ENSEMBL_TRANSCRIPT_FILE,$cachingDir);
	my $localXref = cachedGet($ftpServer,$ensemblPath.'/'.ENSEMBL_XREF_FILE,$cachingDir);
	
	Carp::croak("FATAL ERROR: Unable to fetch files from $ensemblPath (host $ensemblHost)")  unless(defined($localSeqRegion) && defined($localGenes) && defined($localTranscripts) && defined($localXref));
	
	$ftpServer->disconnect()  if($ftpServer->can('disconnect'));
	$ftpServer->quit()  if($ftpServer->can('quit'));
	$ftpServer = undef;
	
	print "Connecting to $gencode_ftp_base...\n";
	my $gencodeHost = $gencode_ftp_base->host();
	$ftpServer = Net::FTP::AutoReconnect->new($gencodeHost,Debug=>0) || Carp::croak("FTP connection to server ".$gencodeHost." failed: ".$@);
	$ftpServer->login(ANONYMOUS_USER,ANONYMOUS_PASS) || Carp::croak("FTP login to server $gencodeHost failed: ".$ftpServer->message());
	$ftpServer->binary();
	
	my $gencodePath = $gencode_ftp_base->path();
	
	my $localGTF = cachedGet($ftpServer,$gencodePath.'/'.$gencode_gtf_file,$cachingDir);
	
	Carp::croak("FATAL ERROR: Unable to fetch files from $gencodePath (host $gencodeHost)")  unless(defined($localGTF));
		
	$ftpServer->disconnect()  if($ftpServer->can('disconnect'));
	$ftpServer->quit()  if($ftpServer->can('quit'));
	$ftpServer = undef;
	
	# This hash holds the seq_region_id -> name correspondence
	my %regionId = ();
	
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
	
	print "Parsing ",$localSeqRegion,"\n";
	if(open(my $SEQREG,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$localSeqRegion)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> \%regionId,
			TabParser::TAG_CALLBACK => \&parseSeqRegions,
			# The columns are 
			#	seq_region_id
			#	name
			TabParser::TAG_FETCH_COLS => [1,2],
		);
		TabParser::parseTab($SEQREG,%config);
		
		close($SEQREG);
	} else {
		print STDERR "ERROR: Unable to parse Sequence Regions file ".$localSeqRegion."\n";
	}
	
	print "Parsing ",$localGenes,"\n";
	if(open(my $ENSG,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$localGenes)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,\%ENShash,\%ENSintHash,'gene'],
			TabParser::TAG_CALLBACK => \&parseENS,
			# The columns are 
			#	seq_region_id
			#	chromosome_start
			#	chromosome_end
			#	display_xref_id
			#	ensembl gene id
			#	ensembl gene id version
			#	description
			TabParser::TAG_FETCH_COLS => [3,4,5,7,14,15,10],
#			TabParser::TAG_NEG_FILTER => [[1 => 'LRG_gene']],
		);
		TabParser::parseTab($ENSG,%config);
		
		close($ENSG);
	} else {
		print STDERR "ERROR: Unable to parse EnsEMBL Genes file ".$localGenes."\n";
	}
	
	print "Parsing ",$localTranscripts,"\n";
	if(open(my $ENST,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$localTranscripts)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> [\%regionId,\%ENShash,\%ENSintHash,'transcript'],
			TabParser::TAG_CALLBACK => \&parseENS,
			# The columns are 
			#	seq_region_id
			#	chromosome_start
			#	chromosome_end
			#	display_xref_id
			#	ensembl transcript id
			#	ensembl transcript id version
			#	description
			TabParser::TAG_FETCH_COLS => [3,4,5,7,13,14,10],
#			TabParser::TAG_NEG_FILTER => [[8 => 'LRG_gene']],
		);
		TabParser::parseTab($ENST,%config);
		
		close($ENST);
	} else {
		print STDERR "ERROR: Unable to open EnsEMBL Transcripts file ".$localTranscripts."\n";
	}

	print "Parsing ",$localXref,"\n";
	if(open(my $XREF,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$localXref)) {
		my %config = (
			TabParser::TAG_CONTEXT	=> \%ENSintHash,
			TabParser::TAG_CALLBACK => \&parseXREF,
			# The columns are 
			#	display_xref_id
			#	dbprimary_acc
			#	display_label
			#	description
			TabParser::TAG_FETCH_COLS => [0,2,3,5],
		);
		TabParser::parseTab($XREF,%config);
		
		close($XREF);
	} else {
		print STDERR "ERROR: Unable to open EnsEMBL XREF file ".$localXref."\n";
	}

	
	my @GencodeObjects = ();
	
	print "Parsing ",$localGTF,"\n";
	if(open(my $GTF,'-|',BP::Loader::CorrelatableConcept::GUNZIP,'-c',$localGTF)) {
		my %config = (
			TabParser::TAG_COMMENT	=>	'#',
			TabParser::TAG_CONTEXT	=> [\%ENShash,\@GencodeObjects],
			TabParser::TAG_CALLBACK => \&parseGTF,
		);
		TabParser::parseTab($GTF,%config);
		
		close($GTF);
	} else {
		print STDERR "ERROR: Unable to open EnsEMBL XREF file ".$localGTF."\n";
	}
	
	# Let's parse the model
	my $modelFile = $ini->val($BP::Loader::Mapper::SECTION,'model');
	# Setting up the right path on relative cases
	$modelFile = File::Spec->catfile(File::Basename::dirname($iniFile),$modelFile)  unless(File::Spec->file_name_is_absolute($modelFile));

	print "Parsing model $modelFile...\n";
	my $model = undef;
	eval {
		$model = BP::Model->new($modelFile);
	};
	
	if($@) {
		Carp::croak('ERROR: Model parsing and validation failed. Reason: '.$@);
	}
	print "\tDONE!\n";
	
	my %storageModels = ();
	
	# Setting up the loader storage model(s)
	Carp::croak('ERROR: undefined destination storage model')  unless($ini->exists($BP::Loader::Mapper::SECTION,'loaders'));
	my $loadModelNames = $ini->val($BP::Loader::Mapper::SECTION,'loaders');
	
	my @loadModels = ();
	foreach my $loadModelName (split(/,/,$loadModelNames)) {
		unless(exists($storageModels{$loadModelName})) {
			$storageModels{$loadModelName} = BP::Loader::Mapper->newInstance($loadModelName,$model,$ini);
			push(@loadModels,$loadModelName);
		}
	}
	# Now, do we need to push the metadata there?
	if(!$ini->exists($BP::Loader::Mapper::SECTION,'metadata-loaders') || $ini->val($BP::Loader::Mapper::SECTION,'metadata-loaders') eq 'true') {
		foreach my $mapper (@storageModels{@loadModels}) {
			$mapper->storeNativeModel();
		}
	}
	
	my $concept = $model->getConceptDomain('external')->conceptHash->{'gencode'};
	my $corrConcept = BP::Loader::CorrelatableConcept->new($concept);
	
	# Storage of the data
	foreach my $mapper (@storageModels{@loadModels}) {
		$mapper->setDestination($corrConcept);
		
		my $destination = $mapper->getInternalDestination();
		
		my $entorp = $mapper->validateAndEnactEntry(values(%ENShash),@GencodeObjects);
		use Data::Dumper;
		
		my $bulkData = $mapper->_bulkPrepare($entorp);
		$mapper->_bulkInsert($destination,$bulkData);
		
		$bulkData = undef;
		$entorp = undef;
		$destination = undef;
		$mapper->freeDestination();
	}
} else {
	print STDERR "Usage: $0 iniFile cachingDir\n"
}
