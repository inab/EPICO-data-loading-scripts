#!/usr/bin/perl

use warnings;
use strict;

use BP::DCCLoader::Parsers::EnsemblGTParser;
use BP::DCCLoader::WorkingDir;

use BP::Loader::Tools;

use Carp;

use Net::FTP::AutoReconnect;

use TabParser;

use URI;

package BP::DCCLoader::Parsers::GencodeGTFParser;

use constant {
	GENCODE_FTP_BASE_TAG	=> 'gencode-ftp-base-uri',
	GENCODE_GTF_FILE_TAG	=> 'gencode-gtf',
};

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
	
	my($p_ENShash,$p_Gencode,$p_PAR,$chroCV) = @{$payload};
	
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
	my $local = 1;
	
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
	
	if($feature eq 'gene' || $feature eq 'transcript') {
		if(exists($p_ENShash->{$ensemblId})) {
			$p_regionData = $p_ENShash->{$ensemblId};
			$local = undef;
		}
	}
	
	if($local) {
		my $term = $chroCV->getTerm($chro);
		if($term) {
			$p_regionData = {
				#$fullEnsemblId,
				'feature_cluster_id'	=> $attributes{'gene_id'},
				'chromosome'	=> $term->key(),
				'chromosome_start'	=> $chromosome_start,
				'chromosome_end'	=> $chromosome_end,
				'symbol'	=> [$fullEnsemblId,$ensemblId],
				'feature'	=> $feature
			};
		#} else {
		#	print STDERR "DEBUG CHRO: $chro\n";
		}
	}
	
	if($p_regionData) {
		foreach my $ensFeature (@{$p_ensFeatures}) {
			push(@{$p_regionData->{symbol}},$attributes{$ensFeature})  if(exists($attributes{$ensFeature}));
		}
		
		# Last, save it!!!
		if($local) {
			push(@{$p_Gencode},$p_regionData);
			$p_PAR->{$ensemblId} = $p_regionData  if($feature eq 'gene' || $feature eq 'transcript');
		}
	}
	
	1;
}

sub getGencodeCoordinates($$$;$) {
	my($model,$workingDir,$ini,$testmode) = @_;
	
	my $gencode_ftp_base = undef;
	my $gencode_gtf_file = undef;
	
	# Check the needed parameters for the construction
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,GENCODE_FTP_BASE_TAG)) {
		$gencode_ftp_base = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,GENCODE_FTP_BASE_TAG);
	} else {
		Carp::croak("Configuration file must have '".GENCODE_FTP_BASE_TAG."' in '".BP::DCCLoader::Parsers::DCC_LOADER_SECTION."' section");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,GENCODE_GTF_FILE_TAG)) {
		$gencode_gtf_file = $ini->val(BP::DCCLoader::Parsers::DCC_LOADER_SECTION,GENCODE_GTF_FILE_TAG);
	} else {
		Carp::croak("Configuration file must have '".GENCODE_GTF_FILE_TAG."' in '".BP::DCCLoader::Parsers::DCC_LOADER_SECTION."' section");
	}
	
	# Now, let's patch the properies of the different remote resources, using the properties inside the model
	eval {
		$model->annotations->applyAnnotations(\($gencode_ftp_base,$gencode_gtf_file));
	};
	
	if($@) {
		Carp::croak("$@ (does not exist in model)");
	}
	
	# And translate these to URI objects
	$gencode_ftp_base = URI->new($gencode_ftp_base);
	
	# Defined outside
	my $ftpServer = undef;
	
	print "Connecting to $gencode_ftp_base...\n";
	my $gencodeHost = $gencode_ftp_base->host();
	$ftpServer = Net::FTP::AutoReconnect->new($gencodeHost,Debug=>0) || Carp::croak("FTP connection to server ".$gencodeHost." failed: ".$@);
	$ftpServer->login(BP::DCCLoader::WorkingDir::ANONYMOUS_USER,BP::DCCLoader::WorkingDir::ANONYMOUS_PASS) || Carp::croak("FTP login to server $gencodeHost failed: ".$ftpServer->message());
	$ftpServer->binary();
	
	my $gencodePath = $gencode_ftp_base->path();
	
	my $localGTF = $workingDir->cachedGet($ftpServer,$gencodePath.'/'.$gencode_gtf_file);
	
	Carp::croak("FATAL ERROR: Unable to fetch files from $gencodePath (host $gencodeHost)")  unless(defined($localGTF));
		
	$ftpServer->disconnect()  if($ftpServer->can('disconnect'));
	$ftpServer->quit()  if($ftpServer->can('quit'));
	$ftpServer = undef;
	
	my $chroCV = $model->getNamedCV('EnsemblChromosomes');
	my $p_ENShash = BP::DCCLoader::Parsers::EnsemblGTParser::getEnsemblCoordinates($model,$workingDir,$ini,$testmode);
	
	print "Parsing ",$localGTF,"\n";
	my @Gencode = ();
	my %PAR = ();
	if(open(my $GTF,'-|',BP::Loader::Tools::GUNZIP,'-c',$localGTF)) {
		my %config = (
			TabParser::TAG_COMMENT	=>	'#',
			TabParser::TAG_CONTEXT	=> [$p_ENShash,\@Gencode, \%PAR, $chroCV],
			TabParser::TAG_CALLBACK => \&parseGTF,
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($GTF,%config);
		
		close($GTF);
	} else {
		Carp::croak("ERROR: Unable to open EnsEMBL XREF file ".$localGTF);
	}
	
	return wantarray ? (\@Gencode,\%PAR,$p_ENShash) : \@Gencode;
}

1;
