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

use Log::Log4perl;

my $LOG;
sub BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::WARN, layout => '%-5p - %m%n' } );

	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
}

use constant {
	GENCODE_FTP_BASE_TAG	=> 'gencode-ftp-base-uri',
	GENCODE_GTF_FILE_TAG	=> 'gencode-gtf',
};

use constant HAVANA_NS	=> 'HAVANA';
my %FeaturesNSMap = (
	'gene_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'gene_name'	=> 'HGNC',
	'havana_gene'	=> HAVANA_NS,
	'transcript_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'transcript_name'	=> BP::DCCLoader::Parsers::EnsemblGTParser::DESCRIPTION_NS,
	'havana_transcript'	=> HAVANA_NS,
	'exon_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'protein_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'ccdsid'	=>	'CCDS',
);

my %MainFeatures = (
	'gene' => undef,
	'transcript' => undef,
	'exon' => undef
);

my %FeaturesKeys = (
	'gene' => 'gene_id',
	'transcript' => 'transcript_id',
	'Selenocysteine' => 'transcript_id',
	'UTR' => 'transcript_id',
	'exon' => 'exon_id',
	'CDS' => 'exon_id',
	'start_codon' => 'exon_id',
	'stop_codon' => 'exon_id',
);

my %FeaturesClusterKeys = (
	'gene' => ['gene_id'],
	'transcript' => ['gene_id','transcript_id'],
	'Selenocysteine' => ['gene_id','transcript_id'],
	'UTR' => ['gene_id','transcript_id'],
	'exon' => ['gene_id','transcript_id','exon_id'],
	'CDS' => ['gene_id','transcript_id','exon_id'],
	'start_codon' => ['gene_id','transcript_id','exon_id'],
	'stop_codon' => ['gene_id','transcript_id','exon_id'],
);

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
		$frame, # frame
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
	if(exists($FeaturesKeys{$feature})) {
		$ensIdKey = $FeaturesKeys{$feature};
	} else {
		$ensIdKey = 'transcript_id';
	}
	
	my $fullEnsemblId = $attributes{$ensIdKey};
	my @ensemblIds = ($fullEnsemblId);
	my $ensemblId = undef;
	if(exists($MainFeatures{$feature})) {
		$ensemblId = substr($fullEnsemblId,0,index($fullEnsemblId,'.'));
		push(@ensemblIds,$ensemblId);
	} else {
		$fullEnsemblId .= '.'.$feature;
		$ensemblId = $fullEnsemblId;
	}
	if(exists($p_ENShash->{$ensemblId})) {
		$p_regionData = $p_ENShash->{$ensemblId};
		$local = undef;
	} else {
		my $term = $chroCV->getTerm($chro);
		if($term) {
			#print "DEBUG: $chro $feature\n"  unless(exists($FeaturesClusterKeys{$feature}));
			my @feature_cluster_ids = exists($FeaturesClusterKeys{$feature}) ? @attributes{@{$FeaturesClusterKeys{$feature}}} : ( $ensIdKey );
			$p_regionData = {
				#$fullEnsemblId,
				'feature_cluster_id'	=> \@feature_cluster_ids,
				'feature_ns'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
				'feature_id'	=> $fullEnsemblId,
				'feature'	=> $feature,
				'coordinates'	=> [
					{
						'feature_ns'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
						'feature_id'	=> $fullEnsemblId,
						'chromosome'	=> $term->key(),
						'chromosome_start'	=> ($chromosome_start+0),
						'chromosome_end'	=> ($chromosome_end+0),
						'chromsome_strand'	=> (($chromosome_strand eq '-') ? -1 : 1),
					}
				],
				'symbol'	=> [{
					'ns'	=>	BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
					'name'	=>	\@ensemblIds,
				}],
			};
			
			$p_regionData->{'coordinates'}[0]{'frame'} = $frame  unless($frame eq '.');
		#} else {
		#	print STDERR "DEBUG CHRO: $chro\n";
		}
	}
	
	if($p_regionData) {
		FEATURES:
		foreach my $feature (keys(%attributes)) {
			my $ns = exists($FeaturesNSMap{$feature}) ? $FeaturesNSMap{$feature} : $feature;
			my $val = $attributes{$feature};
			# Inefficient, but effective
			my $p_name = undef;
			foreach my $p_sym (@{$p_regionData->{'symbol'}}) {
				if($p_sym->{'ns'} eq $ns) {
					$p_name = $p_sym->{'name'};
					
					last;
				}
			}
			
			unless(defined($p_name)) {
				$p_name = [];
				push(@{$p_regionData->{'symbol'}},{'ns' => $ns, 'name' => $p_name});
			}
			
			# Inefficient, but effective
			foreach my $sym (@{$p_name}) {
				if($sym eq $val) {
					#$val = undef;
					
					next FEATURES;
				}
			}
			push(@{$p_name},$val);
		}
		
		# Last, save it!!!
		if($local) {
			push(@{$p_Gencode},$p_regionData);
			$p_ENShash->{$ensemblId} = $p_regionData;
			$p_PAR->{$ensemblId} = $p_regionData  if(exists($MainFeatures{$feature}));
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
	
	$LOG->info("Connecting to $gencode_ftp_base...");
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
	
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	my $p_ENShash = BP::DCCLoader::Parsers::EnsemblGTParser::getEnsemblCoordinates($model,$workingDir,$ini,$testmode);
	
	$LOG->info("Parsing ".$localGTF);
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
