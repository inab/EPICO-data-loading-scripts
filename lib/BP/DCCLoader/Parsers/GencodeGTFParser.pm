#!/usr/bin/perl

use warnings;
use strict;
use experimental 'smartmatch';

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

use constant {
	GENE_ID => 'gene_id',
	TRANSCRIPT => 'transcript',
	TRANSCRIPT_ID => 'transcript_id'
};

use constant {
	SEEN_ID_KEY => 'se',
	GENE_CACHE_KEY => 'gc',
	FEATURES_LIST_KEY => 'feat'
};

my %AttributesNSMap = (
	GENE_ID()	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'gene_name'	=> 'HGNC',
	'havana_gene'	=> HAVANA_NS,
	TRANSCRIPT_ID()	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'transcript_name'	=> BP::DCCLoader::Parsers::EnsemblGTParser::DESCRIPTION_NS,
	'havana_transcript'	=> HAVANA_NS,
	'exon_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'protein_id'	=> BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
	'ccdsid'	=>	'CCDS',
);

my %MainFeatures = (
	'gene' =>  {
		'gene_name'	=> undef,
		'havana_gene'	=> undef,
	},
	TRANSCRIPT() => {
		'transcript_name'	=> undef,
		'havana_transcript'	=> undef,
	},
	'exon' => undef
);

my %FeaturesKeys = (
	'gene' => GENE_ID,
	TRANSCRIPT() => TRANSCRIPT_ID,
	'Selenocysteine' => TRANSCRIPT_ID,
	'UTR' => TRANSCRIPT_ID,
	'exon' => 'exon_id',
	'CDS' => 'exon_id',
	'start_codon' => 'exon_id',
	'stop_codon' => 'exon_id',
);

my %FeaturesClusterKeys = (
	'gene' => [GENE_ID],
	TRANSCRIPT() => [GENE_ID,TRANSCRIPT_ID],
	'Selenocysteine' => [GENE_ID,TRANSCRIPT_ID],
	'UTR' => [GENE_ID,TRANSCRIPT_ID],
	'exon' => [GENE_ID,TRANSCRIPT_ID,'exon_id'],
	'CDS' => [GENE_ID,TRANSCRIPT_ID,'exon_id'],
	'start_codon' => [GENE_ID,TRANSCRIPT_ID,'exon_id'],
	'stop_codon' => [GENE_ID,TRANSCRIPT_ID,'exon_id'],
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
	
	my($p_ENShash,$p_APPRISannot,$chroCV,$batchProcessor,$p_batchCache) = @{$payload};
	
	$p_batchCache->{GENE_ID()} = ''  unless(exists($p_batchCache->{GENE_ID()}));
	
	my %attributes = ();
	
	#print STDERR "DEBUG: ",join(" ",@_),"\n"  unless(defined($attributes_str));
	my @tokens = split(/\s*;\s*/,$attributes_str);
	foreach my $token (@tokens) {
		my($key,$value) = split(/\s+/,$token,2);
		
		# Removing double quotes
		$value =~ tr/"//d;
		
		$attributes{$key} = $value;
	}
	
	# Micro-cache management
	my $currentGeneId = $attributes{GENE_ID()};
	my $p_currentCache;
	my $p_currentBatchData;
	if($currentGeneId eq $p_batchCache->{GENE_ID()}) {
		$p_currentCache = $p_batchCache->{GENE_CACHE_KEY()};
		$p_currentBatchData = $p_batchCache->{FEATURES_LIST_KEY()};
	} else {
		# Run the batchProcessor just here
		$batchProcessor->($p_batchCache->{FEATURES_LIST_KEY()})  if(defined($batchProcessor) && exists($p_batchCache->{FEATURES_LIST_KEY()}));
		
		$p_batchCache->{GENE_ID()} = $currentGeneId;
		$p_currentCache = {};
		$p_currentBatchData = [];
		$p_batchCache->{GENE_CACHE_KEY()} = $p_currentCache;
		$p_batchCache->{FEATURES_LIST_KEY()} = $p_currentBatchData;
	}
	
	my $p_regionData = undef;
	
	my $ensIdKey = undef;
	my $p_ensFeatures = undef;
	if(exists($FeaturesKeys{$feature})) {
		$ensIdKey = $FeaturesKeys{$feature};
	} else {
		$ensIdKey = TRANSCRIPT_ID;
	}
	
	my $fullEnsemblId = $attributes{$ensIdKey};
	my @ensemblIds;
	my $ensemblId = undef;
	if(exists($MainFeatures{$feature})) {
		$ensemblId = substr($fullEnsemblId,0,index($fullEnsemblId,'.'));
		@ensemblIds = ($fullEnsemblId,$ensemblId);
	} else {
		$fullEnsemblId .= '.'.$feature;
		$ensemblId = $fullEnsemblId;
		@ensemblIds = ($fullEnsemblId);
	}
	
	my $local;
	if(exists($p_currentCache->{$ensemblId})) {
		$p_regionData = $p_currentCache->{$ensemblId};
	} elsif(exists($MainFeatures{$feature}) && exists($p_ENShash->{$ensemblId})) {
		$p_regionData = $p_ENShash->{$ensemblId};
	} else {
		my $term = $chroCV->getTerm($chro);
		if($term) {
			$local = 1;
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
					'domain'	=>	BP::DCCLoader::Parsers::EnsemblGTParser::ENSEMBL_NS,
					'value'	=>	\@ensemblIds,
				}],
			};
			
			$p_regionData->{'coordinates'}[0]{'frame'} = $frame  unless($frame eq '.');
		#} else {
		#	print STDERR "DEBUG CHRO: $chro\n";
		}
	}
	
	if($p_regionData) {
		my $p_fea = exists($MainFeatures{$feature}) ? $MainFeatures{$feature} : {};
		foreach my $attribute (keys(%attributes)) {
			# Pairs of attribute type and attribute namespace
			my @attributePairs = (['attribute',$attribute]);
			if(exists($AttributesNSMap{$attribute})) {
				push(@attributePairs, [ exists($p_fea->{$attribute}) ? 'symbol' : 'attribute' , $AttributesNSMap{$attribute}]);
			}
			
			my $val = $attributes{$attribute};
			
			ATTRIBUTE_PAIR:
			foreach my $p_attributePair (@attributePairs) {
				my($attributeType,$ns) = @{$p_attributePair};
				# Inefficient, but effective
				foreach my $p_sym (@{$p_regionData->{$attributeType}}) {
					if($p_sym->{'domain'} eq $ns) {
						my $p_name = $p_sym->{'value'};
						
						# Inefficient, but effective
						push(@{$p_name},$val)  unless($val ~~ @{$p_name});
						
						next ATTRIBUTE_PAIR;
					}
				}
				
				push(@{$p_regionData->{$attributeType}},{'domain' => $ns, 'value' => [$val]});
			}
		}
		
		# APPRIS enrichment (only once!)
		if($feature eq TRANSCRIPT && exists($p_APPRISannot->{$ensemblId})) {
			push(@{$p_regionData->{'attribute'}},{'domain' => 'APPRIS_PRINCIPAL', 'value' => [$p_APPRISannot->{$ensemblId}]});
			#delete($p_APPRISannot->{$ensemblId});
		}
		
		# Last, save it!!!
		unless(exists($p_currentCache->{$ensemblId})) {
			# Saving it for the micro-cache
			$p_batchCache->{SEEN_ID_KEY()} = {}  unless(exists($p_batchCache->{SEEN_ID_KEY()}));
			$p_batchCache->{SEEN_ID_KEY()}{$ensemblId} = undef;
			$p_currentCache->{$ensemblId} = $p_regionData;
			push(@{$p_currentBatchData},$p_regionData);
		}
		
		# This is needed in order to avoid missing genes, transcripts or exons
		if($local && exists($MainFeatures{$feature})) {
			$p_ENShash->{$ensemblId} = $p_regionData;
		}
	}
	
	1;
}

sub getGencodeCoordinates($$$;$$$) {
	my($model,$workingDir,$ini,$testmode,$p_APPRISannot,$batchProcessor) = @_;
	
	my $gencode_ftp_base = undef;
	my $gencode_gtf_file = undef;
	
	# Check the needed parameters for the construction
	if($ini->exists(BP::DCCLoader::Parsers::EPICO_LOADER_SECTION,GENCODE_FTP_BASE_TAG)) {
		$gencode_ftp_base = $ini->val(BP::DCCLoader::Parsers::EPICO_LOADER_SECTION,GENCODE_FTP_BASE_TAG);
	} else {
		$LOG->logdie("Configuration file must have '".GENCODE_FTP_BASE_TAG."' in '".BP::DCCLoader::Parsers::EPICO_LOADER_SECTION."' section");
	}
	
	if($ini->exists(BP::DCCLoader::Parsers::EPICO_LOADER_SECTION,GENCODE_GTF_FILE_TAG)) {
		$gencode_gtf_file = $ini->val(BP::DCCLoader::Parsers::EPICO_LOADER_SECTION,GENCODE_GTF_FILE_TAG);
	} else {
		$LOG->logdie("Configuration file must have '".GENCODE_GTF_FILE_TAG."' in '".BP::DCCLoader::Parsers::EPICO_LOADER_SECTION."' section");
	}
	
	# Now, let's patch the properies of the different remote resources, using the properties inside the model
	eval {
		$model->annotations->applyAnnotations(\($gencode_ftp_base,$gencode_gtf_file));
	};
	
	if($@) {
		$LOG->logdie("$@ (does not exist in model)");
	}
	
	# And translate these to URI objects
	$gencode_ftp_base = URI->new($gencode_ftp_base);
	
	# Defined outside
	my $ftpServer = undef;
	
	$LOG->info("Connecting to $gencode_ftp_base...");
	my $gencodeHost = $gencode_ftp_base->host();
	$ftpServer = Net::FTP::AutoReconnect->new($gencodeHost,Debug=>0) || $LOG->logdie("FTP connection to server ".$gencodeHost." failed: ".$@);
	$ftpServer->login(BP::DCCLoader::WorkingDir::ANONYMOUS_USER,BP::DCCLoader::WorkingDir::ANONYMOUS_PASS) || $LOG->logdie("FTP login to server $gencodeHost failed: ".$ftpServer->message());
	$ftpServer->binary();
	
	my $gencodePath = $gencode_ftp_base->path();
	
	my($localGTF, $gencodePathReason) = $workingDir->cachedGet($ftpServer,$gencodePath.'/'.$gencode_gtf_file);
	
	$LOG->logdie("FATAL ERROR: Unable to fetch files from $gencodePath (host $gencodeHost). Reason: ".$gencodePathReason)  unless(defined($localGTF));
		
	$ftpServer->disconnect()  if($ftpServer->can('disconnect'));
	$ftpServer->quit()  if($ftpServer->can('quit'));
	$ftpServer = undef;
	
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	my $p_ENShash = BP::DCCLoader::Parsers::EnsemblGTParser::getEnsemblCoordinates($model,$workingDir,$ini,$testmode);
	#my $p_ENShash;
	#use Test::LeakTrace;
	#my $leakedRefs = leaked_info {
	#	$p_ENShash = BP::DCCLoader::Parsers::EnsemblGTParser::getEnsemblCoordinates($model,$workingDir,$ini,$testmode);
	#};
	#foreach my $leak (@{$leakedRefs}) {
	#	print STDERR "DEBUG: Leaked: $leak->[1] : $leak->[2]\n";
	#}
	#exit(1);
	
	$LOG->info("Parsing ".$localGTF);
	if(open(my $GTF,'-|',BP::Loader::Tools::GUNZIP,'-c',$localGTF)) {
		my %batchCache;
		my %config = (
			TabParser::TAG_COMMENT	=>	'#',
			# Fourth parameter is a temporal storage used along the full parsing process
			TabParser::TAG_CONTEXT	=> [$p_ENShash,$p_APPRISannot, $chroCV, $batchProcessor, \%batchCache],
			TabParser::TAG_CALLBACK => \&parseGTF,
		);
		$config{TabParser::TAG_VERBOSE} = 1  if($testmode);
		TabParser::parseTab($GTF,%config);
		
		close($GTF);
		
		# And last steps!
		if(defined($batchProcessor) && exists($batchCache{FEATURES_LIST_KEY()})) {
			# First, the last registered ones
			$batchProcessor->($batchCache{FEATURES_LIST_KEY()});
			# And then, the ones which are not in Gencode
			my @onlyInEnsembl = ();
			my $p_seenEnsIds = $batchCache{SEEN_ID_KEY()};
			foreach my $ensemblId (keys(%{$p_ENShash})) {
				push(@onlyInEnsembl,$p_ENShash->{$ensemblId})  unless(exists($p_seenEnsIds->{$ensemblId}));
			}
			$batchProcessor->(\@onlyInEnsembl);
		}
	} else {
		$LOG->logdie("ERROR: Unable to open EnsEMBL XREF file ".$localGTF);
	}
	
	return $p_ENShash;
}

1;
