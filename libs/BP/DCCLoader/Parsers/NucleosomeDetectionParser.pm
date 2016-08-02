#!/usr/bin/perl

use warnings;
use strict;

use Carp;

package BP::DCCLoader::Parsers::NucleosomeDetectionParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use constant MPI_NUCLEODETECT_METADATA => {
	'program_versions'	=>	[
		{
			'program'	=>	'BWA',
			'version'	=>	'0.5.9',
		},
	],
	'alignment_algorithm'	=>	{
		'name'	=>	'BWA',
		'url'	=>	'http://bio-bwa.sourceforge.net/',
	},
	'other_analysis_algorithm'	=>	[
		{
			'name'	=>	'NucHunter',
			'url'	=>	'http://epigen.molgen.mpg.de/nuchunter/',
		},
	],
};

use constant TSV_EXPERIMENT_COLUMN_START	=>	8;
use constant ANALYSIS_CONCEPT_DOMAIN_NAME	=>	'nucleosome';
use constant ANALYSIS_CONCEPT_NAME	=>	'c';
use constant NUCLEOSOME_GROUP_ID	=>	'12c';

use Tie::IxHash;

use Log::Log4perl;

my $LOG;
sub BEGIN {
	Log::Log4perl->easy_init( { level => $Log::Log4perl::WARN, layout => '%-5p - %m%n' } );

	$LOG = Log::Log4perl->get_logger(__PACKAGE__);
}

# This is the empty constructor
sub new(;\%) {
	my($self)=shift;
	my($class)=ref($self) || $self;
	
	$self = $class->SUPER::new(@_)  unless(ref($self));
	
	return $self;
}

# _complexInsertInternal parameters:
#	analysis_id: The analysis_id for each entry
#	p_insertMethod: We feed this method with the prepared entries
#	p_files: The involved fileset
# It does the job
sub _complexInsertInternal($$\@) {
	my($self)=shift;
	
	my($analysis_id, $p_insertMethod, $p_files) = @_;
	
	my $tsvFile;
	my $bedFile;
	foreach my $localFile (@{$p_files}) {
		if($localFile =~ /\.tsv$/) {
			$tsvFile = $localFile;
		} else {
			$bedFile = $localFile;
		}
	}
	
	if(open(my $TSV,'<',$tsvFile)) {
		# Process header
		my $tLine = <$TSV>;
		chomp($tLine);
		
		my @headers = split(/\t/,$tLine);
		# Getting the different positions
		my $readsBaseIndex = TSV_EXPERIMENT_COLUMN_START;
		my $activeBaseIndexRel = (scalar(@headers) - $readsBaseIndex)  / 2;
		
		my $pos = 0;
		my @reads_cols = ();
		foreach my $header (@headers[$readsBaseIndex..($readsBaseIndex + $activeBaseIndexRel - 1)]) {
			my @toks = split(/\./,$header);
			
			# Storing the experiment_id, the column for the reads and the column for the active state
			push(@reads_cols,[$toks[1],$pos,$activeBaseIndexRel+$pos]);
			$pos ++;
		}
		
		if(open(my $BED,'<',$bedFile)) {
			my $model = $self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL};
			my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
			my $termChro = '';
			my $term = undef;
			my $chromosome = undef;
			
			# Discard header
			my $bLine = <$BED>;
			
			while($bLine = <$BED>) {
				chomp($bLine);
				my($chro,$chromosome_start,$chromosome_end,$nucleosome_id) = split(/\t/,$bLine);
				
				# Synchronous read
				$tLine = <$TSV>;
				chomp($tLine);
				
				if($termChro ne $chro) {
					$termChro = $chro;
					$term = $chroCV->getTerm($chro);
					if($term) {
						$chromosome = $term->key();
					} else {
						print STDERR "\tdiscarding entries from unknown chromosome $chro\n";
						$chromosome = undef;
					}
				}
				
				if(defined($chromosome)) {
					my($t_nucleosome_id,undef,$position,$z_score,$total_reads,$noise,$p_value,$fuzziness,@readsByExp) = split(/\t/,$tLine);
					
					#$LOG->error("JAJA ".$bLine);
					#$LOG->error("JOJO ".join(",",@readsByExp));
					
					my @reads_by_experiment = map { {'experiment_id' => $_->[0], 'num_reads' => int($readsByExp[$_->[1]]), 'is_active' => (($readsByExp[$_->[2]] == '1')? boolean::true : boolean::false) } } @reads_cols;
					
					my %entry = (
						'analysis_id'	=>	$analysis_id,
						'nucleosome_id'	=>	$nucleosome_id,
						'chromosome'	=>	$chromosome,
						'chromosome_start'	=>	$chromosome_start+1,	# Bed holds the data 0-based
						'chromosome_end'	=>	$chromosome_end+0,	# Bed holds the end coordinate as exclusive, so it does not change
						'position'	=>	$position+1,	# TSV holds the position 0-based
						'z_score'	=>	$z_score + 0e0,
						'total_reads'	=>	$total_reads + 0,
						'noise'	=>	$noise + 0e0,
						'p_value'	=>	$p_value + 0e0,
						'fuzziness'	=>	$fuzziness + 0e0,
						'reads_by_experiment'	=>	\@reads_by_experiment,
					);
					
					$p_insertMethod->(\%entry);
				}
			}
			
			close($BED);
		} else {
			$LOG->warn("Unable to open $bedFile for analysis $analysis_id");
		}
		
		close($TSV);
	} else {
		$LOG->warn("Unable to open $tsvFile for analysis $analysis_id");
	}
}

use constant {
	NUCDET_PATH	=>	'secondary_analysis/Nucleosome_detection_and_histone_annotation',
	NUCDET_MD5	=>	'md5s.txt',
};

sub __getShortPath($) {
	my $path = shift;
	
	my $rpos = rindex($path,'/');
	if($rpos!=-1) {
		$path = substr($path,$rpos+1);
	}
	
	return $path;
}

sub GetAnalysisMetadata($\%\%) {
	my($metadataParser,$p_anal,$p_primary_anal) = @_;
	
	# First, let's get the filelist
	my $nucdet_path = join('/',$metadataParser->getIndexPath(),NUCDET_PATH());
	$LOG->info("Analyzing files from ".$nucdet_path.'...');
	my $p_nuclDetectFiles = $metadataParser->dataServerListing($nucdet_path);
	
	if(defined($p_nuclDetectFiles)) {
		my %pairs = ();
		tie(%pairs,'Tie::IxHash');
		
		my $publicIndexPayload = $metadataParser->getPublicIndexPayload();
		
		my $workingDir = $publicIndexPayload->{workingDir};
		my $bpDataServer = $publicIndexPayload->{bpDataServer};
		
		# And process it (getting the MD5 file BTW)
		my $md5filePath = undef;
		foreach my $nuclDetectFile (@{$p_nuclDetectFiles}) {
			if(!defined($md5filePath) && index($nuclDetectFile,NUCDET_MD5()) != -1) {
				$md5filePath = $nuclDetectFile;
			} elsif($nuclDetectFile =~ /\/([^\/]+)\.(bed|tsv)$/) {
				my $baseAnal = $1;
				my $ext = $2;
				
				$pairs{$baseAnal}{$ext} = $nuclDetectFile;
			}
		}
		
		# Now, let's fetch the MD5 file
		my %MD5SUM = ();
		if(defined($md5filePath)) {
			my($localMD5file, $reason) = $workingDir->cachedGet($bpDataServer,$md5filePath);
			
			if(defined($localMD5file)) {
				if(open(my $MD5file,'<',$localMD5file)) {
					while(my $line=<$MD5file>) {
						chomp($line);
						my @toks = split(/ +/,$line,3);
						$MD5SUM{$toks[1]} = $toks[0];
					}
					
					close($MD5file);
				}
			} else {
				$LOG->logdie("Error while fetching MD5 file $md5filePath . Reason: ".$reason);
			}
		} else {
			$LOG->warn("No MD5 file found");
		}
		
		my $instance = __PACKAGE__->new({
			BP::DCCLoader::Parsers::AbstractInsertionParser::K_INI		=>	$metadataParser->getConfigIni(),
			BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL	=>	$metadataParser->getModel(),
			BP::DCCLoader::Parsers::AbstractInsertionParser::K_WORKINGDIR	=>	$publicIndexPayload->{workingDir},
			BP::DCCLoader::Parsers::AbstractInsertionParser::K_TESTMODE	=>	$publicIndexPayload->{testmode},
		});
		my $grch_version = $metadataParser->getGRChVersion();
		my $ensembl_version = $metadataParser->getEnsemblVersion();
		my $gencode_version = $metadataParser->getGencodeVersion();
		
		# Computing correspondence between sample name and sample id
		my %sName2Id = ();
		foreach my $p_sample (values(%{$publicIndexPayload->{samples}})) {
			$sName2Id{$p_sample->{sample_name}} = $p_sample->{sample_id};
		}
		
		# Once it is processed, let's consider only those which have both tsv and bed
		foreach my $baseAnal (keys(%pairs)) {
			my $p_baseAnal = $pairs{$baseAnal};
			
			if(scalar(keys(%{$p_baseAnal})) >= 2) {
				# Let's fetch them
				my $tsvFile;
				my $bedFile;
				foreach my $path (values(%{$p_baseAnal})) {
					$LOG->info("\t* Fetching $path");
					my $shortPath = __getShortPath($path);
					my($localFile,$reason) = $workingDir->cachedGet($bpDataServer,$path,undef,exists($MD5SUM{$shortPath}) ? $MD5SUM{$shortPath} : undef);
					
					if(defined($localFile)) {
						if($localFile =~ /\.tsv$/) {
							$tsvFile = $localFile;
						} else {
							$bedFile = $localFile;
						}
					} else {
						$LOG->warn("Unable to fetch $path . Reason: ".$reason);
						last;
					}
				}
				
				# And assemble them in the memory structures
				if(defined($tsvFile) && defined($bedFile)) {
					# We need the header, in order to learn about the related experiments
					if(open(my $TSV,'<',$tsvFile)) {
						my $tLine = <$TSV>;
						chomp($tLine);
						
						my @headers = split(/\t/,$tLine);
						# Getting the different positions
						my $readsBaseIndex = TSV_EXPERIMENT_COLUMN_START;
						my $activeBaseIndexRel = (scalar(@headers) - $readsBaseIndex)  / 2;
						
						my @experiment_ids = ();
						my $sample_id = undef;
						my $analysis_id = undef;
						foreach my $header (@headers[$readsBaseIndex..($readsBaseIndex + $activeBaseIndexRel - 1)]) {
							my @toks = split(/\./,$header);
							
							# Storing the experiment_id, the column for the reads and the column for the active state
							push(@experiment_ids,$toks[1]);
							unless(defined($sample_id)) {
								my $sample_name = $toks[0];
								# Correspondence between sample name and id
								$sample_id = $sName2Id{$sample_name};
								
								$analysis_id = $sample_id .'_nucleosome';
							}
						}
						
						# The metadata of the analysis
						my %a_metadata = (
							'analysis_id'	=>	$analysis_id,
							'experiment_id'	=>	\@experiment_ids,
							'analyzed_sample_id'	=>	$sample_id,
							'analysis_group_id'	=>	[ NUCLEOSOME_GROUP_ID() ],
							'data_status'	=>	2,
							'assembly_version'	=>	$grch_version,
							'ensembl_version'	=>	$ensembl_version,
							'gencode_version'	=>	$gencode_version,
						);
						
						@a_metadata{keys(%{MPI_NUCLEODETECT_METADATA()})} = values(%{MPI_NUCLEODETECT_METADATA()});
						
						push(@{$p_anal->{ANALYSIS_CONCEPT_DOMAIN_NAME()}},\%a_metadata);
						
						# And now, pushing the fileset
						
						push(@{$p_primary_anal->{ANALYSIS_CONCEPT_DOMAIN_NAME()}},[
							$analysis_id,
							ANALYSIS_CONCEPT_NAME(),
							$instance,
							[
								{
									'l_file'	=>	$tsvFile,
								},
								{
									'l_file'	=>	$bedFile,
								}
							]
						]);
					} else {
						$LOG->warn("Unable to open file $tsvFile. Reason: ".$!);
					}
				}
			}
		}
	} else {
		$LOG->warn("Path not found $nucdet_path");
	}
}

1;
