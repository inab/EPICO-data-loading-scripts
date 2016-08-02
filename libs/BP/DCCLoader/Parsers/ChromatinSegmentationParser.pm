#!/usr/bin/perl

use warnings;
use strict;

use Carp;

use TabParser;

package BP::DCCLoader::Parsers::ChromatinSegmentationParser;

use base qw(BP::DCCLoader::Parsers::AbstractInsertionParser);

use File::Spec;
use File::Temp qw();

use constant CNIO_CHROMATIN_SEGMENTATION_METADATA => {
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
			'name'	=>	'ChromHMM',
			'url'	=>	'http://compbio.mit.edu/ChromHMM/',
		},
	],
};

use constant EXPERIMENT_CONCEPT_NAME	=>	'cs';
use constant ANALYSIS_CONCEPT_DOMAIN_NAME	=>	'chromseg';
use constant ANALYSIS_CONCEPT_NAME	=>	'c';
use constant CHROMATIN_GROUP_ID	=>	'6';

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

sub _insertInternal($$) {
	my($self)=shift;
	
	my($analysis_id,$p_insertMethod) = @_;
	
	my $model = $self->{BP::DCCLoader::Parsers::AbstractInsertionParser::K_MODEL};
	my $chroCV = $model->getNamedCV('ChromosomesAndScaffolds');
	my $termChro = '';
	my $term = undef;
	my $chromosome = undef;
	
	my %chromHMMBedParserConfig = (
		TabParser::TAG_CALLBACK => sub {
			my(
				$chro,
				$chromosome_start,
				$chromosome_end,
				$state_str,
			) = @_;
			
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
				my $state;
				if($state_str =~ /^E([1-9][0-9]*)$/) {
					$state = $1 + 0;
				}
								
				my %entry = (
					'analysis_id'	=>	$analysis_id,
					'chromosome'	=>	$chromosome,
					'chromosome_start'	=>	$chromosome_start+1,	# Bed holds the data 0-based
					'chromosome_end'	=>	$chromosome_end+0,	# Bed holds the end coordinate as exclusive, so it does not change
					'state'	=>	$state,
				);
				
				$p_insertMethod->(\%entry);
			}
		},
	);
	return \%chromHMMBedParserConfig;
}

use constant {
	CHROMHMM_PATH	=>	'secondary_analysis/Segmentation_of_ChIP-Seq_data',
	CHROMHMM_MD5_PATTERN	=>	qr/^ChromHMM.*\.md5sum$/,
	CHROMHMM_README_PATTERN	=>	qr/^README/,
	CHROMHMM_SEGMENTATION_PATTERN	=>	qr/^SEGMENTATION_.*\.tar\.gz$/,
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
	my $chromhmm_path = join('/',$metadataParser->getIndexPath(),CHROMHMM_PATH());
	$LOG->info("Analyzing files from ".$chromhmm_path.'...');
	my $p_chromHMMFiles = $metadataParser->dataServerListing($chromhmm_path);
	
	if(defined($p_chromHMMFiles)) {
		my $publicIndexPayload = $metadataParser->getPublicIndexPayload();
		
		my $workingDir = $publicIndexPayload->{workingDir};
		my $bpDataServer = $publicIndexPayload->{bpDataServer};
		
		# And process it (getting the MD5 file BTW)
		my $md5filePath = undef;
		my $readmePath = undef;
		my $segmentationPath = undef;
		foreach my $chromHMMFile (@{$p_chromHMMFiles}) {
			my $shortPath =  __getShortPath($chromHMMFile);
			if(!defined($md5filePath) && $shortPath =~ CHROMHMM_MD5_PATTERN()) {
				$md5filePath = $chromHMMFile;
			} elsif(!defined($readmePath) && $shortPath =~ CHROMHMM_README_PATTERN()) {
				$readmePath = $chromHMMFile;
			} elsif(!defined($segmentationPath) && $shortPath =~ CHROMHMM_SEGMENTATION_PATTERN()) {
				$segmentationPath = $chromHMMFile;
			}
		}
		
		# Now, let's fetch the MD5 file
		my %MD5SUM = ();
		if(defined($md5filePath)) {
			$LOG->info("\t* Fetching $md5filePath");
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
		
		# The README, which contains the states' descriptions
		my @states = ();
		if(defined($readmePath)) {
			$LOG->info("\t* Fetching $readmePath");
			my $shortPath = __getShortPath($readmePath);
			my($localReadmeFile, $reason) = $workingDir->cachedGet($bpDataServer,$readmePath,undef,exists($MD5SUM{$shortPath}) ? $MD5SUM{$shortPath} : undef);
			
			if(defined($localReadmeFile)) {
				if(open(my $README,'<',$localReadmeFile)) {
					while(my $line=<$README>) {
						chomp($line);
						if($line =~ /State([1-9][0-9]*)[^ ]+ (.+)$/) {
							push(@states,{'state' => $1 + 0, 'annotation' => $2});
						}
					}
					
					close($README);
				}
			} else {
				$LOG->logdie("Error while fetching README file $readmePath . Reason: ".$reason);
			}
		} else {
			$LOG->logdie("No README file found");
		}
		
		# Preparation of common variables
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
		
		# Computing correspondence between sample id and experiment id for ChIP-Seq experiments
		my %sId2ExpIds = ();
		foreach my $p_exp (@{$publicIndexPayload->{'lab'}{EXPERIMENT_CONCEPT_NAME()}}) {
			#my $p_exp = $publicIndexPayload->{experiments}{$experiment_id};
			#if(!defined($p_exp)) {
			#	$LOG->info("JODER $experiment_id");
			#} elsif(!exists($p_exp->{'library_strategy'})) {
			#	use Data::Dumper;
			#	$LOG->info("JARL",Dumper($p_exp));
			#}
			
			#if($p_exp->{'library_strategy'} eq 'ChIP-Seq') {
				push(@{$sId2ExpIds{$p_exp->{'analyzed_sample_id'}}},$p_exp->{'experiment_id'});
			#}
		}
		
		# The segmentation archive
		if(defined($segmentationPath)) {
			$LOG->info("\t* Fetching $segmentationPath");
			my $shortPath = __getShortPath($segmentationPath);
			my($localSegmentationArchive, $reason) = $workingDir->cachedGet($bpDataServer,$segmentationPath,undef,exists($MD5SUM{$shortPath}) ? $MD5SUM{$shortPath} : undef);
			
			if(defined($localSegmentationArchive)) {
				my $tempDir = $workingDir->createTempDir();
				# Avoiding to lose the temporary directory
				$instance->{__tempChromDir__} = $tempDir;
				my $tempDirName = $tempDir->dirname();
				
				$LOG->info("\t* Temporary extracting $localSegmentationArchive");
				if(system('tar','-x','-z','-C',$tempDirName,'-f',$localSegmentationArchive)==0) {
					if(opendir(my $SEGDIR,$tempDirName)) {
						while(my $bedRelPath = readdir($SEGDIR)) {
							if($bedRelPath =~ /^(.+)_[1-9][0-9]*_Blueprint.*\.bed$/) {
								my $sample_id = $1;
								my $p_experiment_ids = [];
								if(exists($sName2Id{$sample_id})) {
									# First, the proper sample_id
									$sample_id = $sName2Id{$sample_id};
									
									# Then, filter by ChIP-Seq experiments
									$p_experiment_ids = $sId2ExpIds{$sample_id}  if(exists($sId2ExpIds{$sample_id}));
								}
								
								my $analysis_id = $sample_id . '_chromseg';
								
								# The metadata of the analysis
								my %a_metadata = (
									'analysis_id'	=>	$analysis_id,
									'experiment_id'	=>	$p_experiment_ids,
									'analyzed_sample_id'	=>	$sample_id,
									'analysis_group_id'	=>	[ CHROMATIN_GROUP_ID() ],
									'data_status'	=>	2,
									'assembly_version'	=>	$grch_version,
									'ensembl_version'	=>	$ensembl_version,
									'gencode_version'	=>	$gencode_version,
								);
								
								@a_metadata{keys(%{CNIO_CHROMATIN_SEGMENTATION_METADATA()})} = values(%{CNIO_CHROMATIN_SEGMENTATION_METADATA()});
								
								push(@{$p_anal->{ANALYSIS_CONCEPT_DOMAIN_NAME()}},\%a_metadata);
								
								# And now, pushing the fileset
								
								push(@{$p_primary_anal->{ANALYSIS_CONCEPT_DOMAIN_NAME()}},[
									$analysis_id,
									ANALYSIS_CONCEPT_NAME(),
									$instance,
									{
										'l_file'	=>	File::Spec->catfile($tempDirName,$bedRelPath),
									}
								]);
							}
						}
						
						closedir($SEGDIR);
					} else {
						$LOG->logdie("Unable to open directory $tempDirName. Reason: ".$!);
					}
				} elsif ($? == -1) {
					$LOG->logdie("failed to execute: $!");
				} elsif ($? & 127) {
					$LOG->error(sprintf("child died with signal %d, %s coredump\n", ($? & 127),  ($? & 128) ? 'with' : 'without'));
					exit 1;
				} else {
					$LOG->error(sprintf("child exited with value %d\n", $? >> 8));
					exit 1;
				}
			} else {
				$LOG->logdie("Error while fetching SEGMENTATION archive $segmentationPath . Reason: ".$reason);
			}
		} else {
			$LOG->logdie("No SEGMENTATION archive found");
		}
	} else {
		$LOG->warn("Path not found $chromhmm_path");
	}
}

1;
