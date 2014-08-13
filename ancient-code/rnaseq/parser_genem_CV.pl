#!/usr/bin/perl


use constants;
use XML::LibXML;
use utf8;
use strict;

my $file = shift;
my $release = shift;


# TODO : canviar path si és necessari
my $path = "/data2/INB/blueprint/ftp/"; # my $path = "/data2/INB/users/rroyo/roadmap/releases/";
my $path_IHEC = $path . $release . "/blueprint/metadata_dump/experiments";

binmode(STDOUT, ":utf8");

my $cmd;

open(EXP_M,">exp_m.txt");
open(EXP_M2,">exp_m.txt.all");
binmode(EXP_M2,":utf8");

open(LAB,">lab_mrna.txt");
open(LAB2,">lab_mrna.txt.all");
binmode(LAB2, ":utf8");

#TODO : modificar segons el model :: DB = exp_m
my $gene_m_str = "analysis_id\texperiment_id\tanalysis_group_id\tdata_status\tassembly_version\tensembl_version\tgencode_version\tprogram_versions\talignment_algorithm\tother_analysis_algorithm\tnotes\tanalysis_group_id_term\tdata_status_term\tassembly_version_term\n";

#TODO : modificar segons el model :: DB = lab_mrna
my $gene_lab_str = "experiment_id\tanalyzed_sample_id\texperiment_type\tlibrary_strategy\texperimental_group_id\tfeatures\traw_data_repository\traw_data_accession\tplatform\tseq_coverage\textraction_protocol\texperimental_group_id_term\traw_data_repository_term\tplatform_term\n";

print EXP_M $gene_m_str;

print LAB $gene_lab_str;

#lam_mrna mRNA-seq


open(F,$file);

my $str2;

while(<F>){
	my @arr = split('\t',$_);

        my $experiment_id = $arr[AN_ID];
	
	my $sfile = $path_IHEC . "/" . substr($experiment_id,0,6) . "/" . $experiment_id . ".xml";
#print $sfile;

        my $str_xml = `cat $sfile`;
        my $ihec;

        my $parser = XML::LibXML->new();
        my $doc;

        #TODO: afegir el nou atribut library_strategy
        my $library_strategy;


        if (-e $sfile ){
                eval { $doc = $parser->parse_string( $str_xml ) };
                die("couldn't parse XML '$str_xml' because\n\t$@") if $@;

                my $xml = $doc->getDocumentElement();

				#TODO: afegir el nou atribut library_strategy
                $library_strategy = $xml->getElementsByLocalName("LIBRARY_STRATEGY")->[0]->textContent;

                my $attributes = $xml->getElementsByLocalName("EXPERIMENT_ATTRIBUTE");    # get the mobyData block      

                foreach my $attr (@{$attributes}){

                        my @nodes = $attr->getChildrenByTagName("TAG");
                        my $tag = $nodes[0]->textContent;

                        @nodes = $attr->getChildrenByTagName("VALUE");
                        my $value = $nodes[0]->textContent;

                        if ($value eq "NA" || $value eq "unknown"){
                                $value = "#VO#";
                        }
#  print $tag . "--" . $value . "\n";

                        $ihec->{$tag} = $value;
                }
        } else {
                warn "NO EXISTE $sfile\n";
        }


	my $analysis_id = "";
	my $analysis_group_id = "";
    my $assembly_version = "";
    my $program_versions = "";
	my $alignment_algorithm = "";
    my $other_analysis_algorithm = "";

	my $type = $arr[EXP_TYPE];
	if ($type eq "RNA_GENE_QUANT_CRG"){

		$analysis_id = $arr[AN_ID] . ".exp_CRG";
		$analysis_group_id = "15b"; # Centre for Genomic Regulation (R. Guigo)
		$assembly_version = "1"; #GRCh37
		$program_versions = "#VO#";
		$alignment_algorithm = "GEMTools library;http://github.com/gemtools";
		$other_analysis_algorithm = "Samtools;http://samtools.sourceforge.net,Cufflinks;http://cufflinks.cbcb.umd.edu/manual.html#cufflinks_input,Flux Capacitor;http://sammeth.net/confluence/display/FLUX/Home";

	} elsif ($type eq "RNA_GENE_QUANT_CBR"){

		# release 201403
		$analysis_id = $arr[AN_ID] . ".exp_CBR";
        $analysis_group_id = "3b"; # Department of Haematology, University of Cambridge (W Ouwehand)
 		$assembly_version = "8"; # Grch37.P10 
		$program_versions = "#VO#";
		$alignment_algorithm = "Bowtie 0.12.8;http://bowtie-bio.sourceforge.net/index.shtml"; # Bowtie 0.12.8
		$other_analysis_algorithm = "MMSEQ 1.0.5;https://github.com/eturro/mmseq";
    }

	# release 201403
	my $ensembl_version = "70"; # Ensembl70
    my $gencode_version = "15"; # Gencode v15

    my $data_status = "2"; #Raw data and results available
	my $seq_coverage = "#VO#";
    my $raw_data_repository = "1";#EGA
	my $raw_data_accession = "#VO#";
	my $platform = "60";#Illumina HiSeq
	my $notes = "#VO#";

	my $analysis_group_id_term = "Centre for Genomic Regulation (R. Guigo)";
    my $data_status_term = "Raw data and results available";
    my $assembly_version_term = "GRCh37";
    my $raw_data_repository_term = "EGA";
    my $platform_term = "Illumina HiSeq";

	#TODO : modificar
	$str2 = "$analysis_id\t$experiment_id\t$analysis_group_id\t$data_status\t$assembly_version\t$ensembl_version\t$gencode_version\t$program_versions\t$alignment_algorithm\t$other_analysis_algorithm\t$notes\t$analysis_group_id_term\t$data_status_term\t$assembly_version_term\n";

	print EXP_M2 $str2;

        my $analyzed_sample_id = $arr[AN_SAMPLE_ID];
        my $experiment_type;

        my $experimental_group_id = "12d"; 

        my $features = "";
        my $extraction_protocol = "#VO#";
        my $experimental_group_id_term = "Max Planck Institute for Molecular Genetics, Dpt Vertebrate Genomics (M-L Yaspo)";

        foreach my $key (keys %{$ihec}){
        #HA DE SERUN HASH!!! print key=value;

                if ($key eq "EXTRACTION_PROTOCOL") {
                        $extraction_protocol = $ihec->{$key};
                } elsif ($key eq "EXPERIMENT_TYPE") {
                        $experiment_type = $ihec->{$key};
                } else {
                        $features .= ";" if ($features);
                        $features .= "$key=".$ihec->{$key};
                }

#print "$key => ".$ihec->{$key}."\n";

        }

		$features =~ s/\n/,/g;
		$features =~ s/\t/ /g;

		#TODO : modificar 
		$str2 = "$experiment_id\t$analyzed_sample_id\t$experiment_type\t$library_strategy\t$experimental_group_id\t$features\t$raw_data_repository\t$raw_data_accession\t$platform\t$seq_coverage\t$extraction_protocol\t$experimental_group_id_term\t$raw_data_repository_term\t$platform_term\n";

        print LAB2 $str2;
}

close(F);


`cat exp_m.txt.all |sort -u >> exp_m.txt`;


`cat lab_mrna.txt.all | sort -u >> lab_mrna.txt`;



