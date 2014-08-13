#!/usr/bin/perl

use strict;
use constants;
use XML::LibXML;
use utf8;

binmode(STDOUT, ":utf8");

my $file = shift;
my $release = shift;

my $cmd;

my $path = "/data2/INB/blueprint/ftp/";
my $path_IHEC = $path . $release . "/blueprint/metadata_dump/experiments";

open(RREG_M,">rreg_m.txt");
open(RREG_M2,">rreg_m.txt.all");
binmode(RREG_M2,":utf8");

open(LAB,">lab_chro.txt");
open(LAB2,">lab_chro.txt.all");
binmode(LAB2, ":utf8");

#TODO : modificar segons el model (part 1)
my $gene_m_str="analysis_id\texperiment_id\tanalysis_group_id\tdata_status\tassembly_version\tensembl_version\tgencode_version\tprogram_versions\talignment_algorithm\tother_analysis_algorithm\tnotes\tanalysis_group_id_term\tdata_status_term\tassembly_version_term\n";

#TODO : modificar segons el model (part 2)
my $lab_chro_str = "experiment_id\tanalyzed_sample_id\texperiment_type\tlibrary_strategy\texperimental_group_id\tfeatures\traw_data_repository\traw_data_accession\tplatform\tseq_coverage\textraction_protocol\texperimental_group_id_term\traw_data_repository_term\tplatform_term\n";

print RREG_M $gene_m_str;

print LAB $lab_chro_str;


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
                        #print $tag . "--" . $value . "\n";

                        $ihec->{$tag} = $value;
                }
        } else {
                warn "NO EXISTE $sfile\n";
        }


    my $analysis_id = $arr[AN_ID] . ".ds_hotspot";
    my $analysis_group_id = "8";#European Bioinformatics Institute (P. Flicek)
    my $data_status = "2";#Raw data and results available
    my $assembly_version = "1";#GRCh37

    # release 20140
    my $ensembl_version = "70"; # Ensembl70
    my $gencode_version = "15"; # Gencode v15

   	my $analyzed_sample_id = $arr[AN_SAMPLE_ID];

	my $alignment_algorithm = "BWA;http://bio-bwa.sourceforge.net/";
	my $program_versions = "BWA:0.5.9;samtools:0.1.18;phantompeakqualtools:1.1;spp:1.11;macs2:2.0.10.20120913"; #"BWA:0.5.9;Hotspot:v3;samtools:0.1.18;Bedops:1.2.5";
	my $other_analysis_algorithm = "samtools:http://samtools.sourceforge.net/;Hotspot:http://www.uwencode.org/proj/hotspot-ptih/;Bedops:http://code.google.com/p/bedops/";
	my $assembly_version_term = "GRCh37";	
	my $data_status_term = "Raw data and results available";
	my $analysis_group_id_term = "European Bioinformatics Institute (P. Flicek)";
	my $control_sample_id = "#VO#";
	my $notes = "#VO#";
	my $seq_coverage = "#VO#";
        my $raw_data_repository = "1";#EGA
	my $raw_data_repository_term = "EGA";
        my $raw_data_accession = "#VO#";
        my $platform = "60";#Illumina HiSeq
	my $platform_term = "Illumina HiSeq";


	#TODO : modificar segons el model (igual que - part 1 - )
	$str2 = "$analysis_id\t$experiment_id\t$analysis_group_id\t$data_status\t$assembly_version\t$ensembl_version\t$gencode_version\t$program_versions\t$alignment_algorithm\t$other_analysis_algorithm\t$notes\t$analysis_group_id_term\t$data_status_term\t$assembly_version_term\n";

	print RREG_M2 $str2;

        my $experiment_type;

        my $experimental_group_id = "1"; #"Radboud University Nijmegen (H. Stunnenberg)";

####rnaseq      my $experimental_group_id = "12d"; #Max Planck Institute for Molecular Genetics, Dpt Vertebrate Genomics (M-L Yaspo)

        my $features = "";
        my $extraction_protocol = "#VO#";
        my $experimental_group_id_term = "Radboud University Nijmegen (H. Stunnenberg)";

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
        }

        #print $features . "\n";

		#TODO : modificar segons el model (igual que - part 2 - )
		my $str = "$experiment_id\t$analyzed_sample_id\t$experiment_type\t$library_strategy\t$experimental_group_id\t$features\t$raw_data_repository\t$raw_data_accession\t$platform\t$seq_coverage\t$extraction_protocol\t$experimental_group_id_term\t$raw_data_repository_term\t$platform_term\n";

        print LAB2 $str;

}

close(F);


`cat rreg_m.txt.all |sort -u >> rreg_m.txt`;


`cat lab_chro.txt.all | sort -u >> lab_chro.txt`;


	
