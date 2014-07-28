#!/usr/bin/perl


use constants;
use strict;


my $file = shift;
my $release = shift;

my $path = "/data2/INB/blueprint/ftp/$release"; # my $path = "/data2/INB/users/rroyo/roadmap/releases/" . $release;


my $cmd;

open(EXP_G,">exp_g.txt");

#TODO : modificar segons el model (part 1)
my $gene_str = "analysis_id\tgene_stable_id\tchromosome\tchromosome_start\tchromosome_end\tchromosome_strand\tquality_score\tprobability\tnormalized_read_count\traw_read_count\tnormalized_expression_level\treference_sample\tis_annotated\tprobeset_id\tnotes\n";


print EXP_G $gene_str;

open(F,$file);

while(<F>){
	my @arr = split('\t',$_);

	my $res = $path ."/". $arr[FILE];	


	print "$res\n";
	my $type = $arr[EXP_TYPE];
	my $analysis_id;

	if ($type eq "RNA_GENE_QUANT_CRG"){
		open(IN, $res) || die "can't open pipe to $res";
		$analysis_id = $arr[AN_ID].".exp_CRG";
	} elsif ($type eq "RNA_GENE_QUANT_CBR"){
		open(IN, "gunzip -c $res |") || die "can't open pipe to $res";	
		$analysis_id = $arr[AN_ID].".exp_CBR";
	}
	while (<IN>) {

		my $str2;
		my @line = split('\t',$_);
		
		my $gene_stable_id = "";

		my $chromosome_start = $line[START];
		my $chromosome_end = $line[END];
		
		my $quality_score = "#VO#";
		my $probability = "#VO#";
		my $normalized_read_count = "#VO#";
		my $raw_read_count = "";
		my $normalized_expression_level = "#VO#";
		my $reference_sample = "#VO#";
		my $is_annotated = "true";
		my $probeset_id = "#VO#";
		my $validation_status = "Not tested";
		my $validation_platform = "#VO#";
		my $notes = "#VO#";

		my $chromosome = substr($line[CHR],3,length($line[CHR]));
		if ($chromosome eq "M"){
			$chromosome = "MT";
		}
	
		my $chromosome_strand = $line[STRAND];
		if ($chromosome_strand eq '+'){
			$chromosome_strand = 1;
		} elsif ($chromosome_strand eq '-'){
			$chromosome_strand = -1;
		} else {
			print "STRAND $chromosome_strand\n";
			$chromosome_strand = "#VO#";
		}

		my $info = $line[GENE];
		if ($info =~ /gene_id\s+"(.+)";\s+transcript_ids\s+".+";\s+RPKM\s+(.+);\s+reads\s+(.+);/){

			# Parsegem RNA_GENE_QUANT_CRG
			$gene_stable_id = $1;
			$normalized_read_count = $2;
			$raw_read_count = $3;

		} elsif($info =~ /RPKM\s+["]?(.+)["]?\s*;\s*gene_id\s+"(.+)"\s*;.*reads\s+(.+)\s*;.*/) { 

#INFO RPKM 0 ; gene_id "ENSG00000243485.1"  ; length 1247 ; log_mu "-8.74953"  ; reads 0 ; transcript_id "ENSG00000243485.1" 

			$gene_stable_id = $2;
            $normalized_read_count = $1;
			$normalized_read_count =~ s/"//g;
            $raw_read_count = $3;

			#RPKM "0.475452994804434"  ; gene_id "ENSG00000227232.3"  ; length 6707 ; log_mu "1.46481"  ; reads 346 ; transcript_id "ENSG00000227232.3"

			# nous fitxers gff
			# Parsegem RNA_GENE_QUANT_CBR


		} else {
			print "INFO $info\n";
		}


		my @genename = split('\.',$gene_stable_id);
		$gene_stable_id = $genename[0];
#print $gene_id."---\n";;

		#TODO : modificar segons el model (igual que - part 1 - )
		$str2 = "$analysis_id\t$gene_stable_id\t$chromosome\t$chromosome_start\t$chromosome_end\t$chromosome_strand\t$quality_score\t$probability\t$normalized_read_count\t$raw_read_count\t$normalized_expression_level\t$reference_sample\t$is_annotated\t$probeset_id\t$notes\n";

		print EXP_G $str2;
	}
	close(IN);
}

close(F);






	
