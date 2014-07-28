#!/ussdata_donor.txt.allr/bin/perl

use constants;
use strict;

my $file = shift;
my $release = shift;

my $cmd;

my $path = "/data2/INB/blueprint/ftp/$release"; #my $path = "/data2/INB/users/rroyo/roadmap/releases/" . $release;


open(PDNAP,">pdna_p.txt");

#TODO: modificar segons model
my $pdna_p_str="analysis_id\tprotein_dna_interaction_id\tchromosome\tchromosome_start\tchromosome_end\trank\tprotein_stable_id\tlog10_pvalue\tlog10_qvalue\tidr\tnotes\n";

#print PDNAM $pdna_m_str;
print PDNAP $pdna_p_str;

open(F,$file);

while(<F>){
	my @arr = split('\t',$_);


	my $analysis_id = $arr[AN_ID] .".pdna";

	my $exp_type = $arr[EXP_TYPE];

	my $res = "$path/".$arr[FILE];	
	print "$res\n";
	
	open(IN, "gunzip -c $res |") || die "can't open pipe to $res";
	my $num = 1;
	while (<IN>) {
        	if (/^#/ || /^chr\s/ || /^\s/){
                	next;
        	}
		my @line = split('\t',$_);
		my $str2;
		my $chromosome = $line[CHR];
		$chromosome = substr($chromosome,3,length($chromosome)-3);
		if ($chromosome eq "M"){
			$chromosome = "MT";
		}

		my $chromosome_start = $line[START];
		my $chromosome_end = $line[END];
		my $rank = "#VO#";
		my $protein_stable_id = $arr[PROT_STABLE_ID];
		my $log10_pvalue = "";
		my $log10_qvalue = "";
		my $idr = "#VO#";
		my $notes = "#VO#";
	
		if ($exp_type eq "CS_MACS"){
			my $protein_dna_interaction_id = $line[PROT_INT_ID];
			chomp($protein_dna_interaction_id);
			$log10_pvalue = $line[PVALUE];
			$log10_qvalue = $line[QVALUE];

			#TODO: modificar segons model
			$str2 = "$analysis_id\t$protein_dna_interaction_id\t$chromosome\t$chromosome_start\t$chromosome_end\t$rank\t$protein_stable_id\t$log10_pvalue\t$log10_qvalue\t$idr\t$notes\n";

			print PDNAP $str2;

		} elsif ($exp_type eq "CS_BROAD_MACS"){
			my $sample_name = $arr[6];
			#bwa_filtered.20130401_macs2_peak_1
			my $protein_dna_interaction_id = $sample_name.".".$arr[PROT_STABLE_ID].".bwa_filtered.20130401_macs2_broad_peak_".$num;
			$num++;
			$log10_pvalue = $line[PVALUE_BROAD];
			$log10_qvalue = "#VO#";

			#TODO: modificar segons model
			$str2= "$analysis_id\t$protein_dna_interaction_id\t$chromosome\t$chromosome_start\t$chromosome_end\t$rank\t$protein_stable_id\t$log10_pvalue\t$log10_qvalue\t$idr\t$notes\n";

			print PDNAP $str2;
	        } else {
        	}
	}
	$num = 1;
	close(IN);
}

close(F);






	
