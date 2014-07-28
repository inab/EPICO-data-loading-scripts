#!/ussdata_donor.txt.allr/bin/perl

use constants;
use strict;

my $file = shift;
my $release = shift;

my $cmd;

my $path = "/data2/INB/blueprint/ftp/" . $release;



#PDNA

open(RREGP,">rreg_p.txt");

#TODO : modificar segons el model (part 1)
my $rreg_p_str = "analysis_id\tregulatory_region_id\tchromosome\tchromosome_start\tchromosome_end\tz_score\tnotes\n";

#print PDNAM $pdna_m_str;
print RREGP $rreg_p_str;

open(F,$file);

while(<F>){
	my @arr = split('\t',$_);


	my $analysis_id = $arr[AN_ID] .".ds_hotspot";

	my $exp_type = $arr[EXP_TYPE];
	my $name = $arr[SAMPLE_NAME];

	my $res = "$path/".$arr[FILE];	
	print "$res\n";
	
	open(IN, "gunzip -c $res |") || die "can't open pipe to $res";
	while (<IN>) {
        	if (/^#/ || /^chr\s/ || /^\s/){
                	next;
        	}
		my @line = split('\t',$_);
		my $str2;
		my $chromosome = $line[CHR];
		$chromosome = substr($chromosome,3,length($chromosome)-3);
		my $chromosome_start = $line[START];
		my $chromosome_end = $line[END];

		my $validation_status ="Not tested";
		my $validation_platform = "#VO#";
		my $notes = "#VO#";

		my $z_score = $line[ZSCORE_PEAKS];
	
		my $regulatory_region_id = $name."_20130228_".$chromosome."_".$line[END]."_".$line[START];

		chomp($z_score);

		#TODO : modificar segons el model (igual que - part 1 - )
		$str2 = "$analysis_id\t$regulatory_region_id\t$chromosome\t$chromosome_start\t$chromosome_end\t$z_score\t$notes\n";
		print RREGP $str2;
	}
	close(IN);
}

close(F);






	
