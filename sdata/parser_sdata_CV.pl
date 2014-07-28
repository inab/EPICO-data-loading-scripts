#!/ussdata_donor.txt.allr/bin/perl
use strict;
use XML::LibXML;
use constants;

my $file = shift;
my $release = shift;

my $cmd;

#TODO: modificar el $path
my $path = "/data2/INB/blueprint/ftp/"; # my $path = "/data2/INB/users/rroyo/roadmap/releases/";
my $path_IHEC = $path . $release . "/blueprint/metadata_dump/samples";


#DONOR
#TODO : modificar segons el model (part 1)
my $donor_str = "donor_id\tdonor_sex\tdonor_region_of_residence\tdonor_ethnicity\tnotes\tdonor_sex_term\tdonor_region_of_residence_term\n";

open(DON,">sdata_donor.txt");
print DON $donor_str;
open(DON2,">sdata_donor.txt.all");


#SPECIMEN
#TODO : modificar segons el model (part 2)
my $spec_str = "specimen_id\ttissue_type\ttissue_depot\tspecimen_term\tcollection_method\tdonor_min_age_at_specimen_acquisition\tdonor_max_age_at_specimen_acquisition\tdonor_health_status\tdonor_disease\tdonor_disease_text\tspecimen_processing\tspecimen_processing_other\tspecimen_storage\tspecimen_storage_other\tspecimen_biomaterial_provider\tspecimen_biomaterial_id\tspecimen_available\tdonor_id\tnotes\tspecimen_term_term\tspecimen_processing_term\tspecimen_storage_term\n";

open(SPE,">sdata_specimen.txt");
print SPE $spec_str;
open(SPE2,">sdata_specimen.txt.all");


#SAMPLE
#TODO : modificar segons el model (part 3)
my $sample_str ="sample_id\tpurified_cell_type\tculture_conditions\tmarkers\tanalyzed_sample_type\tanalyzed_sample_type_other\tanalyzed_sample_interval\tspecimen_id\tnotes\tpurified_cell_type_term\tanalyzed_sample_type_term\n";

open(SAM,">sdata_sample.txt");
print SAM $sample_str;
open(SAM2,">sdata_sample.txt.all");


open(F,$file);

while(<F>){
	my @arr = split('\t',$_);

    my $sample_id = $arr[AN_SAMPLE_ID];

    my $sfile = $path_IHEC . "/" . substr($sample_id,0,6) . "/" . $sample_id . ".xml";

    my $str_xml = `cat $sfile`;
	my $ihec;

    my $parser = XML::LibXML->new();
    my $doc;

	#TODO: afegir el nou atribut library_strategy
	my $library_strategy;

    eval { $doc = $parser->parse_string( $str_xml ) };
    die("couldn't parse XML '$str_xml' because\n\t$@") if $@;

    my $xml = $doc->getDocumentElement();

	#TODO: afegir el nou atribut library_strategy
    #### mon : no ho te!! :: $library_strategy = $xml->getElementsByLocalName("LIBRARY_STRATEGY")->[0]->textContent;

    my $attributes = $xml->getElementsByLocalName("SAMPLE_ATTRIBUTE");    # get the mobyData block      

    foreach my $attr (@{$attributes}){

        #print $attr->toString() . "\n";
        no strict 'refs';

        my @nodes = $attr->getChildrenByTagName("TAG");
        my $tag = $nodes[0]->textContent;

        @nodes = $attr->getChildrenByTagName("VALUE");
        my $value = $nodes[0]->textContent;

        if ($value eq "NA" || $value eq "unknown" || $value eq "Unknown" || $value eq "-"){
        	$value = "#VO#";
        }

    	$ihec->{$tag} = $value;
    }


	### en queden de repetits!!
	my $donor_sex = $arr[SEX];
    my $donor_sex_term;

	# mon
	if ($donor_sex eq "Unknown" || $donor_sex eq "-"){
		$donor_sex = "#VO#";
		$donor_sex_term = "#VO#";
	}

	if ($donor_sex eq "Female"){
		$donor_sex = "f";
		$donor_sex_term = "female";
	} elsif ($donor_sex eq "Male"){
		$donor_sex = "m";
		$donor_sex_term = "male";
	}

	my $donor_region_of_residence = $arr[REGION];
	my $donor_region_of_residence_term;

	# mon
    if ($donor_region_of_residence eq "Unknown" || $donor_region_of_residence eq "-"){
        $donor_region_of_residence = "#VO#";
		$donor_region_of_residence_term = "#VO#";
    }

	#my $donor_region_of_residence_term;
	if ($donor_region_of_residence eq "East Anglia"){
		$donor_region_of_residence = "ALIAS:EAL";
		$donor_region_of_residence_term = "East Anglia: United Kingdom region composed of the administrative counties of Norfolk to the north, Suffolk to the south, Cambridgeshire and Essex to the west.";
	}
	my $donor_id = $arr[DONOR_ID];
	my $donor_ethnicity = $ihec->{DONOR_ETHNICITY} || "#VO#";
	my $notes = "#VO#";
	
	#TODO : modificar segons el model (igual que - part 1 - )
	my $str = "$donor_id\t$donor_sex\t$donor_region_of_residence\t$donor_ethnicity\t$notes\t$donor_sex_term\t$donor_region_of_residence_term\n";

	print DON2 $str;

	my $donor_min_age_at_specimen_acquisition;
	my $donor_max_age_at_specimen_acquisition;
	if ($arr[AGE] =~ /(\d+)\s+-\s+(\d+)/){
		$donor_min_age_at_specimen_acquisition = $1;
		$donor_max_age_at_specimen_acquisition = $2;
	} else {
		$donor_min_age_at_specimen_acquisition = "#VO#";
		$donor_max_age_at_specimen_acquisition = "#VO#";
	}
	my $specimen_type;
    my $specimen_type_term;
	my $speciment_type_aux;

	$speciment_type_aux = $ihec->{SPECIMEN_TYPE};

    my $tissue_type = $ihec->{TISSUE_TYPE} || "#VO#";
	if($tissue_type eq "#VO#") { $tissue_type = $arr[18]; } # Agafem la columna TISSUE del public.results.index

	my $tissue_depot = $ihec->{TISSUE_DEPOT} || "#VO#";

	my $specimen_term = "#VO#";
	my $specimen_term_term = "#VO#";

    if ($tissue_type eq "Peripheral blood"){
		$specimen_term = "UBERON:0013756";
		$specimen_term_term = "Venous blood";

    } elsif ($tissue_type eq "Cord blood"){
		$specimen_term = "UBERON:0012168";
        $specimen_term_term = "umbilical cord blood";

    } elsif ($tissue_type eq "Tonsil"){
		$specimen_term = "UBERON:0002372";
        $specimen_term_term = "tonsil";

    } elsif ($tissue_type eq "Bone marrow"){
		$specimen_term = "UBERON:0002371";
        $specimen_term_term = "bone marrow";

    } else {
       print $arr[FILE] ."++++++++++++++++++++++++++++++++++++++\n";
    }

	#if ($arr[FILE] =~ /Peripheral_blood/){
	#	$specimen_type = "per_blood";
	#	$specimen_type_term = "peripheral blood";
	#} elsif ($arr[FILE] =~ /Cord_blood/){
	#	$specimen_type = "cord_blood";
	#	$specimen_type_term = "cord blood";
	#} elsif ($arr[FILE] =~ /Tonsil/){
    #    $specimen_type = "c_other";
    #    $specimen_type_term = "normal control (other)";
    #} elsif ($arr[FILE] =~ /Cell_Line/){
    #    $specimen_type = "c_other";
    #    $specimen_type_term = "normal control (other)";
    #} elsif ($arr[FILE] =~ /Bone_marrow/){
    #    $specimen_type = "6";
    #    $specimen_type_term = "bone marrow";
    #} else {
	#	print $arr[FILE] ."++++++++++++++++++++++++++++++++++++++\n";
	#}
	#my $specimen_type_other = "#VO#";
	#if ($arr[FILE] =~ /Tonsil/){
	#	$specimen_type_other = "tonsil"; # mon
	#} elsif ($arr[FILE] =~ /Cell_Line/) {
	#	 $specimen_type_other = "cell line"; # mon
	#}

	my $specimen_processing_other = "#VO#";
	my $specimen_source = "#VO#";
	my $specimen_storage_other = "#VO#";
	#my $specimen_biobank = $ihec->{BIOMATERIAL_PROVIDER}; 	# mon
	#my $specimen_biobank_id = "#VO#"; 						# mon
	my $specimen_biomaterial_provider = $ihec->{BIOMATERIAL_PROVIDER};
	my $specimen_biomaterial_id = "#VO#";
	my $specimen_available = "#VO#";

	my $specimen_processing = $arr[SPEC_PROC];
    my $specimen_processing_term;
	if ($specimen_processing eq "fresh"){
		$specimen_processing = 8;
		$specimen_processing_term = "fresh";
	} elsif ($specimen_processing eq "-"){
		$specimen_processing = "#VO#";
		$specimen_processing_term = "#VO#";
	}
	
	my $specimen_storage = "7";
	my $specimen_id = $donor_id."_spec";

	my $collection_method = $ihec->{COLLECTION_METHOD} || "#VO#";
#	my $specimen_collection_protocol = $ihec->{COLLECTION_METHOD} || "#VO#";

	my $donor_health_status = $ihec->{DONOR_HEALTH_STATUS};

	if($donor_health_status eq "") {
		$donor_health_status = $arr[DONOR_HEALTH_STATUS];
		#DONOR_HEALTH---columna 22 (mirar mail) de public.results.index];
	}


	my $donor_disease = "#VO#";

	if($donor_health_status eq "") {
		$donor_health_status = "#VO#";
		$donor_disease = "#VO#";
	} elsif($donor_health_status eq "Healthy") {
        $donor_disease = "EFO_0000761";
    } elsif($donor_health_status eq "Multiple myeloma") {
        $donor_disease = "EFO_0001378";
    } elsif($donor_health_status eq "Presumed healthy") {
        $donor_disease = "EFO_0000761 (Presumed)";
    } elsif($donor_health_status eq "Probable tonsillitis") {
        $donor_disease = "EFO_0000001 (Probable)";
    }


	my $donor_disease_text = $ihec->{DISEASE};

	my $specimen_storage_term = "other";


	#TODO : modificar segons el model (igual que - part 2 - )
	$str = "$specimen_id\t$tissue_type\t$tissue_depot\t$specimen_term\t$collection_method\t$donor_min_age_at_specimen_acquisition\t$donor_max_age_at_specimen_acquisition\t$donor_health_status\t$donor_disease\t$donor_disease_text\t$specimen_processing\t$specimen_processing_other\t$specimen_storage\t$specimen_storage_other\t$specimen_biomaterial_provider\t$specimen_biomaterial_id\t$specimen_available\t$donor_id\t$notes\t$specimen_term_term\t$specimen_processing_term\t$specimen_storage_term\n";

	print SPE2 $str;

	my $analyzed_sample_type = "n_blood";
	my $analyzed_sample_type_other = "#VO#";
	my $analyzed_sample_interval = "#VO#";

	my $purified_cell_type = `grep "$arr[CELL]" cell_types.txt | cut -f 2`;
	chomp($purified_cell_type);

	if (!$purified_cell_type){
		print "CELLTYPE ++++++++++++++++++ ".$arr[CELL]."\n";
	}
	
	my $culture_conditions = $ihec->{CULTURE_CONDITIONS} || "#VO#";
	my $markers = $ihec->{MARKERS} || "#VO#";

	my $purified_cell_type_term = $arr[CELL];
	my $analyzed_sample_type_term = "Normal blood";

	#TODO : modificar segons el model (igual que - part 3 - )
	$str = "$sample_id\t$purified_cell_type\t$culture_conditions\t$markers\t$analyzed_sample_type\t$analyzed_sample_type_other\t$analyzed_sample_interval\t$specimen_id\t$notes\t$purified_cell_type_term\t$analyzed_sample_type_term\n";

	print SAM2 $str;

}

close(F);

close(DON);
close(SPE);
close(SAM);
close(DON2);
close(SPE2);
close(SAM2);

`cat sdata_donor.txt.all |sort -u >> sdata_donor.txt`;
`cat sdata_specimen.txt.all |sort -u >> sdata_specimen.txt`;
`cat sdata_sample.txt.all |sort -u >> sdata_sample.txt`;




	
