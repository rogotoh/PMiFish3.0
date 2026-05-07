#!/usr/bin/perl
use strict;
use warnings;

#PMiFish ver.3.0

#Setting
my ($algo, $denoise, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Algorithm\s*=\s*(\d+)/){$algo = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^Type\s*=\s*(\S+)/){$type = $1;}
}
close(SET);
unless($denoise){$denoise = "yes";}
unless($algo){$algo = 1;}
unless($type){$type = 1;}

if($denoise =~ /^yes$/i){
	if($algo == 1){
	#usearch/vseasrch
		system "perl ./Scripts/1_1_Preprocessing.pl";
		system "perl ./Scripts/2_1_Find_unique_in_each_sample.pl";
		if($type == 2){
			system "perl ./Scripts/2_4_0_Uniques.pl";
			system "perl ./Scripts/2_4_1_Pool_Denoise.pl";
			system "perl ./Scripts/2_4_2_Separate_chimera.pl";
		}else{
			system "perl ./Scripts/2_2_Denoise.pl";
			system "perl ./Scripts/2_3_Separate_chimera.pl";
			system "perl ./Scripts/2_4_ASVs.pl";
		}
		system "perl ./Scripts/2_5_Rarefaction.pl";
		system "perl ./Scripts/2_6_OTU_Clustering.pl";
		system "perl ./Scripts/2_7_LULU.pl";
		system "perl ./Scripts/3_1_Usearch_global.pl";
		system "perl ./Scripts/4_1_Annotation.pl";
		system "perl ./Scripts/5_1_Fasta_for_Phylogenetic_Analysis.pl";
		system "perl ./Scripts/5_2_Summary_Table.pl";
		system "perl ./Scripts/5_3_Fasta_classified_by_family.pl";
		system "perl ./Scripts/6_1_Portal.pl";
	}else{
	#DADA2
		system "perl ./Scripts/D_1_Trim_Primer_for_DADA2.pl";
		system "perl ./Scripts/D_2_DADA2.pl";
		system "perl ./Scripts/D_3_Table_converter.pl";
		system "perl ./Scripts/2_5_Rarefaction.pl";
		system "perl ./Scripts/2_6_OTU_Clustering.pl";
		system "perl ./Scripts/2_7_LULU.pl";
		system "perl ./Scripts/3_1_Usearch_global.pl";
		system "perl ./Scripts/4_1_Annotation.pl";
		system "perl ./Scripts/5_1_Fasta_for_Phylogenetic_Analysis.pl";
		system "perl ./Scripts/5_2_Summary_Table.pl";
		system "perl ./Scripts/5_3_Fasta_classified_by_family.pl";
		system "perl ./Scripts/6_1_Portal.pl";
	}
}else{
	system "perl ./Scripts/1_1_Preprocessing.pl";
	system "perl ./Scripts/2_1_Find_unique_in_each_sample.pl";
	system "perl ./Scripts/2_2_Denoise.pl";
	system "perl ./Scripts/2_3_Separate_chimera.pl";
	system "perl ./Scripts/2_4_ASVs.pl";
	system "perl ./Scripts/2_5_Rarefaction.pl";
	system "perl ./Scripts/2_6_OTU_Clustering.pl";
	system "perl ./Scripts/2_7_LULU.pl";
	system "perl ./Scripts/3_1_Usearch_global.pl";
	system "perl ./Scripts/4_1_Annotation.pl";
	system "perl ./Scripts/5_1_Fasta_for_Phylogenetic_Analysis.pl";
	system "perl ./Scripts/5_2_Summary_Table.pl";
	system "perl ./Scripts/5_3_Fasta_classified_by_family.pl";
	system "perl ./Scripts/6_1_Portal.pl";
}
