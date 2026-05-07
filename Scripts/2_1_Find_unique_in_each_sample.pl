#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($temporary);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
}
close(SET);
unless($temporary){$temporary = "YES";}

#Usearch options
my ($minuniquesize);
open (SET, "<", "Options_usearch.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^\t*-minuniquesize\s*(\d+)/){$minuniquesize = $1;}
}
close(SET);
unless($minuniquesize){$minuniquesize = 1;}

#Usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1; last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}


#Get data names
opendir (DIR, ".\/Results\/1_3_Quality_filter") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_filtered.fa/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "              2_1_Find_unique_in_each_sample                \n";
print "============================================================\n";

#Find_unique_reads
mkdir ".\/Results\/2_1_Find_unique_in_each_sample";
foreach(sort keys %file){
	print "\n\n$_\n";
	my $file = $_;
	my $command;
	if($vv){
		$command = ".\/Tools\/$usearch --fastx_uniques \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" --sizeout --relabel Uniq --fastaout \".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa\" --minuniquesize $minuniquesize";
	}else{
		$command = ".\/Tools\/$usearch -fastx_uniques \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -sizeout -relabel Uniq -fastaout \".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa\" -minuniquesize $minuniquesize";
	}
	system $command;
	unless(-f ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa"){
		open (DATA, ">", ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa") or die("error:$!");
		close(DATA);
	}
	if($temporary =~ /yes/i){
		unlink ".\/Results\/1_3_Quality_filter\/${file}_filtered.fa";
	}
}
