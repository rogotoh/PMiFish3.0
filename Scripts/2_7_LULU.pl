#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

#2025/01/28

#Setting
my ($cluster, $identity, $denoise, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Curation\s*=\s*(\S+)/){$cluster = $1;}
	elsif($_ =~ /^pcid\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
}
close(SET);
unless($cluster){$cluster = "no";}
if($cluster =~ /^no$/i){exit;}
unless($identity){print "Error: Please check the Post-Clustering Identity(pcid) value in Setting.txt.\n"; exit;}
unless($denoise){$denoise = "no";}

my ($usearch, $vv);
#userch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}


my $as_us;
unless($denoise =~ /^yes$/i){
	$as_us = "Unique";
}else{
	if($type == 1){
		$as_us = "ZOTU";
	}else{
		$as_us = "ASV";
	}
}

#data input
my ($data, $dseq, %data, $filename, $tablename);
if(-f ".\/Results\/2_6_OTU_Clustering\/OTU.fa"){
	$filename = ".\/Results\/2_6_OTU_Clustering\/OTU.fa";
	$tablename = ".\/Results\/2_6_OTU_Clustering\/OTU_table.tsv";
}elsif(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	$filename = ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa";
	$tablename = ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv";
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	$filename = ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa";
	$tablename = ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv";
}else{
	exit;
}
open (DATA, "<", $filename) or die("error:$!");
mkdir ".\/Results\/2_7_LULU";
open (OUT, ">", ".\/Results\/2_7_LULU/temp.fa") or die("error:$!");

while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	unless($_ =~ /^\n/){
		if($_ =~ /^>(.+)/){
			$data = $1;
			$_ =~ s/\;size.+//;
			print OUT "$_\n";
		}else{
			$data{$data} = $_;
			print OUT "$_\n";
		}
	}
}
close(DATA);
close(OUT);

print "============================================================\n";
print "                         2_7_LULU                           \n";
print "============================================================\n";

my $command;
if($vv){
	$command = ".\/Tools\/$usearch --usearch_global \".\/Results\/2_7_LULU/temp.fa\" --db \".\/Results\/2_7_LULU/temp.fa\" --self --id $identity --iddef 1 --maxaccepts 0 --query_cov 0.9 --userout \".\/Results\/2_7_LULU/match_list.txt\" --userfields query+target+id";
}else{
	$command = ".\/Tools\/$usearch -usearch_global \".\/Results\/2_7_LULU/temp.fa\" -db \".\/Results\/2_7_LULU/temp.fa\" -self -id $identity -maxaccepts 0 -query_cov 0.9 -userout \".\/Results\/2_7_LULU/match_list.txt\" -userfields query+target+id -strand plus";
}
system $command;

open (OUT, ">", ".\/Results\/2_7_LULU\/LULU.R") or die("error:$!");

print OUT <<"EOS";
library(lulu)

otutab = read.table("$tablename", row.names=1)
matchlist = read.table("./Results/2_7_LULU/match_list.txt", header=FALSE, as.is=TRUE)
curated_result = lulu(otutab, matchlist)
write.csv(curated_result\$curated_table, ".\/Results\/2_7_LULU\/temp.csv")
write.csv(curated_result\$otu_map, ".\/Results\/2_7_LULU\/otu_map.csv")

EOS

$command = "Rscript .\/Results\/2_7_LULU\/LULU.R";
system $command;
close(OUT);

#Table
open (DATA, "<", ".\/Results\/2_7_LULU\/temp.csv") or die("error:$!");
open (OUT, ">", ".\/Results\/2_7_LULU\/curated_list.tsv") or die("error:$!");
my $count = 0;
my (%narabi, %check, %size);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	$_ =~ s/\"//g;
	$_ =~ s/,/\t/g;
	if($count == 0){
		$_ =~ s/X\d+\.//g;
		print OUT "$_\n";
		$count++;
	}else{
		$_ =~ /(\d+)/;
		$narabi{$1} = $_;
		$_ =~ /^([^\t]+)/;
		my $tempname = $1;
		$check{$1}++;
		$size{$tempname} = 0;
		while($_ =~ /\t(\d+)/g){$size{$tempname} += $1;}
	}
}
close(DATA);
foreach(sort {$a <=> $b} keys %narabi){print OUT "$narabi{$_}\n";}
close(OUT);

unlink(".\/Results\/2_7_LULU\/temp.csv");

#fasta
open (DATA, "<", ".\/Results\/2_7_LULU\/temp.fa") or die("error:$!");
open (OUT, ">", ".\/Results\/2_7_LULU\/curated.fa") or die("error:$!");
$count = 0;
my $curated = 0;
my $discarded = 0;
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>(.+)/){
		if($check{$1}){
			print OUT "$_\;size=$size{$1}\;\n";
			$count = 1;$curated++;
		}else{
			$discarded++;
			$count = 0;
		}
	}else{
		if($count == 1){print OUT "$_\n";}
	}
}
print "\nCurated = $curated Discarded = $discarded\n";
close(DATA);
close(OUT);
unlink(".\/Results\/2_7_LULU\/temp.fa");

my $dest_dir = "./Results/2_7_LULU/";
foreach my $file (glob("lulu.log*")) {
    move($file, $dest_dir) or die ("error:$!");
}
