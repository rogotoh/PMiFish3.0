#!/usr/bin/perl
use strict;
use warnings;

#2025/01/28

#Setting
my ($cluster, $identity, $depth, $denoise, $type, $type2, $dis);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Clustering\s*=\s*(\S+)/){$cluster = $1;}
	elsif($_ =~ /^Cluster_Identity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\d+)/){$depth = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Algorithm2\s*=\s*(\d+)/){$type2 = $1;}
	elsif($_ =~ /^Difference\s*=\s*(\d+)/){$dis = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
}
close(SET);
unless($cluster){$cluster = "no";}
if($cluster =~ /^no$/i){exit;}
unless($identity){print "Error: Please check the Cluster_Identity value in Setting.txt.\n"; exit;}
unless($denoise){$denoise = "no";}
unless($type){$type = 1;}
unless($type2){$type2 = 1;}
unless($dis){$dis = 1;}

#SWARM check
my ($sw, $usearch, $vv);
if($type2 == 2){
	opendir (DIR, ".\/Tools") or die ("error:$!");
	my @tool = readdir DIR;
	foreach (@tool) {
		if ($_ =~ /(swarm.*)/){$sw = $1;last;}
	}
	closedir DIR;
	unless($sw){print "Error: Please ensure the SWARM executable file is placed in the Tools directory.\n"; exit;}
}else{
	#userch check
	opendir (DIR, ".\/Tools") or die ("error:$!");
	my @tool = readdir DIR;
	foreach (@tool) {
		if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
		elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
	}
	closedir DIR;
	unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}
}

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
my ($data, $dseq, %data);
if(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	open (DATA, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa") or die("error:$!");
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa") or die("error:$!");
}else{
	exit;
}
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	unless($_ =~ /^\n/){
		if($_ =~ /^>(.+)/){
			$data = $1;
		}else{
			$data{$data} = $_;
		}
	}
}
close(DATA);

my $motoid = $identity;
$identity = $identity/100;

my (%count, $command);
if(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	print "============================================================\n";
	print "                    2_6_OTU_Clustering                      \n";
	print "============================================================\n";
	mkdir ".\/Results\/2_6_OTU_Clustering";
	if($sw){
		print "Clusterized by SWARM at $dis difference(s)\n";
		$command = ".\/Tools\/$sw \".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa\" -z -u \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" -d $dis ";
	}elsif($vv){
		print "Clusterized by VSEARCH at ${motoid}\% identity\n";
		$command = ".\/Tools\/$usearch --cluster_smallmem \".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa\" --id $identity --uc \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" --usersort";
	}else{
		print "Clusterized by USEARCH at ${motoid}\% identity\n";
		$command = ".\/Tools\/$usearch -cluster_smallmem \".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa\" -id $identity -uc \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" -sortedby size -quiet";
	}
	system $command;
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	print "============================================================\n";
	print "                     2_6_OTU_Clustering                     \n";
	print "============================================================\n";
	mkdir ".\/Results\/2_6_OTU_Clustering";
	if($sw){
		print "Clusterized by SWARM at $dis difference(s)\n";
		$command = ".\/Tools\/$sw \".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa\" -z -u \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" -d $dis";
	}elsif($vv){
		print "Clusterized by VSEARCH at ${motoid}\% identity\n";
		$command = ".\/Tools\/$usearch --cluster_smallmem \".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa\" --id $identity --uc \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" --usersort";
	}else{
		print "Clusterized by USEARCH at ${motoid}\% identity\n";
		$command = ".\/Tools\/$usearch -cluster_smallmem \".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa\" -id $identity -uc \".\/Results\/2_6_OTU_Clustering\/OTU.uc\" -sortedby size -quiet";
	}
	system $command;
}else{
	exit;
}
print "\n";

open (TEMP, "<", ".\/Results\/2_6_OTU_Clustering\/OTU.uc") or die("error:$!");
my (@cluster, %name);
while(<TEMP>){
	chomp($_);
	$_ =~ s/\r//g;
	my @temp = split(/\t/, $_);
	if($temp[0] eq "S"){
		push(@cluster, $temp[8]);
		$temp[8] =~ /size=(\d+)/;
		$count{$temp[8]} = $1;
		push(@{$name{$temp[8]}}, $temp[8]);
	}elsif($temp[0] eq "H"){
		$temp[8] =~ /size=(\d+)/;
		$count{$temp[9]} += $1;
		push(@{$name{$temp[9]}}, $temp[8]);
	}
}
close(TEMP);

#sort by size
my %narabi;
foreach(@cluster){
	push(@{$narabi{$count{$_}}}, $_);
}

open (OUT1, ">", ".\/Results\/2_6_OTU_Clustering\/OTU.fa") or die("error:$!");
open (OUT2, ">", ".\/Results\/2_6_OTU_Clustering\/OTU_list.txt") or die("error:$!");

#Get integrated sequence table
open (OUT3, ">", ".\/Results\/2_6_OTU_Clustering\/OTU_table.tsv") or die("error:$!");
if(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	open (DATA, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv") or die("error:$!");
}
my $retu = 0;
my (@file, %is);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($retu == 0){
		@file = split(/\t/, $_);
		print OUT3 "$_\n";
	}else{
		my @temp = split(/\t/, $_);
		my $length = @file;
		for(my $n = 1; $n < $length; $n++){
			$is{$temp[0]}{$file[$n]} = $temp[$n];
		}
	}
	$retu++;
}
close(DATA);
shift @file;

#output
my $num = 0;
foreach(sort {$b <=> $a} keys %narabi){
	my $narabi = $_;
	foreach(sort @{$narabi{$narabi}}){
		$num++;
		my $name = $_;
		my $temp_name = $name;
		$temp_name =~ s/${as_us}\d+//;
		$temp_name =~ s/size=\d+/size=$count{$name}/;
		print OUT1 ">OTU$num$temp_name\n$data{$name}\n";
		print OUT2 ">OTU$num$temp_name\n";
		print OUT3 "OTU$num";
		foreach(@{$name{$name}}){
			print OUT2 "\t$_\n";
		}
		foreach(@file){
			my $filename = $_;
			my $read_num = 0;
			foreach(@{$name{$name}}){
				my $temp_name = $_;
				$temp_name =~ s/\;size=\d+\;//;
				if($is{$temp_name}{$filename}){
					$read_num += $is{$temp_name}{$filename};
				}
			}
			print OUT3 "\t$read_num";
		}
		print OUT3 "\n";
	}
}
print "The number of OTUs = $num\n\n";
close(OUT1);
close(OUT2);
close(OUT3);
