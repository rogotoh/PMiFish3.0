#!/usr/bin/perl
use strict;
use warnings;

#2025/01/29 

#Setting
my ($db, $rarefy, $depth, $denoise, $algo, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$algo = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^Type\s*=\s*(\S+)/){$type = $1;}
}
close(SET);
unless($depth){print "Error: Please check the Depth value in Setting.txt.\n"; exit;}
unless($denoise){$denoise = "no";}
if($denoise =~ /^no$/i){exit;}
unless($algo){$algo = 1;}
unless($type){$type = 1;}

#Get data names
if(-e ".\/Results\/2_3_Separate_chimera"){
	opendir (DIR, ".\/Results\/2_3_Separate_chimera") or die ("error:$!");
}else{
	opendir (DIR, ".\/Results\/2_1_Find_unique_in_each_sample") or die ("error:$!");
}
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_zotu_nonchimeras.fa/){$file{$1}++;}
	if ($_ =~ /(.+)_uniques.fa/){$file{$1}++;}
}
closedir DIR;

#usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}

#vsearch editing
if($vv){
	foreach(sort keys %file){
		my $file = $_;
		my $data = "\n";
		if(-f ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa"){
			open (DATA, "<", ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa") or die("error:$!");
			while(<DATA>){
				chomp($_);
				$_ =~ s/\r//g;
				if($_ =~ /^>/){
					$data .= "\n$_\;\n";
				}else{
					$data .= "$_";
				}
			}
			close(DATA);
			$data =~ s/^\n\n//;
			open (OUT, ">", ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa") or die("error:$!");
			print OUT "$data\n";
			close(OUT);
		}
		$data = "\n";
		if(-f ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa"){
			open (DATA, "<", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
			while(<DATA>){
				chomp($_);
				$_ =~ s/\r//g;
				if($_ =~ /^>/){
					$data .= "\n$_\;\n";
				}else{
					$data .= "$_";
				}
			}
			close(DATA);
			$data =~ s/^\n\n//;
			open (OUT, ">", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
			print OUT "$data\n";
			close(OUT);
		}
	}
}

#depth filter when Denoise option was not selected.
unless(-e ".\/Results\/2_3_Separate_chimera"){
	foreach(sort keys %file){
		my ($data, $dseq, @lead);
		my $file = $_;
		open (DATA, "<", ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa") or die("error:$!");
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			unless($_ =~ /^\n/){
				if($_ =~ /^>(.+\;)/){
					if($data){push(@lead, "$data\n$dseq"); $data = $1; undef($dseq);}
					else{$data = $1;}
				}else{
					if($dseq){$dseq = $dseq . $_;}
					else{$dseq = $_;}
				}
			}
		}
		close(DATA);
		if($data){push(@lead, "$data\n$dseq");}
		
		my $sizecount = 0;
		foreach(@lead){
			if($_ =~ /size=(\d+)/){$sizecount += $1;}
		}
		my $tempdepth = $depth;
		if($depth =~ /(\S+)%/){$tempdepth = int($1/100 * $sizecount);}
		
		open (OUT, ">", ".\/Results\/2_1_Find_unique_in_each_sample/${file}_uniques.fa") or die("error:$!");
		foreach(@lead){
			if($_ =~ /size=(\d+)/ and $1 >= $tempdepth){print OUT ">$_\n";}
		}
		close(OUT);
	}
}

my $as_us;
$as_us = "Unique";
print "============================================================\n";
print "                       2_4_0_Uniques                         \n";
print "============================================================\n";

#Get data
mkdir ".\/Results\/2_4_0_${as_us}s";
my (%seq, %count, %table);
foreach(sort keys %file){
	my $file = $_;
	if(-f ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa"){
		open (DATA, "<", ".\/Results\/2_3_Separate_chimera\/${file}_zotu_nonchimeras.fa") or die("error:$!");
	}else{
		open (DATA, "<", ".\/Results\/2_1_Find_unique_in_each_sample/${file}_uniques.fa") or die("error:$!");
	}
	my ($data, $dseq);
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		unless($_ =~ /^\n/){
			if($_ =~ /^>(.+\;)/){
				my $title = $1;
				if($data){
					push (@{$seq{$dseq}}, "${file}\;$data"); 
					$data =~ /size=(\d+)/;
					$table{$dseq}{$file} = $1;
					unless($count{$dseq}){$count{$dseq} = $1;}
					else{$count{$dseq} += $1;}
					$data = $title; 
					undef($dseq);
				}else{$data = $title;}
			}else{
				if($dseq){$dseq = $dseq . $_;}
				else{$dseq = $_;}
			}
		}
	}
	close(DATA);
	if($data){
		push (@{$seq{$dseq}}, "${file}\;$data");
		$data =~ /size=(\d+)/;
		$table{$dseq}{$file} = $1;
		unless($count{$dseq}){$count{$dseq} = $1;}
		else{$count{$dseq} += $1;}
	}
	undef($data);
	undef($dseq);
}

open (OUT1, ">", ".\/Results\/2_4_0_${as_us}s/${as_us}s_seq.fa") or die("error:$!");
open (OUT2, ">", ".\/Results\/2_4_0_${as_us}s/${as_us}s_seq_list.txt") or die("error:$!");
open (OUT3, ">", ".\/Results\/2_4_0_${as_us}s/${as_us}s_table.tsv") or die("error:$!");

my %narabi;
foreach(keys %seq){
	my $temp = $_;
	push(@{$narabi{$count{$temp}}}, $temp);
}
my $otu_num = 0;
#row titles
foreach(sort keys %file){
	print OUT3 "\t$_";
}
print OUT3 "\n";

foreach(sort {$b <=> $a} keys %narabi){
	my $kazu = $_;
	foreach(@{$narabi{$kazu}}){
		my $seq = $_;
		$otu_num++;
		print OUT1 ">$as_us$otu_num\;size=$kazu\;\n$_\n";
		print OUT2 ">$as_us$otu_num\;size=$kazu\;\n";
		foreach(@{$seq{$seq}}){
			print OUT2 "\t$_\n";
		}
		print OUT3 "$as_us$otu_num";
		foreach(sort keys %file){
			unless($table{$seq}{$_}){print OUT3 "\t0";}
			else{print OUT3 "\t$table{$seq}{$_}";}
		}
		print OUT3 "\n";
	}
}
print "\n${as_us}s across all samples= $otu_num\n\n";
close(OUT1);
close(OUT2);
close(OUT3);
