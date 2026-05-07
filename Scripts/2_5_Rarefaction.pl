#!/usr/bin/perl
use strict;
use warnings;
use List::Util;

#Setting
my ($rarefy, $denoise, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Rarefaction\s*=\s*(\S+)/){$rarefy = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
}
close(SET);
if($rarefy =~ /^no$/i){exit;}
unless($denoise){$denoise = "no";}
unless($type){$type = 1;}


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


print "============================================================\n";
print "                      2_5_Rarefaction                       \n";
print "============================================================\n";

#Rarefy
mkdir ".\/Results\/2_5_Rarefaction";

#ASVs_table
my (@data, $sample, @asv);
open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv") or die("error:$!");
my $count = 0;
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	$count++;
	if($count == 1){
		$sample = $_;
	}else{
		my @temp = split(/\t/, $_);
		$temp[0] =~ /${as_us}(\d+)/;
		my $num = $1;
		push(@asv, $num);
		for(my $n = 1; $n < @temp; $n++){
			if($count == 2){$data[$n-1] = "$num\t" x $temp[$n];}
			else{$data[$n-1] .= "$num\t" x $temp[$n];}
		}
	}
}
close(DATA);

#rarefaction
open (OUT, ">", ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv") or die("error:$!");
print OUT "$sample\n";

$count = 0;
my @filename = split(/\t/, $sample);
my (%read, %total);
foreach(@data){
	my $lead = $_;
	$lead =~ s/\t$//;
	my @lead = split(/\t/, $lead);
	my $length = @lead;
	srand(123);
	@lead = List::Util::shuffle @lead;
	@lead = splice (@lead, 0, $rarefy);
	if($length < $rarefy){
		print "$filename[$count+1]\t$length Reads\n\tWarning: Reads are less than $rarefy (Rarefaction)\n";
	}else{
		print "$filename[$count+1]\t$rarefy Reads\n";
	}
	foreach(@lead){
		$read{$count}{$_}++;
		$total{$_}++;
	}
	$count++;
}

#table_output
foreach(@asv){
	my $asv = $_;
	if($total{$asv}){
		print OUT "${as_us}$asv";
		for(my $n = 0; $n < @data; $n++){
			if($read{$n}{$asv}){
				print OUT "\t$read{$n}{$asv}";
			}else{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
}
close(OUT);

#fasta output
open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa") or die("error:$!");
my $check = 0;
my (%narabi, $name, $size);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>${as_us}(\d+)/){
		$check = 0;
		if($total{$1}){
			$name = ">${as_us}$1\;size=$total{$1}\;\n";
			$check = 1;
			$size = $total{$1};
		}
	}else{
		if($check == 1){
			$name .= "$_\n";
			push(@{$narabi{$size}}, $name);
		}
	}
}
close(DATA);

open (OUT, ">", ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa") or die("error:$!");
foreach(sort {$b <=> $a} keys %narabi){
	my $temp = $_;
	foreach(@{$narabi{$temp}}){
		print OUT "$_";
	}
}
close(OUT);

