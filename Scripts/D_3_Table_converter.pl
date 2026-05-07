#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($depth, $size, $temporary, $compress);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /^Compress\s*=\s*(\S+)/){$compress = $1;}
}
close(SET);

unless(-d ".\/Results"){mkdir ".\/Results";}
mkdir ".\/Results\/2_4_ASVs";

#data
open (DATA,"<", ".\/Results\/ASV_table_original.csv") or die ("error:$!");
open (OUT1, ">", ".\/Results\/2_4_ASVs/ASVs_seq.fa") or die("error:$!");
open (OUT2, ">", ".\/Results\/2_4_ASVs/ASVs_table.tsv") or die("error:$!");
my $count = 0;
my (%data, @seq, %sample, @sample, %length, %depth);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//;
	$_ =~ s/\"//g;
	if($count == 0){
		$count = 1;
		my @temp = split(/,/, $_);
		for(my $n = 1; $n < @temp; $n++){
			if(length($temp[$n]) < $size){$length{$temp[$n]}++;} #length filter
			$data{$temp[$n]} = 0;
			push(@seq, $temp[$n]);
		}
	}else{
		my @temp = split(/,/, $_);
		push(@sample, $temp[0]);
		for(my $n = 1; $n < @temp; $n++){
			unless($temp[$n] < $depth){ #depth filter
				$data{$seq[$n-1]} += $temp[$n];
				$sample{$seq[$n-1]}{$temp[0]} = $temp[$n];
			}else{
				$sample{$seq[$n-1]}{$temp[0]} = 0;
			}
		}
	}
}
close(DATA);

foreach(@sample){
	print OUT2 "\t$_";
}
print OUT2 "\n";


my %narabi;
foreach(keys %data){
	push(@{$narabi{$data{$_}}}, $_);
}

my $unique = 1;
foreach(sort {$b <=> $a} keys %narabi){
	my $temp = $_;
	foreach(@{$narabi{$temp}}){
		unless($length{$_}){
			my $seq = $_;
			unless($temp < $depth){
				print OUT1 ">ASV$unique\;size=$temp\;\n$_\n";
				print OUT2 "ASV$unique";
				foreach(@sample){
					print OUT2 "\t$sample{$seq}{$_}";
				}
				print OUT2 "\n";
				$unique++;
			}
		}
	}
}
close(OUT1);
close(OUT2);

print "============================================================\n";
print "                         2_4_ASVs                           \n";
print "============================================================\n";
print "\nASVs across all samples = $unique\n\n";

if($temporary =~ /yes/i){
	my @sakujo = glob(".\/Results\/1_1_Primer_Trimmed_fastq\/*.fq");
	unlink(@sakujo);
	@sakujo = glob(".\/Results\/1_2_Filtered\/*.fastq.gz");
	unlink(@sakujo);
}

if($compress =~ /yes/i){
	my @sakujo = glob(".\/Run\/*.fastq");
	print "\nCompressing fastq file to gz file...\n";
	foreach(@sakujo){
		my $comand = "gzip \"$_\"";
		system $comand;
	}
}
