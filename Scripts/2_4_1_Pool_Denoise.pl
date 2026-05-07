#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($depth, $size, $dictionary, @setting, $correct, $denoise);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	push(@setting, $_);
}
close(SET);
unless($depth){print "Error: Please check the Depth value in Setting.txt.\n"; exit;}
unless($size){print "Error: Please check the Length value in Setting.txt.\n"; exit;}
$correct = "yes";
unless($denoise){$denoise = "no";}
if($denoise =~ /^no$/i){exit;}

#Usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1; last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH executable file is placed in the Tools directory.\n"; exit;}

#Usearch options
my ($alpha);
open (SET, "<", "Options_usearch.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^\t*-unoise_alpha\s*(\S+)/){
		unless($vv){$alpha = $1; last;}
		$alpha = $1;
	}
}
close(SET);
unless($alpha){$alpha = 2.0;}


#Get data names

print "============================================================\n";
print "                    2_4_1_Pool_Denoise                      \n";
print "============================================================\n";

#Denoise
mkdir ".\/Results\/2_4_1_Pool_Denoise";
open (DATA, "<", ".\/Results\/2_4_0_Uniques\/Uniques_seq.fa") or die("error:$!");
my $count = 0;
my $sizecount = 0;
while(<DATA>){if($_ =~ /^>.+size=(\d+)/){$count++; $sizecount += $1;}}
close(DATA);

if($count < 1){
	open (DATA, ">", ".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.txt") or die("error:$!");
	close(DATA);
	open (DATA, ">", ".\/Results\/2_4_1_Pool_Denoise\/zotu.fa") or die("error:$!");
	close(DATA);
	print "No sequences.\n";
}
my $tempdepth = $depth;
if($depth =~ /(\S+)%/){$tempdepth = int($1/100 * $sizecount);}
	
my $command;
if($vv){
	$command = ".\/Tools\/$usearch --cluster_unoise \".\/Results\/2_4_0_Uniques\/Uniques_seq.fa\" --minsize $tempdepth --centroids \".\/Results\/2_4_1_Pool_Denoise\/zotu.fa\" --unoise_alpha $alpha";
}else{
	$command = ".\/Tools\/$usearch -unoise3 \".\/Results\/2_4_0_Uniques\/Uniques_seq.fa\" -minsize $tempdepth -ampout \".\/Results\/2_4_1_Pool_Denoise\/zotu.fa\" -tabbedout \".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.txt\" -unoise_alpha $alpha";
}
system $command;
	
open (CSV, "<", ".\/Results\/2_4_0_Uniques\/Uniques_table.tsv") or die("error:$!");
my (%table, @sample_n);
$count = 0;
while(<CSV>){
	chomp($_);
	$_ =~ s/\r//g;
	my @temp = split(/\t/, $_);
	if($count == 0){
		@sample_n = @temp;
		shift @sample_n;
		$count++;
	}else{
		for(my $n = 1; $n < @temp; $n++){
			$table{$sample_n[$n-1]}{$temp[0]} = $temp[$n];
		}
	}
}
close(CSV);

my @uniq_n;
my %count;

if($vv){
	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/zotu.fa") or die("error:$!");
	while(<DATA>){
		if($_ =~ />(Unique\d+)/){
			push(@uniq_n, $1);
		}
	}
	close(DATA);
}else{
	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.txt") or die("error:$!");
	while(<DATA>){
		if($_ =~ /amp/){
			$_ =~ /(Unique\d+)/; my $name = $1;
			$_ =~ /size=(\d+)/; my $size = $1;
			$count{$name} = $size;
			push(@uniq_n, $name);
		}elsif($_ =~ /dqt/){
			$_ =~ /top=(Unique\d+)/; my $name = $1;
			$_ =~ /^(Unique\d+)/; my $noize = $1;
			$_ =~ /size=(\d+)/; my $size = $1;
			$count{$name} += $size;
			foreach(@sample_n){
				if($table{$_}){
					$table{$_}{$name} += $table{$_}{$noize};
				}
			}
		}elsif($_ =~ /chfilter/){last;}
	}
	close(DATA);
}

open (TSV, ">", ".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.tsv") or die("error:$!");
foreach(@sample_n){
	print TSV "\t$_";
}
print TSV "\n";
foreach(@uniq_n){
	my $temp = $_;
	print TSV "$temp";
	foreach(@sample_n){
		print TSV "\t$table{$_}{$temp}"
	}
	print TSV "\n";
}
close(TSV);

unless($vv){
	my @tempdata;
	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/zotu.fa") or die("error:$!");
	while(<DATA>){
		if($_ =~ /^>/){
			$_ =~ /(Unique\d+)/; my $name = $1;
			if($count{$name}){$_ =~ s/size=\d+/size=$count{$name}/;}
			push(@tempdata, $_);
		}else{push(@tempdata, $_);}
	}
	close(DATA);

	open (OUT, ">", ".\/Results\/2_4_1_Pool_Denoise\/zotu.fa") or die("error:$!");
	foreach(@tempdata){print OUT "$_";}
	close(OUT);
}