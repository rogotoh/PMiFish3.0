#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($depth, $size, $dictionary, @setting, $denoise);
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
opendir (DIR, ".\/Results\/2_1_Find_unique_in_each_sample") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_uniques.fa/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "                        2_2_Denoise                         \n";
print "============================================================\n";

#Denoise
mkdir ".\/Results\/2_2_Denoise";
foreach(sort keys %file){
	my $file = $_;
	print "\n\n$_\n";
	
	open (DATA, "<", ".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa") or die("error:$!");
	my $count = 0;
	my $sizecount = 0;
	while(<DATA>){if($_ =~ /^>.+size=(\d+)/){$count++; $sizecount += $1;}}
	close(DATA);
	
	if($count < 1){
		open (DATA, ">", ".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt") or die("error:$!");
		close(DATA);
		open (DATA, ">", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		close(DATA);
		print "$file is no sequences.\n";
		next;
	}
	my $tempdepth = $depth;
	if($depth =~ /(\S+)%/){$tempdepth = int($1/100 * $sizecount);}
	
	my $command;
	if($vv){
		$command = ".\/Tools\/$usearch --cluster_unoise \".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa\" --minsize $tempdepth --centroids \".\/Results\/2_2_Denoise\/${file}_zotu.fa\" --unoise_alpha $alpha";
	}else{
		$command = ".\/Tools\/$usearch -unoise3 \".\/Results\/2_1_Find_unique_in_each_sample\/${file}_uniques.fa\" -minsize $tempdepth -ampout \".\/Results\/2_2_Denoise\/${file}_zotu.fa\" -tabbedout \".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt\" -unoise_alpha $alpha";
	}
	system $command;
	
	#correct_error
	unless($vv){
		open (DATA, "<", ".\/Results\/2_2_Denoise\/${file}_unoise3_result.txt") or die("error:$!");
		my %count;
		while(<DATA>){
			if($_ =~ /amp/){
				$_ =~ /(Uniq\d+)/; my $name = $1;
				$_ =~ /size=(\d+)/; my $size = $1;
				$count{$name} = $size;
			}elsif($_ =~ /dqt/){
				$_ =~ /top=(Uniq\d+)/; my $name = $1;
				$_ =~ /size=(\d+)/; my $size = $1;
				$count{$name} += $size;
			}elsif($_ =~ /chfilter/){last;}
		}
		close(DATA);
		
		open (DATA, "<", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		my @tempdata;
		while(<DATA>){
			if($_ =~ /^>/){
				$_ =~ /(Uniq\d+)/; my $name = $1;
				if($count{$name}){$_ =~ s/size=\d+/size=$count{$name}/;}
				push(@tempdata, $_);
			}else{push(@tempdata, $_);}
		}
		close(DATA);
		
		open (OUT, ">", ".\/Results\/2_2_Denoise\/${file}_zotu.fa") or die("error:$!");
		foreach(@tempdata){print OUT "$_";}
		close(OUT);
	}
}
