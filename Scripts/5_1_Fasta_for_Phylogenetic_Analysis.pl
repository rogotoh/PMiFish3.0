#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($primer, $separate, $denoise, $type);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
}
close(SET);
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($primer =~ /^no$/i){
	unless(-f ".\/DataBase\/$primer"){print "Error: Primer file not found or primer file name mismatch in Setting.txt.\n"; exit;}
}
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

#Get primer sequences
my $countp = 0;
my (%forward, %reverse, $namep, @forward);
my ($primerF, $primerR);
my $primer_num;

unless($primer =~ /^no$/i){
	open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
	while(<DATA2>){
		if($_ =~ /\#/){next;}
		$_ =~ s/\r//g;
		chomp($_);
		if($_ =~ /Forward/){$countp = 1;next;}
		if($_ =~ /Reverse/){$countp = 2;next;}
		if($countp == 1){
			if($_ =~ /^>(.+)/){$forward{$1}++; push(@forward, $1);}
		}
		if($countp == 2){
			if($_ =~ /^>(.+)/){$reverse{$1}++;}
		}
	}
	close(DATA2);
	unless(%forward){
		print "Error: No primer sequence found in $primer. Please check the primer file.\n";
		exit;
	}
	unless(%reverse){
		print "Error: No primer sequence found in $primer. Please check the primer file.\n";
		exit;
	} 
	foreach(keys %forward){
		unless($reverse{$_}){
			print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
			exit;
		}
	}
	foreach(keys %reverse){
		unless($forward{$_}){
			print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
			exit;
		}
	}
	$primer_num = keys %forward;
	if($primer_num == 1){$separate = 0;}
}


#Get data names
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
my @file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;

print "============================================================\n";
print "            5_1_Fasta_for_Phylogenetic_Analysis             \n";
print "============================================================\n";

#Annotation
mkdir ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis";
if($separate){
	my (%merged, %all);
	foreach(sort keys %forward){
		my $temp_primer = $_;
		foreach(sort keys %file){
			if($_ =~ /$temp_primer$/){
				my $file = $_;
				open (DATA, "<", "./Results/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
				my $fname;
				while(<DATA>){
					chomp($_);
					if($_ =~ /^>(.+)/){$fname = $1 . "_$file";}
					else{$merged{$_}{$fname}++; push(@{$all{$temp_primer}}, ">$fname\n$_\n");}
				}
				close(DATA);
			}
		}
	}

	foreach(sort keys %all){
		my $temp_primers = $_;
		open (OUT, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/${temp_primers}_representative_seqs.fas") or die("error:$!");
		my $temp = $all{$temp_primers};
		my @temp = @$temp;
		foreach(sort @temp){
			unless($_ =~ /^>Nohit/){
				print OUT "$_";
			}
		}
		close(OUT);
	}

	foreach(sort keys %all){
		my $temp_primer = $_;
		open (OUT2, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/${temp_primer}_merged_seq.fas") or die("error:$!");
		open (OUT3, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/${temp_primer}_merged_list.txt") or die("error:$!");
		my (@narabi1, @narabi2, %haplo);
		foreach(keys %merged){
			my $seq = $_;
			my $hash = $merged{$_};
			my %names = %$hash;
			my @keys = keys %names;
			my @narabikae;
			my @temp;
			foreach(@keys){
				if($_ =~ /$temp_primer/){
					push(@temp, $_);
				}
			}
			my $t = @temp;
			@keys = @temp;
			unless(@keys){next;}
			if($keys[0] =~ /Nohit/){next;}
			foreach(sort keys %file){
				my $itiji = $_;
				foreach(@keys){
					if($_ =~ /$itiji/){push (@narabikae, $_); last;}
				}
			}
			my $num = @keys;
			if($num > 1){
				my ($temp, $temp2);
				if($keys[0] =~ /(.+)_OTU\d+_\d+_reads/){
					$temp = $1;
				}elsif($keys[0] =~ /(.+)_${as_us}\d+_\d+_reads/){
					$temp = $1;
				}
				$temp =~ s/_otu\d+//;
				$haplo{$temp}++;
				if($haplo{$temp} > 1){
					push(@narabi1, ">${temp}_h$haplo{$temp}_from_${num}_sites\n$seq\n");
					$temp2 = ">${temp}_h$haplo{$temp}_from_${num}_sites\n";
				}else{
					push(@narabi1, ">${temp}_from_${num}_sites\n$seq\n");
					$temp2 = ">${temp}_from_${num}_sites\n";
				}
				foreach(@narabikae){$temp2 = $temp2 . "\t$_\n";}
				push(@narabi2, $temp2);
			}else{
				push(@narabi1, ">$keys[0]\n$seq\n");
			}
		}
		foreach(sort @narabi1){print OUT2 "$_";}
		foreach(sort @narabi2){print OUT3 "$_";}
		my $count = @narabi1;
		print "${temp_primer}: $count Mereged Sequences\n";
		close(OUT2);
		close(OUT3);
	}
}else{
	my (%merged, %all);
	foreach(sort keys %file){
		my $file = $_;
		open (DATA, "<", "./Results/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
		my $fname;
		while(<DATA>){
			chomp($_);
			if($_ =~ /^>(.+)/){$fname = $1 . "_$file";}
			else{$merged{$_}{$fname}++; $all{$fname} = $_;}
		}
		close(DATA);
	}

	open (OUT, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/all_representative_seqs.fas") or die("error:$!");
	foreach(sort keys %all){
		if($_ =~ /Nohit/){next;}
		print OUT ">$_\n$all{$_}\n";
	}
	close(OUT);

	open (OUT2, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq.fas") or die("error:$!");
	open (OUT3, ">", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_list.txt") or die("error:$!");

	my (@narabi1, @narabi2, %haplo);
	foreach(keys %merged){
		my $seq = $_;
		my $hash = $merged{$_};
		my %names = %$hash;
		my @keys = keys %names;
		my @narabikae;
		foreach(sort keys %file){
			my $itiji = $_;
			foreach(@keys){
				if($_ =~ /$itiji/){push (@narabikae, $_); last;}
			}
		}
		my $num = @keys;
		if($keys[0] =~ /Nohit/){next;}
		if($num > 1){
			my ($temp, $temp2);
			if($keys[0] =~ /(.+)_OTU\d+_\d+_reads/){
				$temp = $1;
			}elsif($keys[0] =~ /(.+)_${as_us}\d+_\d+_reads/){
				$temp = $1;
			}
			$temp =~ s/_otu\d+//;
			$haplo{$temp}++;
			if($haplo{$temp} > 1){
				push(@narabi1, ">${temp}_h$haplo{$temp}_from_${num}_sites\n$seq\n");
				$temp2 = ">${temp}_h$haplo{$temp}_from_${num}_sites\n";
			}else{
				push(@narabi1, ">${temp}_from_${num}_sites\n$seq\n");
				$temp2 = ">${temp}_from_${num}_sites\n";
			}
			foreach(@narabikae){$temp2 = $temp2 . "\t$_\n";}
			push(@narabi2, $temp2);
		}else{
			push(@narabi1, ">$keys[0]\n$seq\n");
		}
	}
	foreach(sort @narabi1){print OUT2 "$_";}
	foreach(sort @narabi2){print OUT3 "$_";}
	my $count = @narabi1;
	print "$count Mereged Sequences\n";

	close(OUT2);
	close(OUT3);
}
