#!/usr/bin/perl
use strict;
use warnings;

my @files;
opendir (DIR, "./") or die ("error:$!");
while (readdir DIR) {
    if ($_ =~ /([\w\W]+).gb$/){push (@files, $1);}
}
closedir DIR;

unless(@files){print "Error: Please put gb file in this directory\n"; exit;}

my $number = 0;
my (%narabi,%species, %taxid);
foreach (@files){
	my $file = $_;
	open (DATAFILE, "<", "${file}.gb") or die("error:$!");
	my $count = 0;
	my $definit = 0;
	my $reference = 0;
	my ($country, $source);
	my ($seq, $name, $id);
	while (<DATAFILE>){
		if ($_ =~ /^\/\//){
			if($count == 1){
				#unless($id =~ /NC_/){
					$seq =~ tr/atgc/ATGC/;
					$number++;
					$species{$name}++;
					$narabi{$name}{$id} = ">gb|${id}|${name}\n$seq\n";
				#}
			}
			undef($seq);
			undef($name);
			undef($country);
			$count = 0;
		}
		if ($_ =~ /ACCESSION\s+(\D+\d+)/){$id = $1;}
		if ($name and $reference){
			unless($_ =~ /\;/){
				$_ =~ /\s+(.+)\n/;
				$name = $name . " $1";
			}
			$reference = 0;
		}
		if ($_ =~ /ORGANISM\s+(.+)\n/){$name = $1;$name =~ s/\s/_/g; $reference = 1;}
		if ($_ =~ /DEFINITION\s+(.+)\n/){$definit = $1;}
		if ($_ =~ /SOURCE\s+(.+)\n/){$source = $1;}
		if ($_ =~ /\/db_xref=\"taxon:(.+)\"/){$taxid{$1}++;}
		if ($count == 1){
			$_ =~ s/[\s\d]//g;
			unless($seq){$seq = $_;}
			else{$seq = $seq.$_;}
		}
		if ($_ =~ /ORIGIN/){$count++;}
	}
}

open (OUT, ">", "database.fas") or die("error:$!");
open (OUT2, ">", "species_list.txt") or die("error:$!");
foreach (sort {$a cmp $b} keys %narabi){
	my $sp = $_;
	my $temp = $narabi{$sp};
	my $datasu = 0;
	foreach (sort {$a cmp $b} keys %$temp){
		print OUT "$narabi{$sp}{$_}";
		$datasu++;
	}
	print OUT2 "$sp\t$datasu\n";
}
close(OUT);
close(OUT2);

#taxid list
open (OUT3, ">", "taxid_list.txt") or die("error:$!");
foreach (sort {$a <=> $b} keys %taxid){
	print OUT3 "$_\n";
}
close(OUT3);

my $sp_num = keys %species;
print "$number sequences, $sp_num species\n";
print "database.fas was created\n\n";
