#!/usr/bin/perl
use strict;
use warnings;

#You need names.dmp/nodes.dmp download from below URL(taxdump.tar.gz)
#https://ftp.ncbi.nih.gov/pub/taxonomy/

#dmp file check
opendir (DIR, ".\/") or die ("error:$!");
my @tool = readdir DIR;
my ($node_dmp,$name_dmp);
foreach (@tool) {
	if($_ =~ /names.dmp/){$name_dmp = 1;}
	if($_ =~ /nodes.dmp/){$node_dmp = 1;}
}
closedir DIR;
unless($node_dmp){print "Error: Please put the \"nodes.dmp\" file in this directory\n"; exit;}
unless($name_dmp){print "Error: Please put the \"names.dmp\" file in this directory\n"; exit;}

my $target = $ARGV[0];
unless($target){
	print "Please enter taxonomic name\n";
	print "ex) perl Taxonomic_Rank.pl \"Mammalia\"\n";
	exit;
}
$target = ucfirst($target);

#Data Set
print "Now processing...\n";

#names.dmp
my (%name, $targetid, %synonym);
open (DATAFILE, "<", "names.dmp") or die("error:$!");
while (<DATAFILE>){
	if($_ =~ /^(\d+)\t\|\t([^\t]+)\t\|\t.*\t\|\t([^\t]+)\t\|/){
		my $num = $1;
		my $name = $2;
		my $sn = $3;
		if($sn eq "scientific name"){
			$name{$num} = $name;
			if($name eq $target){$targetid = $num;}
		}
		elsif($sn eq "synonym"){
			$synonym{$name} = $num;
		}
	}
}
close(DATAFILE);
if($synonym{$target}){
	print "\nThe taxonomic name \"$target\" you entered is a synonym\n";
	print "The valid name in NCBI is $name{$synonym{$target}}\n";
	$targetid = $synonym{$target};
	$target = $name{$synonym{$target}};
}

unless($targetid){print "Nothing Taxon_ID!\n"; exit;}
else{print "\n$target Taxon_ID = $targetid\; ";}

#nodes.dmp
my (%kankei, %rank, @species);
open (DATAFILE, "<", "nodes.dmp") or die("error:$!");
while (<DATAFILE>){
$_ =~ /^(\d+)\t\|\t(\d+)\t\|\t([^\t]+)\t\|/;
$kankei{$1} = $2;
$rank{$1} = $3;
my $tempid = $1;
my $temprank = $3;
if($temprank eq "species"){push(@species, $tempid);}
}
close(DATAFILE);

if($rank{$targetid}){print "Rank = $rank{$targetid}\n\n";}

#output
my $rank_check = 0;
if($rank{$targetid} eq "species"){$rank_check = 1;}
elsif($rank{$targetid} eq "genus"){$rank_check = 2;}
elsif($rank{$targetid} eq "family"){$rank_check = 3;}
elsif($rank{$targetid} eq "order"){$rank_check = 4;}
elsif($rank{$targetid} eq "class"){$rank_check = 5;}
elsif($rank{$targetid} eq "phylum"){$rank_check = 6;}
elsif($rank{$targetid} eq "kingdom"){$rank_check = 7;}

my @output;
my @rank = qw(genus family order class phylum kingdom);
foreach(sort {$a <=> $b} @species){
	if($rank_check == 1){
		unless($_ == $targetid){next;}
	}
	my $id = $_;
	my $taxon;
	my $norank = 0;
	$taxon = "$name{$id}\t$id";
	if($name{$id} eq $target){$norank = 1;}
	my $check = 1;
	my %narabi;
	while($check == 1){
		$id = $kankei{$id};		#parent id
		my $rank = $rank{$id};	#parent rank
		if($id == 1){last;}
		#if($rank eq "no rank"){
		#	$norank++;
		#	if($norank > 1){last;}
		#	else{next;}
		#}
		if($rank eq $rank{$targetid}){
			$norank = 1;
			unless($name{$id} eq $target){last;}
		}
		if($rank eq "genus"){
			if($rank_check == 2){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
		}elsif($rank eq "family"){
			if($rank_check == 3){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
		}elsif($rank eq "order"){
			if($rank_check == 4){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
		}elsif($rank eq "class"){
			if($rank_check == 5){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
		}elsif($rank eq "phylum"){
			if($rank_check == 6){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
		}elsif($rank eq "kingdom"){
			if($rank_check == 7){
				unless($name{$id} eq $target){last;}
			}
			$narabi{$rank} = $name{$id};
			$check = 2;
		}
	}
	if($check == 2){
		foreach(@rank){
			if($narabi{$_}){
				$taxon = "$narabi{$_}\t" . $taxon;
			}else{
				$taxon = "\t" . $taxon;
			}
		}
		if($norank == 1){push(@output, $taxon);}
		#if($taxon =~ /\s*$target\s*/){push(@output, $taxon);}
		#push(@output, $taxon);
		if($rank_check == 1){last;}
	}
}

open (OUT, ">", "Taxonomic_Ranks.tsv") or die("error:$!");
print OUT "Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tTaxonID\n";
my @rank_n = qw(Kingdom phy cls ord fam gen);
my $number = @output;
foreach(sort @output){
	my $out = $_;
	if($_ =~ /\t\t/){
		my @temp = split(/\t/, $_);
		for(my $n = 0; $n < @temp; $n++){
			unless($temp[$n]){
				$temp[$n-1] =~ /(\S+)/;
				$temp[$n] = "$1 $rank_n[$n] Incertae sedis";
			}
		}
		$out = join("\t", @temp);
	}
	print OUT "$out\n";
}
print "Taxonomic_Ranks.tsv was successfully created.\n";
