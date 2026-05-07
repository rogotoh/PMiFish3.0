#!/usr/bin/perl
use strict;
use warnings;

#option setting
my ($file, $forward, $reverse);
for (my $v = 0; $v < @ARGV; $v++){
	if($ARGV[$v] =~ /-in$/){$file = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-f$/){$forward = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-r$/){$reverse = $ARGV[$v+1];}
}
unless($file){print "Error: Please check file name\n"; exit;}
unless($forward){print "Error: Please check forward primer\n"; exit;}
unless($reverse){print "Error: Please check reverse primer\n"; exit;}

my $rev_f = $forward;
$rev_f = reverse($rev_f);
$rev_f =~ tr/[a-z]/[A-Z]/;
$rev_f =~ tr/ATGCURYMKDHBV/TACGAYRKMHDVB/;

my $rev_r = $reverse;
$rev_r = reverse($rev_r);
$rev_r =~ tr/[a-z]/[A-Z]/;
$rev_r =~ tr/ATGCURYMKDHBV/TACGAYRKMHDVB/;

$forward = &degenerate($forward);
$reverse = &degenerate($reverse);
$rev_f = &degenerate($rev_f);
$rev_r = &degenerate($rev_r);

open (DATAFILE, "<", $file) or die ("error:$!");
open (OUT, ">", "amplified_region.fas") or die ("error:$!");

my $name;
while(<DATAFILE>){
	$_ =~ s/\r//;
	chomp($_);
	if ($_ =~ /^>/){$name = $_;}
	else{
		if($_=~ /$forward(.+)$rev_r/){
			print OUT "$name\n$1\n";
		}elsif($_ =~ /$reverse(.+)$rev_f/){
			my $temp = $1;
			my $reverse = reverse($temp);
			$reverse =~ tr/ATGC/TACG/;
			print OUT "$name\n$reverse\n";
		}
	}
}
close(DATAFILE);
close(OUT);
print "amplified_region.fas was created\n\n";

#Degenerate
sub degenerate {
	my ($dege) = @_;
	$dege =~ s/B/[CGT]/g; $dege =~ s/D/[AGT]/g; $dege =~ s/H/[ACT]/g; $dege =~ s/K/[GT]/g; 
	$dege =~ s/M/[AC]/g; $dege =~ s/N/[ACGT]/g; $dege =~ s/R/[AG]/g; $dege =~ s/S/[CG]/g;
	$dege =~ s/V/[ACG]/g; $dege =~ s/W/[AT]/g; $dege =~ s/Y/[CT]/g;
	$dege =~ s/\./\.\*/g;
	return($dege);
}
