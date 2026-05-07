#!/usr/bin/perl
use strict;
use warnings;

#perl cluster_maker.pl -in database.fas -id cut_line(ex.97)
#delete the same sequence of the same species
#option setting
my ($file, $identity, $dis);
for (my $v = 0; $v < @ARGV; $v++){
	if($ARGV[$v] =~ /-in$/){$file = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-id$/){$identity = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-d$/){$dis = $ARGV[$v+1];}
}
unless($file){print "Error: Please check file name\n"; exit;}
unless($dis){
	unless($identity){print "Error: Please input cluster indentity value (-id <number>)\n"; exit;}
}
unless($identity){
	unless($dis){print "Error: Please input difference value (-d <number>)\n"; exit;}
}

#userch check
opendir (DIR, "..\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv, $sw);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
	elsif($_ =~ /(swarm.*)/){$sw = $1;}
}
closedir DIR;
if($identity){
	unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}
}
if($dis){
	unless($sw){print "Error: Please ensure the SWARM executable file is placed in the Tools directory.\n"; exit;}
}
#file reading
$file =~ /([^\.|^\\]+)\.fa/;
my $filename = $1;
unless($filename){$filename = $file;}

open(DATA, "<", $file) or die("error:$!");
my ($name, %data, $original, %modosi, %seq);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){
		$original = $_;
		$_ =~ s/.+\|//;
		$name = $_;
	}else{
		#delete the same sequence of the same species
		unless($seq{$name}{$_}){
			$seq{$name}{$_}++;
			$data{$name}{$original} = $_;
		}
	}
}
close(DATA);

open (OUT, ">", ".\/${filename}_clusterized.fas") or die("error:$!");
print "Now processing...\n";

foreach(sort {$a cmp $b} keys %data){
	my $spname = $_;
	my $oriname = $data{$_};
	my $num = keys %$oriname;
	if($num > 1){
		open (TEMP, ">", ".\/temp.fas") or die("error:$!");
		my $count = 0;
		foreach(sort keys %$oriname){
			$count++; 
			print TEMP ">Uniq$count;size=$num;\n$data{$spname}{$_}\n";
			$modosi{"Uniq$count;size=$num;"} = $_;
			$num -= 1;
		}
		close(TEMP);
		
		my $command;
		if($sw){
			$command = "..\/Tools\/$sw \".\/temp.fas\" -d $dis -u \".\/temp.uc\" -z";
		}elsif($vv){
			my $tempid = $identity/100;
			$command = "..\/Tools\/$usearch --cluster_smallmem \".\/temp.fas\" --id $tempid --uc \".\/temp.uc\" --usersort";
		}else{
			my $tempid = $identity/100;
			$command = "..\/Tools\/$usearch -cluster_smallmem \".\/temp.fas\" -id $tempid -uc \".\/temp.uc\" -sortedby size -quiet";
		}
		system $command;
		
		open (TEMP, "<", ".\/temp.uc") or die("error:$!");
		my $sp = 0;
		my (%otus, %out);
		while(<TEMP>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			if($temp[0] eq "S"){
				$sp++;
				$out{"$modosi{$temp[8]}_cluster$sp"} = $data{$spname}{$modosi{$temp[8]}};
				$otus{$temp[8]} = "cluster$sp";
			}elsif($temp[0] eq "H"){
				$out{"$modosi{$temp[8]}_$otus{$temp[9]}"} = $data{$spname}{$modosi{$temp[8]}};
			}
		}
		close(TEMP);
		
		foreach(sort keys %out){
			my $out = $_;
			if($sp == 1){
				$out =~ s/_cluster\d+//;
				print OUT "$out\n$out{$_}\n";
			}else{
				print OUT "$_\n$out{$_}\n";
			}
		}
		unlink ".\/temp.fas";
		unlink ".\/temp.uc";
	}else{
		my @keys = keys %$oriname;
		print OUT "$keys[0]\n$data{$spname}{$keys[0]}\n";
	}
}
close(OUT);
print "${filename}_clusterized.fas was created\n\n";
