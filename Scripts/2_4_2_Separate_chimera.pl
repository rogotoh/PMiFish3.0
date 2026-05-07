#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($temporary, $denoise);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
}
close(SET);
unless($temporary){$temporary = "YES";}
unless($denoise){$denoise = "no";}
if($denoise =~ /^no$/i){exit;}

#Usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH executable file is placed in the Tools directory.\n"; exit;}


#Get data names

print "============================================================\n";
print "                  2_4_2_Separate_chimera                    \n";
print "============================================================\n";
mkdir ".\/Results\/2_4_2_Separate_chimera";

my (@logtxt);
if(-f ".\/Results\/log.txt"){
	open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
	while(<LOG>){
		push(@logtxt, $_);
	}
	close(LOG);
}

#Separate chimeras
mkdir ".\/Results\/2_4_ZOTUs";
if($vv){
	my %uniq_n;
	my %seq_n;
	my $command = ".\/Tools\/$usearch --uchime3_denovo \".\/Results\/2_4_1_Pool_Denoise\/zotu.fa\" --nonchimeras \".\/Results\/2_4_2_Separate_chimera\/zotu_nonchimeras.fa\"";
	system $command;
	
	open (DATA, "<", ".\/Results\/2_4_2_Separate_chimera\/zotu_nonchimeras.fa") or die("error:$!");
	open (OUT, ">", ".\/Results\/2_4_ZOTUs\/ZOTUs_seq.fa") or die("error:$!");
	my $n = 0;
	while(<DATA>){
		if($_ =~ /^>(Unique\d+)/){
			$uniq_n{$1}++;
			if($n > 0){print OUT "\n";}
			$n++;
			$_ =~ s/(Unique\d+)/ZOTU$n/;
			$seq_n{$1} = "ZOTU$n";
			$_ =~ s/amptype.+//;
			print OUT "$_";
		}else{
			chomp($_);
			$_ =~ s/\r//g;
			print OUT "$_";
		}
	}
	print OUT "\n";
	close(DATA);
	close(OUT);
	
	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.tsv") or die("error:$!");
	open (OUT, ">", ".\/Results\/2_4_ZOTUs\/ZOTUs_table.tsv") or die("error:$!");
	my $count = 0;
	while(<DATA>){
		if($_ =~ /^(Unique\d+)/){
			my $temp = $1;
			if($uniq_n{$temp}){
				$_ =~ s/Unique\d+/$seq_n{$temp}/;
				print OUT "$_";
				$count++;
			}
		}else{
			print OUT "$_";
		}
	}
	close(OUT);
	close(DATA);
	print "============================================================\n";
	print "                         2_4_ZOTUs                          \n";
	print "============================================================\n";
	print "\nZOTUs across all samples = ", $count, "\n\n";
	
	open (DATA, "<", ".\/Results\/2_4_ZOTUs\/ZOTUs_table.tsv") or die("error:$!");
	$count = 0;
	my @sample_n;
	my %read;
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		my @temp = split(/\t/, $_);
		if($count == 0){
			@sample_n = @temp;
			$count++;
		}else{
			for(my $n = 1; $n < @temp; $n++){
				$read{$sample_n[$n]} += $temp[$n];
			}
		}
	}
	close(DATA);
	if(-f ".\/Results\/log.txt"){
		$count = 0;
		open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
		foreach(@logtxt){
			if($count > 0){
				chomp($_);
				$_ =~ s/\r//g;
				my @log = split(/\t/, $_);
				my $fname = $log[0];
				my $first = $log[1];
				if($read{$fname}){
					my $per = sprintf("%.2f", $read{$fname}/$first*100);
					print LOG $_ . "\t$read{$fname} ($per%)\n";
				}else{
					print LOG $_ . "\t0 (0.00%)\n";
				}
			}else{
				print LOG "$_";
			}
			$count++;
		}
		close(LOG);
	}
}else{
	my ($data, $dseq, @lead);
	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/zotu.fa") or die("error:$!");
	open (OUT1, ">", ".\/Results\/2_4_2_Separate_chimera\/zotu_chimeras.fa") or die("error:$!");
	open (OUT2, ">", ".\/Results\/2_4_2_Separate_chimera\/zotu_nonchimeras.fa") or die("error:$!");
	open (OUT3, ">", ".\/Results\/2_4_ZOTUs\/ZOTUs_seq.fa") or die("error:$!");
	open (OUT4, ">", ".\/Results\/2_4_ZOTUs\/ZOTUs_list.txt") or die("error:$!");
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
	
	my $count1 = 0;
	my $count2 = 0;
	my $num = 1;
	my %chimera;
	my %zotu;
	my %narabi;
	foreach(@lead){
		if($_ =~ /chimera/){
			print OUT1 ">$_\n";$count1++;
			$_ =~ /(Unique\d+)/;
			$chimera{$1}++;
		}else{
			$_ =~ /size=(\d+)/;
			$narabi{$1}{$_}++;
			print OUT2 ">$_\n";$count2++;
		}
	}
	foreach(sort {$b <=> $a} keys %narabi){
		my $size = $_;
		foreach(sort {$a cmp $b} keys %{$narabi{$size}}){
			my $name = $_;
			$name =~ s/(Unique\d+)/ZOTU$num/;
			my $uniq = $1;
			$zotu{$uniq} = "ZOTU$num";
			$name =~ s/amptype.+//;
			print OUT3 ">$name\n";
			print OUT4 ">ZOTU$num\n $uniq\n";
			$num++;
		}
	}
	
	print "$count2 non-chimeras\n\t$count1 chimeras\n";
	close(OUT1); close(OUT2); close(OUT3); close(OUT4);

	open (DATA, "<", ".\/Results\/2_4_1_Pool_Denoise\/unoise3_result.tsv") or die("error:$!");
	open (OUT, ">", ".\/Results\/2_4_ZOTUs\/ZOTUs_table.tsv") or die("error:$!");
	my %read_c;
	my @sample_n;
	my %narabi2;
	my $count = 0;
	while(<DATA>){
		my @split = split(/\t/, $_);
		if($count == 0){
			@sample_n = @split;
			for(my $n = 1;$n < @sample_n; $n++){
				$read_c{$split[$n]} = 0;
			}
			$count++;
			print OUT "$_";
			next;
		} 
		if($chimera{$split[0]}){next;}
		else{
			if($zotu{$split[0]}){
				$_ =~ s/^$split[0]\t/$zotu{$split[0]}\t/;
				for(my $n = 1;$n < @split; $n++){
					$read_c{$sample_n[$n]} += $split[$n];
				}
			}
			$_ =~ /ZOTU(\d+)/;
			$narabi2{$1} = $_;
		}
	}
	foreach(sort{$a <=> $b} keys %narabi2){print OUT "$narabi2{$_}";}
	close(DATA);
	close(OUT);
	
	if(-f ".\/Results\/log.txt"){
		$count = 0;
		open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
		foreach(@logtxt){
			if($count > 0){
				chomp($_);
				$_ =~ s/\r//g;
				my @log = split(/\t/, $_);
				my $fname = $log[0];
				my $first = $log[1];
				if($read_c{$fname}){
					my $per = sprintf("%.2f", $read_c{$fname}/$first*100);
					print LOG $_ . "\t$read_c{$fname} ($per%)\n";
				}else{
					print LOG $_ . "\t0 (0.00%)\n";
				}
			}else{
				print LOG "$_";
			}
			$count++;
		}
		close(LOG);
	}
	
	print "============================================================\n";
	print "                         2_4_ZOTUs                          \n";
	print "============================================================\n";
	print "\nZOTUs across all samples = ", $num - 1, "\n\n";

}

