#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

#Setting
my ($primer, $diff, $separate, $trim, $shift);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^MaxDiff\s*=\s*(\S+)/){$diff = $1;}
	elsif($_ =~ /^Shift\s*=\s*(\S+)/){$shift = $1;}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
}
close(SET);

if($diff =~ /^length$/i){$trim = 1;}
unless($diff =~ /\d+/){$diff = 2;}
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($primer =~ /^no$/i){
	unless(-f ".\/DataBase\/$primer"){print "Error: Primer file not found or primer file name mismatch in Setting.txt.\n"; exit;}
}

#usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1; last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}

#Decompressing files
opendir (DIR, ".\/Run") or die ("error:$!");
my @run = readdir DIR;
my $gzip = 0;
foreach (@run) {
	if ($_ =~ /\.gz$/){$gzip++;}
}
closedir DIR;

if($gzip){
	print "Decompressing gz file to fastq file...\n";
	foreach (@run) {
		if ($_ =~ /\.gz$/){
			my $comand = "gzip -d \".\/Run\/$_\"";
			system $comand;
		}
	}
}

#rename
opendir (DIR, ".\/Run") or die ("error:$!");
my @read = readdir DIR;
foreach (@read){
	my $temp;
	$temp = $_;
	$_ =~ s/\r//;
	if($_ =~ /_1.fastq|_1.fq/){
		$temp =~ s/_1/_R1/;
	}elsif($_ =~ /_2.fastq|_2.fq/){
		$temp =~ s/_2/_R2/;
	}
	if($temp){rename ".\/Run\/$_", ".\/Run\/$temp";}
}

#Get data names
opendir (DIR, ".\/Run") or die ("error:$!");
@read = readdir DIR;
my $read_count = 0;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_R1/){$file{$1}++; $read_count++;}
	if ($_ =~ /(.+)_R2/){$file{$1}++; $read_count++;}
}
closedir DIR;
unless(%file){print "Error: None of fastq files in Run file.\n"; exit;}
foreach(keys %file){
	if($file{$_} < 2){
		print "Error: Missing paired-end files in the Run directory! Both R1 and R2 fastq files are required.\n";
		print "Missing either R1 or R2 fastq file for sample: $_.\n";
		exit;
	}
}

#primer no
if($primer =~ /^no$/i){
	unless(-d ".\/Results"){mkdir ".\/Results";}
	mkdir ".\/Results\/1_1_Primer_Trimmed_fastq";
	foreach(sort @read){
		my $file = $_;
		if($_ =~ /_R\d_/){
			my ($filename, $comand);
			$_ =~ /(.+)_R\d/;
			$filename = $1;
			my $rname;
			if($file =~ /_R1/){$rname = "R1";}
			else{$rname = "R2";}
			
			copy(".\/Run\/$file", ".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_${rname}_stripped.fq") or die ("Failed to copy file:$!");
		}
	}
	exit;
}

#Get primer sequences
my $countp = 0;
my (%forward, %reverse, $namep, @primername);
my ($primerF, $primerR);

open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
while(<DATA2>){
	if($_ =~ /\#/){next;}
	$_ =~ s/\r//g;
	chomp($_);
	if($_ =~ /Forward/){
		$countp = 1;
		next;
	}
	if($_ =~ /Reverse/){
		$countp = 2;
		next;
	}
	if($countp == 1){
		if($_ =~ /^>(.+)/){
			$namep = $1;
			push(@primername, $namep);
		}
		elsif($_ =~ /^[a-z]/i){
			$_ =~ tr/[a-z]/[A-Z]/;
			$forward{$namep} = $_;
		}
	}
	if($countp == 2){
		if($_ =~ /^>(.+)/){$namep = $1;}
		elsif($_ =~ /^[a-z]/i){
			$_ =~ tr/[a-z]/[A-Z]/;
			$reverse{$namep} = $_;
		}
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

#primer_set
open (OUT1, ">", ".\/DataBase\/forward.txt") or die("error:$!");
open (OUT2, ">", ".\/DataBase\/reverse.txt") or die("error:$!");
print OUT1 "Forward\n";
foreach(@primername){
	unless($reverse{$_}){
		print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
		exit;
	}
	print OUT1 ">$_\n$forward{$_}\n";
}
print OUT1 "Reverse\n";
foreach(@primername){
	unless($forward{$_}){
		print "Error: Paired primers must have the same name. Please check the primer file: $primer.\n";
		exit;
	}
	print OUT1 ">$_\n\.\n";
}

print OUT2 "Forward\n";
foreach(@primername){
	print OUT2 ">$_\n$reverse{$_}\n";
}
print OUT2 "Reverse\n";
foreach(@primername){
	print OUT2 ">$_\n\.\n";
}
close(OUT1);
close(OUT2);

my $primer_num = keys %forward;
if($primer_num == 1){$separate = 0;}
if($trim and $primer_num == 1){
	foreach(keys %forward){
		$primerF = length($forward{$_}); 
		$primerR = length($reverse{$_});
	}
}elsif($trim and $primer_num > 1){
	print "Error: The setting \"MaxDiff = length\" in Setting.txt can only be used when you select a single primer pair.\n";
	exit; 
}
if($trim){$separate = 0;}


unless(-d ".\/Results"){mkdir ".\/Results";}
mkdir ".\/Results\/1_1_Primer_Trimmed_fastq";

#log
open (LOG, ">", ".\/Results\/log.txt") or die("error:$!");
my (%log, %raw);

my @youbi = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat');
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
if($min < 10){$min = "0$min";}
if($hour < 10){$hour = "0$hour";}
if($mon < 10){$mon = "0$mon";}
if($mday < 10){$mday = "0$mday";}
print LOG "Start: $year/$mon/$mday ($youbi[$wday]) $hour:$min\n";
close(LOG);


print "============================================================\n";
print "                 1_Primer Trimming for DADA2                \n";
print "============================================================\n";
print "Now processing...\n";
foreach(sort @read){
	my $file = $_;
	if($_ =~ /_R\d/){
		my ($filename, $comand);
		$_ =~ /(.+)_R\d/;
		$filename = $1;
		if($file =~ /_R1/){
			if($trim){
				if($vv){
					$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file\" -db \"./DataBase/forward.txt\" -o \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R1_stripped.fq\" -length";
				}else{
					$comand = ".\/Tools\/$usearch -fastx_truncate \".\/Run\/$file\" -stripleft $primerF -stripright 0 -fastqout \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R1_stripped.fq\"";
				}
			}else{
				$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file\" -db \"./DataBase/forward.txt\" -diff $diff -o \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R1_stripped.fq\" -s $separate -l -shift $shift";
			}
			system $comand;
		}else{
			if($trim){
				if($vv){
					$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file\" -db \"./DataBase/reverse.txt\" -o \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R2_stripped.fq\" -length";
				}else{
					$comand = ".\/Tools\/$usearch -fastx_truncate \".\/Run\/$file\" -stripleft $primerR -stripright 0 -fastqout \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R2_stripped.fq\"";
				}
			}else{
				$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file\" -db \"./DataBase/reverse.txt\" -diff $diff -o \".\/Results\/1_1_Primer_Trimmed_fastq\/${filename}_R2_stripped.fq\" -s $separate -l -shift $shift";
			}
			system $comand;
		}
	}
}

#Get data names
undef(%file);
undef(@read);
opendir (DIR, ".\/Results\/1_1_Primer_Trimmed_fastq") or die ("error:$!");
@read = readdir DIR;
foreach (@read) {
	unless($_ =~ /_log.txt/){
		if ($_ =~ /(.+)_R1_.*stripped/){$file{$1}++; $read_count++;}
		if ($_ =~ /(.+)_R2_.*stripped/){$file{$1}++; $read_count++;}
	}
}
closedir DIR;
unless(%file){print "Error: None of fastq files in Run file.\n"; exit;}
foreach(keys %file){
	if($file{$_} < 2){
		print "Error: Missing paired-end files in the Run directory! Both R1 and R2 fastq files are required.\n";
		print "Missing either R1 or R2 fastq file for sample: $_.\n";
		exit;
	}
}


foreach(sort @read){
	unless($_ =~ /_log.txt/){
		my $file = $_;
		my $file2 = $_;
		if($file =~ /_R1_/){
			#R1 input
			open (R1, "<", ".\/Results\/1_1_Primer_Trimmed_fastq\/$file") or die ("error:$!");
			my (%fastq1, $name);
			my $count = 0;
			while(<R1>){
				chomp($_);
				$_ =~ s/\r//;
				if($count == 0){
					$name = $_; $count++;
				}elsif($count == 1){
					$fastq1{$name} = "$_\n"; $count++;
				}elsif($count == 2){
					$fastq1{$name} .= "$_\n";$count++;
				}elsif($count == 3){
					$fastq1{$name} .= "$_\n"; $count = 0;
				}
			}
			close(R1);
			unlink(".\/Results\/1_1_Primer_Trimmed_fastq\/$file");
			
			#R2 input
			$file2 =~ s/_R1_/_R2_/;
			open (R2, "<", ".\/Results\/1_1_Primer_Trimmed_fastq\/$file2") or die ("error:$!");
			my %fastq2;
			$count = 0;
			while(<R2>){
				chomp($_);
				$_ =~ s/\r//;
				if($count == 0){
					$_ =~ s/2:N/1:N/;
					$name = $_; $count++;
				}elsif($count == 1){
					$fastq2{$name} = "$_\n"; $count++;
				}elsif($count == 2){
					$fastq2{$name} .= "$_\n";$count++;
				}elsif($count == 3){
					$fastq2{$name} .= "$_\n"; $count = 0;
				}
			}
			close(R2);
			unlink(".\/Results\/1_1_Primer_Trimmed_fastq\/$file");
			
			#output
			open (OUT1, ">", ".\/Results\/1_1_Primer_Trimmed_fastq\/$file") or die ("error:$!");
			open (OUT2, ">", ".\/Results\/1_1_Primer_Trimmed_fastq\/$file2") or die ("error:$!");
			foreach(sort keys %fastq1){
				if($fastq2{$_}){
					print OUT1 "$_\n$fastq1{$_}";
					my $temp = $_;
					$_ =~ s/1:N/2:N/;
					print OUT2 "$_\n$fastq2{$temp}";
				}
			}
			close(OUT1);
			close(OUT2);
		}
	}
}

unlink(".\/DataBase\/forward.txt");
unlink(".\/DataBase\/reverse.txt");

