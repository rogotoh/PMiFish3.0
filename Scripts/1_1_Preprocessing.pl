#!/usr/bin/perl
use strict;
use warnings;

#2025/02/07 Usearch Options
#2018/12/04 1_2 Strip_primers change to use the Primer_Cleaner.pl
#2018/07/01 1_2_Strip_primers.pl were changed to permit no primer seq
#2018/01/29 change Blast to Usearch_global
#2018/01/20 output a log file at 1_1 step
#2018/01/05 add compress option

#Setting
my ($read_type, $db, $primer, $diff, $shift, $separate, $depth, $size, $identity, $identity2, $dictionary, $family, $temporary, $compress, $trim);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^DB\s*=\s*([\w\W]+)/){$db = $1;}
	elsif($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^MaxDiff\s*=\s*(\S+)/){$diff = $1;}
	elsif($_ =~ /^Shift\s*=\s*(\S+)/){$shift = $1;}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
	elsif($_ =~ /^Depth\s*=\s*(\S+)/){$depth = $1;}
	elsif($_ =~ /^Length\s*=\s*(\d+)/){$size = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Common_name\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
	elsif($_ =~ /^Temporary\s*=\s*(\S+)/){$temporary = $1;}
	elsif($_ =~ /^Compress\s*=\s*(\S+)/){$compress = $1;}
}
close(SET);

$read_type = 2; #paired-end

if($diff =~ /^length$/i){$trim = 1;}
unless($diff =~ /\d+/){$diff = 2;}
unless($shift =~ /\d+/){$shift = 2;}
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($db){print "Error: Please check the database name in Setting.txt.\n"; exit;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($depth){print "Error: Please check the Depth value in Setting.txt.\n"; exit;}
unless($size){print "Error: Please check the Length value in Setting.txt.\n"; exit;}
unless($identity){print "Error: Please check the upper threshold for homology search in Setting.txt.\n"; exit;}
unless($identity2){print "Error: Please check the lower threshold for homology search in Setting.txt.\n"; exit;}
unless($temporary){$temporary = "YES";}
unless($compress){$compress = "NO";}
my @dblist;
if($db =~ /\s/){
	$db =~ s/\s+/ /g;
	@dblist = split(/\s/, $db);
}else{push(@dblist, $db);}
foreach(@dblist){
	unless(-f ".\/DataBase\/$_"){print "Error: Database file not found or database name mismatch in Setting.txt.\n"; exit;}
}
unless($primer =~ /^no$/i){
	unless(-f ".\/DataBase\/$primer"){print "Error: Primer file not found or primer file name mismatch in Setting.txt.\n"; exit;}
}
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "Error: Dictionary file for common names not found or name mismatch in Setting.txt.\n"; exit;}
}
unless(-f ".\/Dictionary\/$family"){
	if($family =~ /^no$/i){undef($family);}
	else{print "Error: Dictionary file for family names not found or name mismatch in Setting.txt.\n"; exit;}
}

#Usearch options
my ($trunctail, $minlen, $maxdiffs, $pctid, $minovlen, $maxee);
open (SET, "<", "Options_usearch.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^\t\*-fastq_trunctail\s*(\d+)/){$trunctail = $1;}
	elsif($_ =~ /^\t\*-fastq_minlen\s*(\d+)/){$minlen = $1;}
	elsif($_ =~ /^\t\*-fastq_maxdiffs\s*(\d+)/){$maxdiffs = $1;}
	elsif($_ =~ /^\t\*-fastq_pctid\s*(\d+)/){$pctid = $1;}
	elsif($_ =~ /^\t\*-fastq_minovlen\s*(\d+)/){$minovlen = $1;}
	elsif($_ =~ /^\t\*-fastq_maxee\s*(\S+)/){$maxee = $1;}
}
close(SET);
unless($trunctail){$trunctail = 2;}
unless($minlen){$minlen = 64;}
unless($maxdiffs){$maxdiffs = 5;}
unless($pctid){$pctid = 90;}
unless($minovlen){$minovlen = 16;}
unless($maxee){$maxee = 1.0;}

#Get primer sequences
my $countp = 0;
my (%forward, %reverse, $namep);
my ($primerF, $primerR);

unless($primer =~ /^no$/i){
	open (DATA2, "<", ".\/DataBase\/$primer") or die("error:$!");
	while(<DATA2>){
		if($_ =~ /\#/){next;}
		$_ =~ s/\r//g;
		chomp($_);
		if($_ =~ /Forward/){$countp = 1;next;}
		if($_ =~ /Reverse/){$countp = 2;next;}
		if($countp == 1){
			if($_ =~ /^>(.+)/){$namep = $1;}
			elsif($_ =~ /^[a-z]/i){$_ =~ tr/[a-z]/[A-Z]/; $forward{$namep} = $_;}
		}
		if($countp == 2){
			if($_ =~ /^>(.+)/){$namep = $1;}
			elsif($_ =~ /^[a-z]/i){
				my $rev = reverse($_);
				$rev =~ tr/ATGCURYMKDHBV/TACGAYRKMHDVB/;
				$reverse{$namep} = $rev;}
			elsif($_ =~ /^\./){$reverse{$namep} = $_;}
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
	my $primer_num = keys %forward;
	if($primer_num == 1){$separate = 0;}
	if($trim and $primer_num == 1){
		foreach(keys %forward){
			$primerF = length($forward{$_}); 
			if($reverse{$_} =~ /^\./){$primerR = 0;}
			else{$primerR = length($reverse{$_});}
		}
	}elsif($trim and $primer_num > 1){
		print "Error: The setting \"MaxDiff = length\" in Setting.txt can only be used when you select a single primer pair.\n";
		exit; 
	}
	if($trim){$separate = 0;}
}

#Database and Usearch check
if(@dblist > 1){
	unless(-f ".\/DataBase\/merged_DB.fas"){
		open (OUT, ">", ".\/DataBase\/merged_DB.fas") or die $!;
		foreach(@dblist){
			my $temp = $_;
			open (DATA, "<", ".\/DataBase\/$temp") or die $!;
			print OUT $_ while(<DATA>);
			close(DATA);
		}
		close(OUT);
	}
	$db = "merged_DB.fas";
}

opendir (DIR, ".\/DataBase") or die ("error:$!");
my @database = readdir DIR;
my $dbcheck = 0;
foreach(@database){
	if ($_ =~ /${db}\.udb/){$dbcheck++;}
}
closedir DIR;

opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv, $v12);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){
		$usearch = $1;
		if($_ =~ /12/){$v12 = 1;}
		last;
	}elsif($_ =~ /(vsearch.*)/){
		$usearch = $1; $vv = 1;
	}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}

unless($dbcheck){
	if($vv){
		my $command = ".\/Tools\/$usearch --makeudb_usearch \".\/DataBase\/$db\" --output \".\/DataBase\/${db}\.udb\"";
		system $command;
	}else{
		my $command = ".\/Tools\/$usearch -makeudb_usearch \".\/DataBase\/$db\" -output \".\/DataBase\/${db}\.udb\"";
		system $command;
	}
}

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
	if($read_type == 2 and $file{$_} < 2){
		print "Error: Missing paired-end files in the Run directory! Both R1 and R2 fastq files are required.\n";
		print "Missing either R1 or R2 fastq file for sample: $_.\n";
		exit;
	}
}
unless(-d ".\/Results"){mkdir ".\/Results";}


print "============================================================\n";
print "                      1_Preprocessing                       \n";
print "============================================================\n";
#assemble pair seqs
mkdir ".\/Results\/1_1_Merge_paird_reads";
unless($primer =~ /^no$/i){mkdir ".\/Results\/1_2_Strip_primers";}
mkdir ".\/Results\/1_3_Quality_filter";
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

foreach(sort keys %file){
	my $file = $_;
	my ($file1, $comand);
	foreach(sort @read){
		if($_ =~ /${file}_R1[\._]/){$file1 = $_;}
	}
	my $raw_read = 0;
	if($read_type == 1){#single read
		open (DATA, "<", ".\/Run\/$file1") or die("error:$!");
		while(<DATA>){
			if($_ =~ /^\@/){$raw_read++;}
		}
		close(DATA);
		$log{$file} = $raw_read;
		$raw{$file} = $raw_read;
		$log{$file} = $log{$file} . "\t-";
	}else{#paried read
		print "============================================================\n";
		print "$file    1_1_Merge_paird_reads                   \n";
		print "============================================================\n";
		if($vv){
			my $file2 = $file1;
			$file2 =~ s/_R1_/_R2_/;
			$comand = ".\/Tools\/$usearch --fastq_mergepairs \".\/Run\/$file1\" --reverse \".\/Run\/$file2\" --fastqout \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" --log \".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt\" --fastq_minlen $minlen --fastq_maxdiffs $maxdiffs --fastq_minovlen $minovlen";
		}else{
			$comand = ".\/Tools\/$usearch -fastq_mergepairs \".\/Run\/$file1\" -fastqout \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -log \".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt\" -fastq_trunctail $trunctail -fastq_minlen $minlen -fastq_maxdiffs $maxdiffs -fastq_pctid $pctid -fastq_minovlen $minovlen";
		}
		system $comand;
		if($compress =~ /yes/i){
			print "\nCompressing fastq file to gz file...\n";
			foreach(@read){
				if($_ !~ /gz$/ and $_ =~ /$file/){
					$comand = "gzip \".\/Run\/$_\"";
					system $comand;
				}
			}
		}
		if(-f ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt"){
			open (DATA, "<", ".\/Results\/1_1_Merge_paird_reads\/${file}_log.txt") or die("error:$!");
			while(<DATA>){
				chomp($_);
				$_ =~ s/\r//g;
				if($_ =~ /(\d+)\s+Pairs/){$log{$file} = $1; $raw{$file} = $1;}
				elsif($_ =~ /(\d+)\s+Read pairs/){$log{$file} = $1; $raw{$file} = $1;}
				if($vv){
					if($_ =~ /(\d+)\s+Merged.+\s*\((\S+)%/){$log{$file} = $log{$file} . "\t$1 ($2%)";}
				}else{
					if($_ =~ /(\d+)\s+Merged.+\s([^\s]+)%/){$log{$file} = $log{$file} . "\t$1 ($2%)";}
				}
				unless($log{$file}){$log{$file} = "0\t0 (0.00%)";}
			}
			close(DATA);
		}
		print "\n\n";
	}
	
	print "============================================================\n";
	print "$file    1_2_Strip_primers                      \n";
	print "============================================================\n";
	unless($primer =~ /^no$/i){
		print "Now processing...\n";
		if($read_type == 1){#single read
			if($trim){
				if($vv or $v12){
					$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file1\" -db \"./DataBase/$primer\" -o \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -length";
				}else{
					$comand = ".\/Tools\/$usearch -fastx_truncate \".\/Run\/$file1\" -stripleft $primerF -stripright $primerR -fastqout \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\"";
				}
			}else{
				$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Run\/$file1\" -db \"./DataBase/$primer\" -diff $diff -o \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -s $separate -l -shift $shift";
			}
			system $comand;
		}else{#paried read
			if($trim){
				if($vv or $v12){
					$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -db \"./DataBase/$primer\" -o \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -length";
				}else{
					$comand = ".\/Tools\/$usearch -fastx_truncate \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -stripleft $primerF -stripright $primerR -fastqout \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\"";
				}
			}else{
				$comand = "perl ./Scripts/Primer_Cleaner.pl -f \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -db \"./DataBase/$primer\" -diff $diff -o \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -s $separate -l -shift $shift";
			}
			system $comand;
		}
		if($temporary =~ /yes/i){unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq";}
	}else{
		print "The step was skipped (You selected \"Primers = No\" in Setting.txt).\n";
	}
	if($primer !~ /^no$/i){
		if($separate){
			opendir (DIR, ".\/Results\/1_2_Strip_primers\/") or die ("error:$!");
			my @reads = readdir DIR;
			my %divide;
			foreach (@reads) {
				if ($_ =~ /(${file}.+)_stripped.fq/){$divide{$1}++;}
			}
			closedir DIR;
			foreach(sort keys %divide){
				my $files = $_;
				print "\n\n";
				print "============================================================\n";
				print "$files    1_3_Quality_filter                      \n";
				print "============================================================\n";
				if($vv){
					$comand = ".\/Tools\/$usearch --fastq_filter \".\/Results\/1_2_Strip_primers\/${files}_stripped.fq\" --fastq_maxee $maxee --fastq_minlen $size --fastaout \".\/Results\/1_3_Quality_filter\/${files}_filtered.fa\" --log \".\/Results\/1_3_Quality_filter\/${files}_log.txt\"";
				}else{
					$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_2_Strip_primers\/${files}_stripped.fq\" -fastq_maxee $maxee -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${files}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${files}_log.txt\"";
				}
				system $comand;
				my $slog;
				open (DATA, "<", ".\/Results\/1_2_Strip_primers\/${files}_log.txt") or die("error:$!");
				while(<DATA>){
					if($_ =~ />\s+(\d+)\s+reads/){$slog = $1;}
				}
				close(DATA);
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${files}_log.txt") or die("error:$!");
				my $log_check = 0;
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($vv){
						if($_ =~ /(\d+)\s+sequences\s+kept/){
							my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$files} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
							$log_check = 1;
							last;
						}
					}else{
						if($_ =~ /(\d+)\s+Filtered/){
							my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$files} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
							$log_check = 1;
							last;
						}
					}
				}
				close(DATA);
				unless($log_check){
					if($slog){
						my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
						$log{$files} = $log{$file} . "\t$slog ($sper%)\t0 (0.00%)";
					}else{
						$log{$files} = $log{$file} . "\t0 (0.00%)\t0 (0.00%)";
					}
				}
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${files}_stripped.fq";}
			}
			undef($log{$file});
		}elsif(-f ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq"){
			print "\n\n";
			print "============================================================\n";
			print "$file    1_3_Quality_filter                      \n";
			print "============================================================\n";
			if($vv){
				$comand = ".\/Tools\/$usearch --fastq_filter \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" --fastq_maxee $maxee --fastq_minlen $size --fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" --log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
			}else{
				$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_2_Strip_primers\/${file}_stripped.fq\" -fastq_maxee $maxee -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
			}
			system $comand;
			unless($trim){
				my $slog;
				open (DATA, "<", ".\/Results\/1_2_Strip_primers\/${file}_log.txt") or die("error:$!");
				while(<DATA>){
					if($_ =~ />\s+(\d+)\s+reads/){$slog = $1;}
				}
				close(DATA);
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
				my $log_check = 0;
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($vv){
						if($_ =~ /(\d+)\s+sequences\s+kept/){
							if($1 == 0){next;}
							my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$file} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
							$log_check = 1;
							last;
						}
					}else{
						if($_ =~ /(\d+)\s+Filtered/){
							my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$file} = $log{$file} . "\t$slog ($sper%)\t$1 ($per%)";
							$log_check = 1;
							last;
						}
					}
				}
				close(DATA);
				unless($log_check){
					if($slog){
						my $sper = sprintf("%.2f", $slog/$raw{$file}*100);
						$log{$file} = $log{$file} . "\t$slog ($sper%)\t0 (0.00%)";
					}else{
						$log{$file} = $log{$file} . "\t0 (0.00%)\t0 (0.00%)";
					}
				}
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq";}
			}else{
				open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
				my $log_check = 0;
				while(<DATA>){
					chomp($_);
					$_ =~ s/\r//g;
					if($vv){
						if($_ =~ /(\d+)\s+sequences\s+kept/){
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$file} = $log{$file} . "\t$1 ($per%)";
							$log_check = 1;
						}
					}else{
						if($_ =~ /(\d+)\s+Filtered/){
							my $per = sprintf("%.2f", $1/$raw{$file}*100);
							$log{$file} = $log{$file} . "\t$1 ($per%)";
							$log_check = 1;
						}
					}
				}
				close(DATA);
				unless($log_check){
					$log{$file} = $log{$file} . "\t0 ($0.00%)";
				}
				if($temporary =~ /yes/i){unlink ".\/Results\/1_2_Strip_primers\/${file}_stripped.fq";}
			}
		}else{
			unless($trim){
				$log{$file} = $log{$file} . "\t0 (0.00%)\t0 (0.00%)";
			}else{
				$log{$file} = $log{$file} . "\t0 ($0.00%)";
			}
		}
	}elsif($primer =~ /^no$/i and -f ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq"){
		print "\n\n";
		print "============================================================\n";
		print "$file    1_3_Quality_filter                      \n";
		print "============================================================\n";
		if($vv){
			$comand = ".\/Tools\/$usearch --fastq_filter \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" --fastq_maxee $maxee --fastq_minlen $size --fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" --log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
		}else{
			$comand = ".\/Tools\/$usearch -fastq_filter \".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq\" -fastq_maxee $maxee -fastq_minlen $size -fastaout \".\/Results\/1_3_Quality_filter\/${file}_filtered.fa\" -log \".\/Results\/1_3_Quality_filter\/${file}_log.txt\"";
		}
		system $comand;
		open (DATA, "<", ".\/Results\/1_3_Quality_filter\/${file}_log.txt") or die("error:$!");
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($vv){
				if($_ =~ /(\d+)\s+sequences\s+kept/){
					my $per = sprintf("%.2f", $1/$raw{$file}*100);
					$log{$file} = $log{$file} . "\t$1 ($per%)";
				}else{
					$log{$file} = $log{$file} . "\t0 (0.00%)";
				}
			}else{
				if($_ =~ /(\d+)\s+Filtered/){
					my $per = sprintf("%.2f", $1/$raw{$file}*100);
					$log{$file} = $log{$file} . "\t$1 ($per%)";
				}else{
					$log{$file} = $log{$file} . "\t0 (0.00%)";
				}
			}
		}
		close(DATA);
		if($temporary =~ /yes/i){unlink ".\/Results\/1_1_Merge_paird_reads\/${file}_assembled_seq.fq";}
	}
	print "\n\n";
}
foreach(sort keys %log){
	if($log{$_}){print LOG "$_\t$log{$_}\n";}
}
close(LOG);
