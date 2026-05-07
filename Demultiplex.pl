#!/usr/bin/perl
use strict;
no strict "refs";
use warnings;

#option setting
my ($scorec, $score, $kensaku, $pair, $unzip);
for (my $v = 0; $v < @ARGV; $v++){
	if($ARGV[$v] =~ /-s$/){$scorec = 1; $score = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-t$/){$pair = $ARGV[$v+1];}
	elsif($ARGV[$v] =~ /-gz$/){$unzip = 1;}
}
unless($pair){$pair = 2;}
unless($unzip){$unzip = 2;}
if($pair == 1 or $pair == 2){}
else{print "-t can only use 1(Single-End) or 2(Pair-End).\n"; exit;}

unless($scorec){$score = 20;}
unless($score =~ /\d+/){print "-s option error!\n"; exit;}
if($score == 0 or $score == 10 or $score == 20 or $score == 30){}
else{print "-s can only use 0 or 10 or 20 or 30.\n"; exit;}

if($score == 0){$kensaku = "a";}
if($score == 10){$kensaku = "\!\"\#\$\%\&\'\(\)\*";}
if($score == 20){$kensaku = "\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/01234";}
if($score == 30){$kensaku = "\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>";}

#Index data
open (DATA, "<", ".\/Run\/non-demultiplexed\/index.csv") or die("error:$!");
my (%index, %r2_check, $index_length, @sample);
while(<DATA>){
	$_ =~ s/\r//g;
	$_ =~ s/\"//g;
	$_ =~ s/ /_/g;
	chomp($_);
	unless($_ =~ /^Sample_Name,/){
		my @temp = split (/,/, $_);
		if($pair == 1){
			$index{$temp[1]} = $temp[0];
			push(@sample, "$temp[0]_R1");
		}else{
			$index{$temp[1]}{$temp[2]} = $temp[0];
			push(@sample, "$temp[0]_R1");
			push(@sample, "$temp[0]_R2");
			$r2_check{$temp[2]}++;
		}
		unless($index_length){$index_length = length($temp[1]);}
	}
}
close(DATA);
if($pair == 1){push(@sample, "Low_quality_QC${score}_R1");}
push(@sample, "Undetermind_QC${score}_R1");
if($pair == 2){push(@sample, "Undetermind_QC${score}_R2");}

#Decompressing files
opendir (DIR, ".\/Run\/non-demultiplexed") or die ("error:$!");
my @run = readdir DIR;
my $gzip = 0;
foreach (@run) {
	if ($_ =~ /\.gz$/){$gzip = 1;}
}
closedir DIR;

#Get data names
my $read_type = 2;
opendir (DIR, ".\/Run\/non-demultiplexed") or die ("error:$!");
my @read = readdir DIR;
my $read_count = 0;
my %file;
my @filename;
foreach (@read) {
	if ($_ =~ /(.+)_R1/){$file{$1}++; $read_count++; push(@filename, $_);}
	if ($_ =~ /(.+)_R2/){$file{$1}++; $read_count++; push(@filename, $_);}
}
closedir DIR;
unless(%file){print "Error: None of fastq files in Run file.\n"; exit;}
foreach(keys %file){
	if($pair == 2 and $file{$_} < 2){
		print "Error: Missing paired-end files in the Run directory! Both R1 and R2 fastq files are required.\n";
		print "Missing either R1 or R2 fastq file for sample: $_.\n";
		exit;
	}
}
if($pair == 1){
	unless($read_count == 1){print "Error: When Single-End was selected, only one fastq file can demultiplex.\n"; exit;}
}
if($pair == 2){
	unless($read_count == 2){print "Error: When Pair-End was selected, only 2 fastq files (R1 and R2) can demultiplex.\n"; exit;}
}

#demultiplex
my (%open_hundle, %read_c);
foreach(@sample){
	if($gzip and $unzip == 1){open ($open_hundle{$_}, "|-", "gzip > .\/Run\/${_}_001.fastq.gz") or die("error:$!");}
	else{open ($open_hundle{$_}, ">", ".\/Run\/${_}_001.fastq") or die("error:$!");}
}
print "============================================================\n";
print "                        Demultiplex                         \n";
print "============================================================\n";

#Single-End
if($pair == 1){
	print "Now processing...\n\n";
	foreach(sort @filename){
		unless($_ =~ /fastq/){next;}
		print "Demultiplex $_ \n\n";
		my $filename = $_;
		
		if($gzip){open (DATA, "gzip -dc .\/Run\/non-demultiplexed\/$filename |") or die("error:$!");}
		else{open (DATA, "<", ".\/Run\/non-demultiplexed\/$filename") or die("error:$!");}
		my $count = 0;
		my $index_check = 0;
		my ($data, $samplename, $read_name, $index);
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($count == 0){
				$data = "$_\n";
				$count++;
			}elsif($count == 1){
				$index = substr($_, -$index_length, $index_length, "");
				if($index{$index}){
					$index_check = 1;
				}
				$data .= "$_\n";
				$count++;
			}elsif($count == 2){
				$data .= "$_\n";
				$count++;
			}elsif($count == 3){
				$data .= "$_\n";
				my $qc = substr($_, -$index_length, $index_length, "");
				if($index_check){
					unless($qc =~ /[$kensaku]/){
						$samplename = $index{$index};
						$read_c{"${samplename}_R1"}++;
						my $fh = $open_hundle{"${samplename}_R1"};
						print $fh "$data";
					}else{
						my $fh = $open_hundle{"Low_quality_QC${score}_R1"};
						print $fh "$data";
						$read_c{"Low_quality_QC${score}_R1"}++;
					}
				}else{
					my $fh = $open_hundle{"Undetermind_QC${score}_R1"};
					print $fh "$data";
					$read_c{"Undetermind_QC${score}_R1"}++;
				}
				undef($data);
				undef($samplename);
				undef($read_name);
				undef($index);
				$count = 0;
				$index_check = 0;
			}
		}
		close(DATA);
	}
	foreach(@sample){
		unless($read_c{$_}){
			print "${_}_001.fastq\t0 reads\n";
		}else{
			print "${_}_001.fastq\t$read_c{$_} reads\n";
		}
		my $fh = $open_hundle{$_};
		close ($fh);
	}
}

#Pair-End
if($pair == 2){
	my (%gokaku, %low, %unde, %r1_index, %r2_index, @unmulti);
	print "Now processing...\n\n";
	foreach(sort @filename){
		unless($_ =~ /fastq/){next;}
		push(@unmulti, $_);
		my $filename = $_;
		my $read_muki;
		if($_ =~ /_R1_/){
			$read_muki = 1;
		}else{
			$read_muki = 2;
			print "Demultiplex $_\n\n";
		}
		
		if($gzip){open (DATA, "gzip -dc .\/Run\/non-demultiplexed\/$filename |") or die("error:$!");}
		else{open (DATA, "<", ".\/Run\/non-demultiplexed\/$filename") or die("error:$!");}
		my $total = 0;
		my $count = 0;
		my $index_check = 0;
		my ($data, $samplename, $read_name, $index);
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($count == 0){
				if($read_muki == 2){$data = "$_\n";}
				$total++;
				$count++;
			}elsif($count == 1){
				if($read_muki == 1){
					$index = substr($_, -$index_length, $index_length, "");
					if($index{$index}){
						$index_check = 1;
						$r1_index{$total} = $index;
					}
				}elsif($read_muki == 2){
					$index = substr($_, 0, $index_length, "");
					if($r2_check{$index}){
						$index_check = 1;
						if($gokaku{$total}){
							if($index{$r1_index{$total}}{$index}){
								$samplename = $index{$r1_index{$total}}{$index};
							}
						}
					}
					$data .= "$_\n";
				}
				$count++;
			}elsif($count == 2){
				if($read_muki == 2){$data .= "$_\n";}
				$count++;
			}elsif($count == 3){
				my $qc;
				if($read_muki == 1){
					$qc = substr($_, -$index_length, $index_length, "");
					if($index_check){
						unless($qc =~ /[$kensaku]/){
							$gokaku{$total}++;
						}else{
							$low{$total}++;
						}
					}else{
						$low{$total}++;
					}
				}else{
					$qc = substr($_, 0, $index_length, "");
					$data .= "$_\n";
					if($index_check){
						unless($qc =~ /[$kensaku]/){
							if($samplename){
								$r2_index{$total} = $index;
								my $fh = $open_hundle{"${samplename}_R2"};
								print $fh "$data";
								$read_c{"${samplename}_R1"}++;
								$read_c{"${samplename}_R2"}++;
							}else{
								$low{$total}++;
								my $fh = $open_hundle{"Undetermind_QC${score}_R2"};
								print $fh "$data";
								$read_c{"Undetermind_QC${score}_R1"}++;
								$read_c{"Undetermind_QC${score}_R2"}++;
							}
						}else{
							$low{$total}++;
							my $fh = $open_hundle{"Undetermind_QC${score}_R2"};
							print $fh "$data";
							$read_c{"Undetermind_QC${score}_R1"}++;
							$read_c{"Undetermind_QC${score}_R2"}++;
						}
					}else{
						$low{$total}++;
						my $fh = $open_hundle{"Undetermind_QC${score}_R2"};
						print $fh "$data";
						$read_c{"Undetermind_QC${score}_R1"}++;
						$read_c{"Undetermind_QC${score}_R2"}++;
					}
				}
				undef($data);
				undef($samplename);
				undef($read_name);
				undef($index);
				$count = 0;
				$index_check = 0;
			}
		}
		close(DATA);
	}

	#R1 output
	foreach(@unmulti){
		unless($_ =~ /fastq/){next;}
		if($_ =~ /_R2_/){next;}
		print "Demultiplex $_\n\n";
		my $filename = $_;
		my $total = 0;
		my $count = 0;
		if($gzip){open (DATA, "gzip -dc .\/Run\/non-demultiplexed\/$filename |") or die("error:$!");}
		else{open (DATA, "<", ".\/Run\/non-demultiplexed\/$filename") or die("error:$!");}
		my $data;
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($count == 0){$data = "$_\n"; $total++; $count++;}
			elsif($count == 1){
				substr($_, -$index_length, $index_length, "");
				$data .= "$_\n"; $count++;
			}elsif($count == 2){$data .= "$_\n"; $count++;}
			elsif($count == 3){
				substr($_, -$index_length, $index_length, "");
				$data .= "$_\n";
				if($r2_index{$total}){
					my $samplename = $index{$r1_index{$total}}{$r2_index{$total}};
					my $fh = $open_hundle{"${samplename}_R1"};
					print $fh "$data";
				}elsif($low{$total}){
					my $fh = $open_hundle{"Undetermind_QC${score}_R1"};
					print $fh "$data";
				}else{
					my $fh = $open_hundle{"Undetermind_QC${score}_R1"};
					print $fh "$data";
				}
				undef($data);
				$count = 0;
			}
		}
		close(DATA);
	}
	my $tr = 0;
	foreach(sort @sample){
		unless($read_c{$_}){
			print "${_}_001.fastq\t0 reads\n";
		}else{
			print "${_}_001.fastq\t$read_c{$_} reads\n";
			$tr += $read_c{$_};
		}
		my $fh = $open_hundle{$_};
		close ($fh);
	}
	print "Total Reads(R1 + R2) = $tr\n";
}
