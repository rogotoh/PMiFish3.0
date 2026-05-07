#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;

#version 2.1 (2025/12/13) each genus ok
#version 2.0 (2025/02/19) each family ok
#version 1.2 (2019/02/13)


#Setting
my ($db, $family);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^DB\s*=\s*([\w\W]+)/){$db = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
}
close(SET);
my @dblist;
unless($db){print "Error!! Check the DB name in Setting.txt.\n"; exit;}
if($db =~ /\s/){
	$db =~ s/\s+/ /g;
	@dblist = split(/\s/, $db);
}else{push(@dblist, $db);}
foreach(@dblist){
	unless(-f ".\/DataBase\/$_"){print "Error: Database file not found or database name mismatch in Setting.txt.\n"; exit;}
}
unless(-f ".\/Dictionary\/$family"){
	if($family =~ /^no$/i){undef($family);}
	else{print "The file of Dictionary for Family name isn't or don't match the name in Setting.txt\n"; exit;}
}

my @primer;
opendir (DIR, ".\/Results\/5_2_Summary_Table") or die ("error:$!");
my @database = readdir DIR;
my $dbcheck = 0;
foreach(@database){
	if ($_ =~ /Representative_seq_(.+)_with_family_name.fas/){push(@primer, $1); $dbcheck++;}
	if ($_ =~ /Representative_seq_with_family_name.fas/){$dbcheck++;}
}
closedir DIR;

unless($dbcheck){print "This scripts need \"Representative_seq_with_family_name.fas\" in 5_2_Summary_Table folder.\n"; exit;}

print "============================================================\n";
print "          5_4_Phylogenetic Analysis with Database           \n";
print "============================================================\n";

my $pri_n = 1;
my $primer;
if(@primer){
	print "Primer List:\n";
	foreach(@primer){
		print " $pri_n $_\n";
		$pri_n++;
	}
	print "Please select the number:";
	$primer = <STDIN>;
	chomp($primer);
	$primer =~ s/\r//;
	if($primer > @primer or $primer == 0){print "Unacceptable number\n"; exit;}
	$primer = $primer[$primer - 1];
}

#import_family_list
my %family;
open(FILE, "<", ".\/Dictionary\/$family") or die ("error:$!");
while(<FILE>){
	chomp($_);
	my @temp = split(/\t/, $_);
	unless($temp[1]){next;}
	unless($temp[1] =~ /^cf$/i){$family{$temp[1]} = $temp[0];}
	else{$temp[2] =~ /([^_]+)/; $family{$1} = $temp[0];}
}
close(FILE);

mkdir ".\/Results\/5_4_PA_with_DB";
my ($familyn, $dataname, %data, @familyn, %fname, %genus, $genus);
if($primer){
	open(FILE, "<", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/${primer}_merged_seq_with_family_name.fas") or die ("error:$!");
}else{
	open(FILE, "<", ".\/Results\/5_1_Fasta_for_Phylogenetic_Analysis\/merged_seq_with_family_name.fas") or die ("error:$!");
}
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ />N_A_/){last;}
	if($_ =~ /^>([^_]+_[^_]+)_([^_]+)/){
		$familyn = $1;
		$genus = $2;
		if($_ =~ /_U[^_]+_([^_]+)/){
			print "$1\n";
			$genus = $1;
		}
		$dataname = $_;
		$fname{$familyn}++;
		$genus{$familyn}{$genus}++;
		if($fname{$familyn} == 1){push(@familyn, $familyn);}
	}else{$data{$familyn}{$dataname} = $_;}
}
close(FILE);
undef($familyn); undef($dataname); undef($genus);

my ($fname, $std);
$pri_n = 1;
if(@familyn){
	#DB
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
	
		print "\nDataBase List:\n";
		foreach(@dblist){
			print " $pri_n $_\n";
			$pri_n++;
		}
		print " $pri_n all\n";
		print "Please select the number:";
		$std = <STDIN>;
		chomp($std);
		$std =~ s/\r//;
		if($std > @dblist + 1 or $std == 0){print "Unacceptable number\n\n"; exit;}
		unless($std == $pri_n){$db = $dblist[$std - 1];}
		else{$db = "merged_DB.fas";}
	}
	$pri_n = 1;
	
	#family
	print "\nFamily List:\n";
	foreach(@familyn){
		print " $pri_n $_\n";
		$pri_n++;
	}
	print " $pri_n all\n";
	print "Please select the number:";
	$std = <STDIN>;
	chomp($std);
	$std =~ s/\r//;
	if($std > @familyn + 1 or $std == 0){print "Unacceptable number\n\n"; exit;}
	unless($std == $pri_n){$fname = $familyn[$std - 1];}
	
	#genus
	if($fname){
		$pri_n = 1;
		print "\nGenus List:\n";
		my @genus = sort {$a cmp $b} keys %{$genus{$fname}};
		foreach(@genus){
			print " $pri_n $_\n";
			$pri_n++;
		}
		print " $pri_n all\n";
		print "Please select the number:";
		$std = <STDIN>;
		chomp($std);
		$std =~ s/\r//;
		if($std > @genus + 1 or $std == 0){print "Unacceptable number\n\n"; exit;}
		unless($std == $pri_n){$genus = $genus[$std - 1];}
	}
}

#Collecting data from database
my (%database);
open(FILE, "<", ".\/Database\/$db") or die ("error:$!");
while(<FILE>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){
		$_ =~ s/__/_/g;
		$_ =~ /.+\|([^_]+)_([^_]+)/;
		my $temp1 = $1;
		my $temp2 = $2;
		if($temp1 and $temp2){
			if($temp1 =~ /^cf$/i){
				unless($family{$temp2}){$familyn = "N_A";}
				else{$familyn = $family{$temp2};}
				$dataname = $_;
			}else{
				unless($family{$temp1}){$familyn = "N_A";}
				else{$familyn = $family{$temp1};}
				$dataname = $_;
			}
		}
	}else{
		$database{$familyn}{$dataname} = $_;
	}
}
close(FILE);

#output
my $fas;
if($primer){$fas = "_${primer}.fas";}
else{$fas = ".fas";}
my $temp_family;
foreach(sort keys %data){
	my $count = 0;
	my $familyname = $_;
	if($fname){
		unless($familyname =~ /$fname/){next;}
	}
	
	if($genus){
		open(FILE, ">", ".\/Results\/5_4_PA_with_DB\/${genus}$fas") or die ("error:$!");
	}else{
		open(FILE, ">", ".\/Results\/5_4_PA_with_DB\/${familyname}$fas") or die ("error:$!");
	}
	my $temp = $data{$_};
	foreach(sort keys %$temp){
		my $filename = $_;
		my $tempname = $_;
		$tempname =~ s/\'//g;
		$tempname =~ s/,//g;
		$tempname =~ s/&//g;
		if($genus){
			if($tempname =~ /_${genus}_/){print FILE "$tempname\n$data{$familyname}{$filename}\n";  $count++;}
		}else{
			print FILE "$tempname\n$data{$familyname}{$filename}\n"; $count++;
		}
	}
	my $temp2 = $database{$familyname};
	
	#Remove duplicates that share both the same name and the same sequence.
	my (@check, %check);
	foreach(sort keys %$temp2){
		my $filename = $_;
		my $cutname = $_;
		$cutname =~ s/.+\|//;
		unless($check{$database{$familyname}{$filename}}){$check{$database{$familyname}{$filename}}{$cutname}++; push(@check, $filename);}
		else{
			unless($check{$database{$familyname}{$filename}}{$cutname}){$check{$database{$familyname}{$filename}}{$cutname}++; push(@check, $filename);}
		}
	}
	foreach(sort @check){
		my $filename = $_; 
		my $tempname = $_;
		$tempname =~ s/\'//g;
		$tempname =~ s/,//g;
		$tempname =~ s/&//g;
		if($genus){
			if($tempname =~ /${genus}_/){print FILE "$tempname\n$database{$familyname}{$filename}\n"; $count++;}
		}else{
			print FILE "$tempname\n$database{$familyname}{$filename}\n"; $count++;
		}
	}
	close(FILE);
	
	if($genus){$familyname = $genus;}
	
	my $length = length($familyname);
	$length = 30 - $length;
	
	my $name = "${familyname}$fas";
	for(my $v = 0; $v < $length;$v++){$name = $name." ";}
	print "\n$name$count sequences\n";
}

print "\n============================================================\n";

#Get data names
opendir (DIR, ".\/Results\/5_4_PA_with_DB") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+).fas/){$file{$1}++;}
}
closedir DIR;

my @check;
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/5_4_PA_with_DB\/${file}.fas") or die("error:$!");
	my $count = 0;
	while(<DATA>){
		if($_ =~ /^>/){$count++;}
	}
	if($count > 3){push(@check, $file);}
}


#Phylogenetic_trees�̍쐬
if(-f ".\/Tools\/muscle_align_nucleotide.mao" and -f ".\/Tools\/infer_NJ_nucleotide.mao"){
	foreach (@check){
		my $file = $_;
		if(-f ".\/Results\/5_4_PA_with_DB\/${file}.nwk"){next;}
		my $command = "megacc12 -a \".\/Tools\/muscle_align_nucleotide.mao\" -d \".\/Results\/5_4_PA_with_DB\/${file}.fas\" -o $file";
		system $command;
		move (".\/${file}.meg", ".\/Results\/5_4_PA_with_DB");
		unlink ".\/${file}_summary.txt";
		
		$command = "megacc12 -a \".\/Tools\/infer_NJ_nucleotide.mao\" -d \".\/Results\/5_4_PA_with_DB\/${file}.meg\" -o $file";
		system $command;
		move (".\/${file}.nwk", ".\/Results\/5_4_PA_with_DB");
		unlink ".\/${file}_consensus.nwk";
		unlink ".\/${file}_summary.txt";
		unlink ".\/${file}_partitions.txt";
	}
}

