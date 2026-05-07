#!/usr/bin/perl
use strict;
use warnings;

#2018/01/05 add qcovs to outfmt
#2025/01/29 add clustering

#Setting
my ($db, $primer, $separate, $identity, $identity2, $denoise, $type, $family);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^DB\s*=\s*([\w\W]+)/){$db = $1;}
	elsif($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^Divide\s*=\s*(\S+)/){$separate = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^Family\s*=\s*(\S+)/){$family = $1;}
}
close(SET);
unless($separate){$separate = 0;}
if($separate =~ /^yes$/i){$separate = 1;}
else{$separate = 0;}
unless($db){print "Error: Please check the database name in Setting.txt.\n"; exit;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($identity){print "Error: Please check the upper threshold for homology search in Setting.txt.\n"; exit;}
unless($identity2){print "Error: Please check the lower threshold for homology search in Setting.txt.\n"; exit;}
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

#Database check
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
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
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


#Usearch_global
mkdir ".\/Results\/3_1_Usearch_global";
$identity2 = $identity2/100;
my ($tempfile, $command);
print "\n";
if(-f ".\/Results\/2_7_LULU\/curated.fa"){
	$tempfile = ".\/Results\/2_7_LULU\/curated.fa";
}elsif(-f ".\/Results\/2_6_OTU_Clustering\/OTU.fa"){
	$tempfile = ".\/Results\/2_6_OTU_Clustering\/OTU.fa";
}elsif(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	$tempfile = ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa";
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	$tempfile = ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa";
}else{
	exit;
}
print "============================================================\n";
print "                     3_1_Usearch_global                     \n";
print "============================================================\n";
if($vv){
	$command = ".\/Tools\/$usearch --usearch_global \"$tempfile\" --db \".\/DataBase\/${db}.udb\" --id $identity2 --maxaccepts 100 --strand plus --userout \".\/Results\/3_1_Usearch_global\/usearch_results.txt\" --userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov";
}else{
	$command = ".\/Tools\/$usearch -usearch_global \"$tempfile\" -db \".\/DataBase\/${db}.udb\" -id $identity2 -maxaccepts 100 -strand plus -userout \".\/Results\/3_1_Usearch_global\/usearch_results.txt\" -userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+qcov";
}
system $command;

print "\n";

#Family
my %family;
unless($family =~ /^no$/i){
	open (DATA, "<", ".\/Dictionary\/$family") or die("error:$!");
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		my @temp = split(/\t/, $_);
		$family{$temp[1]} = $temp[0];
	}
}
close(DATA);

#make table from usearch_results.txt
open (DATA, "<", ".\/Results\/3_1_Usearch_global\/usearch_results.txt") or die("error:$!");
my (%check, %data, %spname, %conserv, %level, %sort, %oriname, %order, %class, %phylum);
my $ff = 0;
my $tsv_c = 0;
unless(%family){$ff = 1;}
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	$_ =~ s/\t+/\t/g;
	my @temp = split(/\t/, $_);
	my $oriname = $temp[1];
	if($temp[1] =~ /^gb\||^id\|/){$temp[1] =~ s/.+\|//;}
	else{$temp[1] =~ s/[^\|]+\|//;}
	if($temp[-1] < 90){next;} #query cover %
	$temp[0] =~ /(.+)\;size=(\d+)/;
	my $id = $1;
	my $num = $2;
	$id =~ /(\d+)/;
	$sort{$1} = $id;
	if($ff == 1){
		if($temp[1] =~ /^k__/){
			$tsv_c =1;
			$temp[1] =~ /p__([^\_]+).*\;c__([^\_]+).*\;o__([^\_]+).*\;f__([^\_]+).*\;g__([^\_]+).*\;/;
			$phylum{$5} = $1;
			$class{$5} = $2;
			$order{$5} = $3;
			$family{$5} = $4;
		}
	}
	unless($check{$temp[0]}){
		$check{$temp[0]}++;
		$data{$id} = "$id\t$temp[1]\t$temp[2]\t$num";
		$spname{$id} = $temp[1];
		$level{$id} = $temp[2];
	}
	if($temp[2] >= $level{$id} - (100-$identity)){
		$temp[1] =~ s/_cluster\d+//;
		if($temp[1] =~ /s__([^\|]+)/){$temp[1] = $1;}
		$temp[1] =~ s/([^_]+_[^_]+)_[^_]+/$1/;
		$conserv{$id}{$temp[1]}++;
		unless($oriname{$id}){$oriname{$id} = "\t$oriname\t$temp[2]\n";}
		else{$oriname{$id} .= "\t$oriname\t$temp[2]\n";}
	}
}
close(DATA);

#Make neighbors.txt
open (OUT, ">", ".\/Results\/3_1_Usearch_global\/neighbors.txt") or die("error:$!");
foreach(sort {$a <=> $b} keys %sort){
	my $tempid = $sort{$_};
	my $tempnum = keys %{$conserv{$tempid}};
	my @split = split(/\t/, $data{$tempid});
	if($tempnum == 1){
		$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tSpecies\t$split[3]";
	}else{
		my (%genus, %f_check, %o_check, %c_check, %p_check);
		foreach(keys %{$conserv{$tempid}}){
			$_ =~ /(^[^_]+)_/;
			$genus{$1}++;
			if($phylum{$1}){$p_check{$phylum{$1}}++;}
			if($class{$1}){$c_check{$class{$1}}++;}
			if($order{$1}){$o_check{$order{$1}}++;}
			if($family{$1}){$f_check{$family{$1}}++;}
			
		}
		if(%genus){
			$tempnum = keys %genus;
			if($tempnum == 1){
				$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tGenus\t$split[3]";
			}else{
				unless(%f_check){
					$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\t>Genus\t$split[3]";print OUT ">$tempid\n$oriname{$tempid}";next;
				}
				$tempnum = keys %f_check;
				if($tempnum == 1){
					$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tFamily\t$split[3]";
				}else{
					unless(%o_check){
						$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\t>Family\t$split[3]";print OUT ">$tempid\n$oriname{$tempid}";next;
					}
					$tempnum = keys %o_check;
					if($tempnum == 1){
						$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tOrder\t$split[3]";
					}else{
						unless(%c_check){
							$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\t>Order\t$split[3]";print OUT ">$tempid\n$oriname{$tempid}";next;
						}
						$tempnum = keys %c_check;
						if($tempnum == 1){
							$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tClass\t$split[3]";
						}else{
							unless(%p_check){
								$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\t>Class\t$split[3]";print OUT ">$tempid\n$oriname{$tempid}";next;
							}
							$tempnum = keys %p_check;
							if($tempnum == 1){
								$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tPhylum\t$split[3]";
							}else{
								$data{$tempid} = "$split[0]\t$split[1]\t$split[2]\tKingdom\t$split[3]";
							}
						}
					}
				}
			}
		}
	}
	print OUT ">$tempid\n$oriname{$tempid}";
}
close(OUT);
		

#table
if(-f ".\/Results\/2_7_LULU\/curated.fa"){
	open (DATA, "<", ".\/Results\/2_7_LULU\/curated_list.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_6_OTU_Clustering\/OTU.fa"){
	open (DATA, "<", ".\/Results\/2_6_OTU_Clustering\/OTU_table.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv"){
	open (DATA, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv"){
	open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv") or die("error:$!");
}
open (OUT1, ">", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment_table.tsv") or die("error:$!");
open (OUT2, ">", ".\/Results\/3_1_Usearch_global\/nohit_list.txt") or die("error:$!");
my $tsv = 0;
my @nohit;
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($tsv == 0){
		print OUT1 "ID\tScientific Name\tIdentity\tLCA Rank\tTotal Reads$_\n";
		$tsv++;
	}else{
		my @temp = split(/\t/, $_);
		if($data{$temp[0]}){
			$_ =~ s/$temp[0]/$data{$temp[0]}/;
			print OUT1 "$_\n";
		}else{
			push(@nohit, $temp[0]);
		}
	}
}
close(DATA);
my $nohit_num = @nohit;
print OUT2 "Nohit = $nohit_num sequences\n";
print "Nohit = $nohit_num sequences\n";
foreach(@nohit){
	print OUT2 "\t$_\n";
}
close(OUT1);
close(OUT2);

#fasta
if(-f ".\/Results\/2_7_LULU\/curated.fa"){
	open (DATA, "<", ".\/Results\/2_7_LULU\/curated.fa") or die("error:$!");
}elsif(-f ".\/Results\/2_6_OTU_Clustering\/OTU.fa"){
	open (DATA, "<", ".\/Results\/2_6_OTU_Clustering\/OTU.fa") or die("error:$!");
}elsif(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	open (DATA, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa") or die("error:$!");
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa") or die("error:$!");
}
open (OUT1, ">", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment.fa") or die("error:$!");
open (OUT2, ">", ".\/Results\/3_1_Usearch_global\/nohit.fa") or die("error:$!");
my $check = 0;
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>/){
		$_ =~ />(.+)\;size/;
		my $temp = $1;
		if($spname{$temp}){
			print OUT1 "$_$spname{$temp}\n";
			$check++;
		}else{
			print OUT2 "$_\n";
		}
	}else{
		if($check){
			print OUT1 "$_\n";
			$check = 0;
		}else{
			print OUT2 "$_\n";
		}
	}
}
close(DATA);
close(OUT1);
close(OUT2);

#divide table
if($separate){
	foreach(sort @forward){
		my $temp_p = $_;
		open (OUT, ">", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment_${_}_table.tsv") or die("error:$!");
		open (DATA, "<", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment_table.tsv") or die("error:$!");
		my $count = 0;
		my %check_c;
		my %fasta;
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			my $narabi;
			if($count == 0){
				for(my $n = 1; $n < @temp; $n++){
					if($temp[$n] =~ /$temp_p/){
						$check_c{$n}++;
					}
				}
			}
			if($count == 0){
				$narabi = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]";
			}else{
				$narabi = "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]";
			}
			my $sample;
			my $total_read = 0;
			for(my $n = 5; $n < @temp; $n++){
				if($check_c{$n}){
					if($count == 0){
						unless($sample){$sample = "\t$temp[$n]";}
						else{$sample .= "\t$temp[$n]";}
					}else{
						$total_read += $temp[$n];
						unless($sample){$sample = "\t$temp[$n]";}
						else{$sample .= "\t$temp[$n]";}
					}
				}
			}
			if($count == 0){
				print OUT "$narabi$sample\n";
			}else{
				if($total_read){
					print OUT "$narabi\t$total_read$sample\n";
					my @temp = split(/\t/, $narabi);
					$fasta{$temp[0]}++;
				}
			}
			$count = 1;
		}
		close(OUT);
		close(DATA);
		
		#divide fasta
		open (OUT, ">", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment_${temp_p}.fa") or die("error:$!");
		open (DATA, "<", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment.fa") or die("error:$!");
		my $seq = 0;
		while(<DATA>){
			chomp($_);
			$_ =~ s/\r//g;
			if($_ =~ /^>/){
				$_ =~ /^>(.+)\;size/;
				if($fasta{$1}){
					print OUT "$_\n";
					$seq = 1;
				}
			}else{
				if($seq){
					print OUT "$_\n";
					$seq = 0;
				}
			}
		}
		close(OUT);
		close(DATA);
	}
}

#Rank tsv
if($tsv_c == 1){
	opendir (DIR, ".\/Results\/3_1_Usearch_global") or die ("error:$!");
	@database = readdir DIR;
	my @tsv;
	foreach(@database){
		chomp($_);
		$_ =~ s/\r//g;
		if ($_ =~ /\.tsv/){
			unless($_ =~ /Rank\.tsv/){push(@tsv, $_);}
		}
	}
	closedir DIR;
	foreach(@tsv){
		my $filename = $_;
		$filename =~ /(.+)\.tsv/;
		my $fname = $1;
		open (DATA, "<", ".\/Results\/3_1_Usearch_global\/$filename") or die("error:$!");
		open (OUT, ">", ".\/Results\/3_1_Usearch_global\/${fname}_Rank.tsv") or die("error:$!");
		my $count = 0;
		while(<DATA>){
			if($count == 0){
				$_ =~ s/Scientific Name/Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/;
				$count = 1;
			}else{
				my @temp = split(/\t/, $_);
				$temp[1] =~ s/[kpcofgs]__//g;
				$temp[1] =~ s/\;/\t/g;
				$_ = join("\t", @temp);
			}
			print OUT "$_";
		}
		close(DATA);
		close(OUT);
	}
}
