#!/usr/bin/perl
use strict;
use warnings;

#2025/12/24 add LCA
#2025/01/29 add clustering
#2018/01/05 add coverage filter (not less than 90%)


#Setting
my ($db, $identity, $identity2, $dictionary, $denoise, $cluster, $algo);
open (SET, "<", "Setting.txt") or die("error:$!");
while(<SET>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^DB\s*=\s*([\w\W]+)/){$db = $1;}
	elsif($_ =~ /^Algorithm\s*=\s*(\d+)/){$algo = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Common_name\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^Clustering\s*=\s*(\S+)/){$cluster = $1;}
}
close(SET);
unless($db){print "Error: Please check the database name in Setting.txt.\n"; exit;}
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
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "Error: Dictionary file for common names not found or name mismatch in Setting.txt.\n"; exit;}
}
unless($denoise){$denoise = "no";}
unless($cluster){$cluster = "no";}
unless($algo){$algo = 1;}

my $as_us;
unless($denoise =~ /^yes$/i){
	$as_us = "Unique";
}else{
	if($algo == 1){
		$as_us = "ZOTU";
	}else{
		$as_us = "ASV";
	}
}

#Usearch check
opendir (DIR, ".\/Tools") or die ("error:$!");
my @tool = readdir DIR;
my ($usearch, $vv);
foreach (@tool) {
	if ($_ =~ /(usearch.*)/){$usearch = $1;last;}
	elsif($_ =~ /(vsearch.*)/){$usearch = $1; $vv = 1;}
}
closedir DIR;
unless($usearch){print "Error: Please ensure the USEARCH/VSEARCH executable file is placed in the Tools directory.\n"; exit;}

#Get sample site name
open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv") or die("error:$!");
my $tsv = 0;
my @file;
while(<DATA>){
	if($tsv == 0){
		chomp($_);
		$_ =~ s/\r//g;
		@file = split(/\t/, $_);
		shift @file;
		$tsv++;
	}else{
		last;
	}
}
close(DATA);

#Get neighbors
my %kinrin;
my $idname;
open (DATA, "<", ".\/Results\/3_1_Usearch_global\/neighbors.txt") or die("error:$!");
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	if($_ =~ /^>(\S+)/){
		$idname = $1;
	}else{
		if($_ =~ /\|k__/){
			$_ =~ /\t([^\|]+)\|([^\t]+)\t([^\t]+)/;
			unless($kinrin{$idname}){$kinrin{$idname} = "\t$1\t$2\t$3\n";}
			else{$kinrin{$idname} .= "\t$1\t$2\t$3\n";}
		}else{
			$_ =~ /\|([^\|]+)\|([^\t]+)\t([^\t]+)/;
			unless($kinrin{$idname}){$kinrin{$idname} = "\t$1\t$2\t$3\n";}
			else{$kinrin{$idname} .= "\t$1\t$2\t$3\n";}
		}
	}
}
close(DATA);

print "============================================================\n";
print "                      4_1_Annotation                       \n";
print "============================================================\n";

#Annotation
mkdir ".\/Results\/4_1_Annotation";

#import dictionary
my (%jp, %genus);
if($dictionary){
	open (DATAFILE, "<", ".\/Dictionary\/$dictionary") or die("error:$!");
	while(<DATAFILE>){
		chomp($_);
		$_ =~ s/\r//g;
		my @temp = split(/\t/, $_);
		if($temp[1]){$jp{$temp[0]} = $temp[1];}
		if($temp[1] and $temp[0] and $temp[0] =~ /^([^ ]+) /){
			unless($genus{$1}){
				my $kensaku = $1;
				if($dictionary =~ /Sname_Jname/){$genus{$kensaku} = "${temp[1]}Ç∆ìØÇ∂ëÆÇÃéÌ";}
				else{$genus{$kensaku} = "The species belongs to the same genus as ${temp[1]}";}
			}
		}
	}
	close(DATAFILE);
}
#import fasta data and number of reads
my ($data, $dseq, %lead, %each_num, @lead);
if(-f ".\/Results\/2_7_LULU\/curated.fa"){
	open (DATA, "<", ".\/Results\/2_7_LULU\/curated.fa") or die("error:$!");
	open (READ, "<", ".\/Results\/2_7_LULU\/curated_list.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_6_OTU_Clustering\/OTU.fa"){
	open (DATA, "<", ".\/Results\/2_6_OTU_Clustering\/OTU.fa") or die("error:$!");
	open (READ, "<", ".\/Results\/2_6_OTU_Clustering\/OTU_table.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa"){
	open (DATA, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_seq_rf.fa") or die("error:$!");
	open (READ, "<", ".\/Results\/2_5_Rarefaction\/${as_us}s_table_rf.tsv") or die("error:$!");
}elsif(-f ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa"){
	open (DATA, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_seq.fa") or die("error:$!");
	open (READ, "<", ".\/Results\/2_4_${as_us}s\/${as_us}s_table.tsv") or die("error:$!");
}
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	unless($_ =~ /^\n/){
		if($_ =~ /^>(.+)\;size/){
			if($data){$lead{$data} = $dseq; $data = $1; undef($dseq);}
			else{$data = $1;}
			push(@lead, $1);
		}else{
			if($dseq){$dseq = $dseq . $_;}
			else{$dseq = $_;}
		}
	}
}
if($data){$lead{$data} = $dseq;}
close(DATA);

while(<READ>){
	chomp($_);
	$_ =~ s/\r//g;
	unless($_ =~ /^\t/){
		my @temp = split(/\t/, $_);
		for(my $n = 0; $n < @file; $n++){
			$each_num{$temp[0]}{$file[$n]} = $temp[$n+1];
		}
	}
}
close(READ);
	

#usearch_results.txt
open (DATA, "<", ".\/Results\/3_1_Usearch_global\/usearch_results.txt") or die("error:$!");
my (%twohit, @blastdata, %num);
while(<DATA>){
	chomp($_);
	$_ =~ s/\r//g;
	$_ =~ s/\t+/\t/g;
	my @temp = split(/\t/, $_);
	if($temp[1] =~ /^gb\||^id\|/){$temp[1] =~ s/.+\|//;}
	if($temp[-1] < 90){next;} #query cover %
	unless($twohit{$temp[0]}){$twohit{$temp[0]} = $temp[1]; push(@blastdata, $_); $num{$temp[0]}++;}
	else{
		if($num{$temp[0]} == 1 and $twohit{$temp[0]} ne $temp[1]){push(@blastdata, $_); $num{$temp[0]}++;}
	}
}
close(DATA);

#Taxonomy_assignment_table
open (DATA, "<", ".\/Results\/3_1_Usearch_global\/Taxonomy_assignment_table.tsv") or die("error:$!");
my (%conserv);
while(<DATA>){
	unless($_ =~ /^ID/){
		my @temp = split(/\t/, $_);
		$conserv{$temp[0]} = $temp[3];
	}
}
close(DATA);

my (%all_kekka, %all_kekka2, %check, %count, %up, %down, $name, %lowhit);
foreach(@blastdata){
	my @temp = split(/\t/, $_);
	unless($temp[0] =~ /\;$/){$temp[0] = $temp[0] . "\;";}
	$check{$temp[0]}++;
	my ($gbn, $ketu);
	if($temp[1] =~ /^gb\||^id\|/){
		$temp[1] =~ /\|(.+)\|/;
		$gbn = $1;
		$temp[1] =~ s/.+\|//;
	}elsif($temp[1] =~ /^([^\|]+)\|/){
		$gbn = $1;
		$temp[1] =~ s/[^\|]+\|//;
	}else{
		$gbn = "na";
	}
	if($check{$temp[0]} == 1 and $temp[2] >= $identity){
		$name = $temp[1];
		$ketu = "$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
		$all_kekka{$temp[1]}{$temp[0]} = $ketu;
		$up{$temp[0]}++;
		$count{$temp[1]}++;
	}elsif($check{$temp[0]} == 1 and $temp[2] < $identity){
		$name = $temp[1];
		$ketu = "$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
		$all_kekka2{$temp[1]}{$temp[0]} = $ketu;
		my $tempname = $temp[0];
		$tempname =~ /(.+)\;size/;
		$lowhit{$temp[1]}{$temp[0]} = $lead{$1};
		$down{$temp[0]}++;
	}elsif($check{$temp[0]} == 2){
		$temp[1] =~ s/_/ /g;
		$ketu = "$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$gbn";
		if($up{$temp[0]}){
			$all_kekka{$name}{$temp[0]} = $all_kekka{$name}{$temp[0]} . "\t$ketu";
		}elsif($down{$temp[0]}){
			$all_kekka2{$name}{$temp[0]} = $all_kekka2{$name}{$temp[0]} . "\t$ketu";
		}
	}
}

#clusterization under $identity seqs
my %all_seiri;
foreach (sort {$a cmp $b} keys %lowhit){
	my $spname = $_;
	my $title = $lowhit{$spname};
	my $seqnum = keys %$title;
	if($seqnum > 1){
		my @keys = keys %$title;
		my (%narabi, %tempname, %modosi);
		my $tempcount = 0;
		foreach(@keys){$_ =~ /size=(\d+)/; $narabi{$1}{$_}++; $tempcount++; $tempname{$_} = "Uniq$tempcount;size=$1;"; $modosi{"Uniq$tempcount;size=$1;"} = $_;}
		undef(@keys);
		open (TEMP, ">", ".\/Results\/4_1_Annotation\/temp.fas") or die("error:$!");
		foreach(sort {$b <=> $a} keys %narabi){
			my $suji = $narabi{$_};
			foreach(sort keys %$suji){print TEMP ">$tempname{$_}\n$lowhit{$spname}{$_}\n";}
		}
		close(TEMP);
		
		my $tempid = $identity/100;
		my $command;
		if($vv){
			$command = ".\/Tools\/$usearch --cluster_smallmem \".\/Results\/4_1_Annotation\/temp.fas\" --id $tempid --uc \".\/Results\/4_1_Annotation\/temp.uc\" --usersort";
		}else{
			$command = ".\/Tools\/$usearch -cluster_smallmem \".\/Results\/4_1_Annotation\/temp.fas\" -id $tempid -uc \".\/Results\/4_1_Annotation\/temp.uc\" -sortedby size -quiet";
		}
		system $command;
		
		open (TEMP, "<", ".\/Results\/4_1_Annotation\/temp.uc") or die("error:$!");
		my $sp = 0;
		my %otus;
		my $s_count = 0;
		while(<TEMP>){
			if($_ =~ /^S\t/){$s_count++;}
		}
		close(TEMP);
		open (TEMP, "<", ".\/Results\/4_1_Annotation\/temp.uc") or die("error:$!");
		while(<TEMP>){
			chomp($_);
			$_ =~ s/\r//g;
			my @temp = split(/\t/, $_);
			if($s_count == 1){
				if($temp[0] eq "S"){$sp++; $all_seiri{$spname}{$modosi{$temp[8]}}++; $otus{$temp[8]} = $spname;}
				elsif($temp[0] eq "H"){$all_seiri{$otus{$temp[9]}}{$modosi{$temp[8]}}++;}
			}else{
				if($temp[0] eq "S"){$sp++; $all_seiri{"${spname}_otu$sp"}{$modosi{$temp[8]}}++; $otus{$temp[8]} = "${spname}_otu$sp";}
				elsif($temp[0] eq "H"){$all_seiri{$otus{$temp[9]}}{$modosi{$temp[8]}}++;}
			}
		}
		close(TEMP);
		unlink ".\/Results\/4_1_Annotation\/temp.fas";
		unlink ".\/Results\/4_1_Annotation\/temp.uc";
	}else{
		my @keys = keys %$title;
		$all_seiri{$spname}{$keys[0]}++;
	}
}

my $body_color;
if($denoise =~ /^yes$/i){
	if($algo ==1){
		$body_color = "EFEFFB";
	}else{
		$body_color = "EFFBEF";
	}
}else{
	unless($cluster =~ /^yes$/i){$body_color = "EEEEEE";}
	else{$body_color = "FBFBEF";}
}

#annotate
foreach(@file){
	print "$_ Annotate...\n";
	my $file = $_;
	
	my (%count, %count2, %kekka, %kekka2, %seiri, %nohit_check);
	foreach(keys %all_kekka){
		my $temp_sp = $_;
		foreach(keys %{$all_kekka{$_}}){
			my $temp_name = $_;
			my $temp_name2 = $_;
			$temp_name =~ s/\;size.+//;
			if($each_num{$temp_name}{$file}){
				unless($count{$temp_sp}){$count{$temp_sp} = $each_num{$temp_name}{$file};}
				else{$count{$temp_sp} += $each_num{$temp_name}{$file};}
				$_ =~ s/\;size=\d+/\;size=$each_num{$temp_name}{$file}/;
				$kekka{$temp_sp}{$_} = $all_kekka{$temp_sp}{$temp_name2};
				$nohit_check{$temp_name}++;
			}
		}
	}
	foreach(keys %all_kekka2){
		my $temp_sp = $_;
		foreach(keys %{$all_kekka2{$_}}){
			my $temp_name = $_;
			my $temp_name2 = $_;
			$temp_name =~ s/\;size.+//;
			if($each_num{$temp_name}{$file}){
				$_ =~ s/\;size=\d+/\;size=$each_num{$temp_name}{$file}/;
				$kekka2{$temp_sp}{$_} = $all_kekka2{$temp_sp}{$temp_name2};
				$nohit_check{$temp_name}++;
			}
		}
	}
	foreach(keys %all_seiri){
		my $temp_sp = $_;
		foreach(keys %{$all_seiri{$_}}){
			my $temp_name = $_;
			my $temp_name2 = $_;
			$temp_name =~ s/\;size.+//;
			if($each_num{$temp_name}{$file}){
				unless($count2{$temp_sp}){$count2{$temp_sp} = $each_num{$temp_name}{$file};}
				else{$count2{$temp_sp} += $each_num{$temp_name}{$file};}
				$_ =~ s/\;size=\d+/\;size=$each_num{$temp_name}{$file}/;
				$seiri{$temp_sp}{$_} = $all_seiri{$temp_sp}{$temp_name2};
			}
		}
	}
	
	my (@nohit, %nohit_read);
	my $total_nohit = 0;
	foreach(@lead){
		if($each_num{$_}{$file}){
			unless($nohit_check{$_}){
				$nohit_read{$_} = $each_num{$_}{$file};
				$total_nohit += $each_num{$_}{$file};
				push(@nohit, $_);
			}
		}
	}
	
	my (%reverse, %reverse2);
	my $high_read = 0;
	my $low_read = 0;
	foreach(keys %count){$reverse{$count{$_}}{$_}++; $high_read += $count{$_};}
	foreach(keys %count2){$reverse2{$count2{$_}}{$_}++; $low_read += $count2{$_};}
	my $species = 0;
	my $species2 = 0;
	$species = keys %count;
	$species2 = keys %count2;
	
	my $total_read = $high_read + $low_read + $total_nohit;
	my ($high_p, $low_p, $nohit_p); 
	if($high_read){$high_p = sprintf("%.1f", ($high_read/$total_read)*100);}
	if($low_read){$low_p = sprintf("%.1f", ($low_read/$total_read)*100);}
	if($total_nohit){$nohit_p = sprintf("%.1f", ($total_nohit/$total_read)*100);}

	my %synonym;
	open (OUT1, ">", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	open (OUT2, ">", ".\/Results\/4_1_Annotation\/${file}_Detail.txt") or die("error:$!");
	open (OUT3, ">", ".\/Results\/4_1_Annotation\/${file}_Representative_seq.fas") or die("error:$!");
	open (OUT4, ">", ".\/Results\/4_1_Annotation\/${file}_neighbors.txt") or die("error:$!");
	open (OUT5, ">", ".\/Results\/4_1_Annotation\/${file}_all_annotated_seq.fas") or die("error:$!");
	print OUT1 "===Results of $file===\n\n";
	if(%reverse){
		if($dictionary){print OUT1 "Identity >= ${identity}%\t$species species\; Total $high_read\/$total_read reads (${high_p}\%)\nScientific_name\tCommon_name\tReads\tIdentity(%)\tConfidence\tLCA Rank\n";}
		else{print OUT1 "Identity >= ${identity}%\t$species species\; Total $high_read\/$total_read reads (${high_p}\%)\nScientific_name\tReads\tIdentity(%)\tConfidence\tLCA Rank\n";}
		print OUT2 "Identity >= ${identity}%\t$species species\; Total $high_read\/$total_read reads (${high_p}\%)\n";
		foreach (sort {$b <=> $a} keys %reverse){
			my $temp = $_;						#Total number of reads
			my $kame = $reverse{$temp};
			foreach(sort {$a cmp $b} keys %$kame){
				my $temp2 = $_;					#scientific name
				my $cut_cluster = $_;
				$cut_cluster =~ s/_cluster\d+//;
				$cut_cluster =~ s/_/ /g;
				$_ =~ s/_/ /g;
				my $wamei;
				#translate
				if($dictionary){
					if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
					else{
						$_ =~ /^([^ ]+) /;
						if($genus{$1}){$wamei = $genus{$1};}
						else{
							if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
							else{$wamei = "No applicable name";}
						}
					}
					print OUT1 "$_\t$wamei\t$temp\t";
				}else{print OUT1 "$_\t$temp\t";}
				
				print OUT2 "$_\t$temp reads\n";
				my $kame2 = $kekka{$temp2};
				my %temp;
				foreach(keys %$kame2){
					$_ =~ /size=(\d+)/;
					$temp{$1} = $_;
				}
				print OUT2 "\tID\tReads\tIdentity(%)\tLOD_score\tConfidence\tLCA Rank\tAlign_len\tMismatch\tAccession No.\t2nd-sp_name\t2nd_Identity(%)\t2nd_Align_len\t2nd_Mismatch\tAccession No.\tSequence\n";
				my $fasta = 0;
				foreach(sort {$b <=> $a} keys %temp){
					$fasta++;
					my($id, $cut);
					my $read = $temp{$_};
					if($read =~ /^OTU/){
						$read =~ /^(OTU\d+)\;size=(\d+)\;/;
						$id = $1;
						$cut = $2;
					}else{
						$read =~ /^(${as_us}\d+)\;size=(\d+)\;/;
						$id = $1;
						$cut = $2;
					}
					my $kazu = 0;
					my $kiri = $kekka{$temp2}{$read};
					while($kiri =~ s/\t//i){$kazu++;}
					my $ori_read = $read;
					$read =~ /(.+)\;size/;
					$read = $1;
					
					#Make Representative, all_annotated_seq
					print OUT5 ">${temp2}_${id}_${cut}_reads\n$lead{$read}\n";
					if($fasta == 1){
						print OUT3 ">${temp2}_${id}_${temp}_reads\n$lead{$read}\n";
					}
					
					#Make neighbors
					if($fasta == 1 and $kinrin{$id}){
						my @bara = split(/\n/, $kinrin{$id});
						print OUT4 ">$temp2\n";
						foreach(@bara){
							$_ =~ /\t[^\t]+\t([^\t]+)/;
							my $koumoku = $_;
							my $spname = $1;
							my $ori = $spname;
							my $cut_cluster = $spname;
							$cut_cluster =~ s/_cluster\d+//;
							$cut_cluster =~ s/_/ /g;
							$spname =~ s/_/ /g;
							my $wamei;
							#translate
							if($dictionary){
								if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
								else{
									$spname =~ /^([^ ]+) /;
									if($genus{$1}){$wamei = $genus{$1};}
									else{
										if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
										else{$wamei = "No applicable name";}
									}
								}
								$koumoku =~ s/$ori/$ori\t$wamei/;
								print OUT4 "$koumoku\n";
							}else{
								print OUT4 "$koumoku\n";
							}
						}
					}
					
					#Make Detail file
					if($kazu > 6){
						my @temp = split(/\t/, $kekka{$temp2}{$ori_read});
						my $lod = log((($temp[1])/($temp[2]+1))/(($temp[6])/($temp[7]+1)));
						$lod = sprintf("%.4f", $lod);
						my $conf;
						if($lod >= 0.9){$conf = "HIGH";}
						elsif($lod >= 0.5){$conf = "MODERATE";}
						else{$conf = "LOW";}
						if($fasta == 1){print OUT1 "$temp[0]\t$conf\t$conserv{$id}\n";}
						print OUT2 "\t$id\t$cut\t$temp[0]\t$lod\t$conf\t$conserv{$id}\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$lead{$read}\n";}
					else{
						my @temp = split(/\t/, $kekka{$temp2}{$ori_read});
						my $lod = "N.A.";
						if($fasta == 1){print OUT1 "$temp[0]\tHIGH\t$conserv{$id}\n";}
						print OUT2 "\t$id\t$cut\t$temp[0]\t$lod\tHIGH\t$conserv{$id}\t$temp[1]\t$temp[2]\t$temp[3]\tN.A.\tN.A.\tN.A.\tN.A.\tN.A.\t$lead{$read}\n";
					}
				}
		 	}
		}
	}else{
		print OUT1 "\nIdentity >= ${identity}%\t0 species\nNo applicable sequence\n";
		print OUT2 "\nIdentity >= ${identity}%\t0 species\nNo applicable sequence\n";
	}

	#lower identity results
	if(%reverse2){
		if($dictionary){print OUT1 "\n${identity}% > Identity > ${identity2}%\t$species2 species\; Total $low_read\/$total_read reads (${low_p}\%)\nScientific_name\tCommon_name\tReads\tIdentity(%)\tConfidence\tLCA Rank\n";}
		else{print OUT1 "\n${identity}% > Identity > ${identity2}%\t$species2 species\; Total $low_read\/$total_read reads (${low_p}\%)\nScientific_name\tReads\tReads\tIdentity(%)\tConfidence\tLCA Rank\n";}
		print OUT2 "\n${identity}% > Identity > ${identity2}%\t$species2 species\; Total $low_read\/$total_read reads (${low_p}\%)\n";
		foreach (sort {$b <=> $a} keys %reverse2){
			my $temp = $_;						#Total number of reads
			my $kame = $reverse2{$temp};
			foreach(sort {$a cmp $b} keys %$kame){
				my $temp2 = $_;					#scientific name
				my $cut_cluster = $_;
				$cut_cluster =~ s/_cluster\d+//;
				$cut_cluster =~ s/_/ /g;
				$_ =~ s/_/ /g;
				my $wamei;
				#translate
				if($dictionary){
					if($jp{$cut_cluster}){
						if($dictionary =~ /Sname_Jname/){$wamei = "$jp{$cut_cluster}Ç∆ãﬂâèÇ»éÌ";}
						else{$wamei = "The species closely rerated to $jp{$cut_cluster}";}
					}else{
						$_ =~ /^([^ ]+) /;
						if($genus{$1}){
							my $itiji = $1;
							if($dictionary =~ /Sname_Jname/){$wamei = "$genus{$itiji}Ç∆ãﬂâèÇ»éÌ";}
							else{$wamei = "The species closely rerated to $genus{$itiji}";}
						}else{
							if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
							else{$wamei = "No applicable name";}
						}
					}
					print OUT1 "U${identity} $_\t$wamei\t$temp\t";
				}else{print OUT1 "U${identity} $_\t$temp\t";}
				
				print OUT2 "U${identity}_$_\t$temp reads\n";
				my $kame2 = $seiri{$temp2};
				my %temp;
				foreach(keys %$kame2){
					$_ =~ /size=(\d+)/;
					$temp{$1} = $_;
				}
				print OUT2 "\tID\tReads\tIdentity(%)\tLOD_score\tConfidence\tLCA Rank\tAlign_len\tMismatch\tAccession No.\t2nd-sp_name\t2nd_Identity(%)\t2nd_Align_len\t2nd_Mismatch\tAccession No.\tSequence\n";
				my $fasta = 0;
				foreach(sort {$b <=> $a} keys %temp){
					$fasta++;
					my $read = $temp{$_};
					my ($id, $cut);
					if($read =~ /^OTU/){
						$read =~ /^(OTU\d+)\;size=(\d+)\;/;
						$id = $1;
						$cut = $2;
					}else{
						$read =~ /^(${as_us}\d+)\;size=(\d+)\;/;
						$id = $1;
						$cut = $2;
					}
					my $kazu = 0;
					my $temp3 = $temp2;
					$temp3 =~ s/_otu\d+//;
					my $kiri = $kekka2{$temp3}{$read};
					while($kiri =~ s/\t//i){$kazu++;}
					my $ori_read = $read;
					$read =~ /(.+)\;size/;
					$read = $1;
					
					#Make Representative, all_annotated_seq
					print OUT5 ">U${identity}_${temp2}_${id}_${cut}_reads\n$lead{$read}\n";
					if($fasta == 1){print OUT3 ">U${identity}_${temp2}_${id}_${temp}_reads\n$lead{$read}\n";}
					
					#Make neighbors
					if($fasta == 1 and $kinrin{$id}){
						my @bara = split(/\n/, $kinrin{$id});
						print OUT4 ">U${identity}_$temp2\n";
						foreach(@bara){
							$_ =~ /\t[^\t]+\t([^\t]+)/;
							my $koumoku = $_;
							my $spname = $1;
							my $ori = $spname;
							my $cut_cluster = $spname;
							$cut_cluster =~ s/_cluster\d+//;
							$cut_cluster =~ s/_/ /g;
							$spname =~ s/_/ /g;
							my $wamei;
							#translate
							if($dictionary){
								if($jp{$cut_cluster}){$wamei = $jp{$cut_cluster};}
								else{
									$spname =~ /^([^ ]+) /;
									if($genus{$1}){$wamei = $genus{$1};}
									else{
										if($dictionary =~ /Sname_Jname/){$wamei = "äYìñëÆñºÅEòañºÇ»Çµ";}
										else{$wamei = "No applicable name";}
									}
								}
								$koumoku =~ s/$ori/$ori\t$wamei/;
								print OUT4 "$koumoku\n";
							}else{
								print OUT4 "$koumoku\n";
							}
						}
					}
					#Make Detail file
					if($kazu > 6){
						my @temp = split(/\t/, $kekka2{$temp3}{$ori_read});
						my $lod = log((($temp[1])/($temp[2]+1))/(($temp[6])/($temp[7]+1)));
						$lod = sprintf("%.4f", $lod);
						my $conf;
						if($lod >= 0.9){$conf = "HIGH";}
						elsif($lod >= 0.5){$conf = "MODERATE";}
						else{$conf = "LOW";}
						if($fasta == 1){print OUT1 "$temp[0]\t$conf\t$conserv{$id}\n";}
						print OUT2 "\t$id\t$cut\t$temp[0]\t$lod\t$conf\t$conserv{$id}\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$lead{$read}\n";}
					else{
						my @temp = split(/\t/, $kekka2{$temp3}{$ori_read});
						my $lod = "N.A.";
						if($fasta == 1){print OUT1 "$temp[0]\tHIGH\t$conserv{$id}\n";}
						print OUT2 "\t$id\t$cut\t$temp[0]\t$lod\tHIGH\t$conserv{$id}\t$temp[1]\t$temp[2]\t$temp[3]\tN.A.\tN.A.\tN.A.\tN.A.\tN.A.\t$lead{$read}\n";
					}
				}
		 	}
		}
	}else{
		print OUT1 "\n${identity}% > Identity > ${identity2}%\t0 species\nNo applicable sequence\n";
		print OUT2 "\n${identity}% > Identity > ${identity2}%\t0 species\nNo applicable sequence\n";
	}
	

	#Output nohit sequence
	my %read_jun;
	foreach(@nohit){
		push(@{$read_jun{$nohit_read{$_}}}, $_);
	}
	
	my $nohit_n = 0;
	unless(@nohit){$nohit_n = 0;}
	else{$nohit_n = @nohit;}
	
	if(@nohit){
		print OUT1 "\nNohit\t$nohit_n sequences\; Total $total_nohit\/$total_read reads (${nohit_p}\%)\n";
		print OUT2 "\nNohit\t$nohit_n sequences\; Total $total_nohit\/$total_read reads (${nohit_p}\%)\n";
		foreach(sort {$b <=> $a} keys %read_jun){
			my $temp = $_;
			foreach(@{$read_jun{$temp}}){
				print OUT1 ">$_\;size=$temp\n$lead{$_}\n";
				print OUT2 ">$_\;size=$temp\n$lead{$_}\n";
				print OUT3 ">Nohit_${_}_${temp}_reads\n$lead{$_}\n";
				print OUT5 ">Nohit_${_}_${temp}_reads\n$lead{$_}\n";
			}
		}
	}else{
		print OUT1 "\nNohit\t0 sequences\n";
		print OUT2 "\nNohit\t0 sequences\n";
	}
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(OUT4);
	close(OUT5);
	
	
	#Make Detail html
	open (OUT, ">", ".\/Results\/4_1_Annotation\/${file}_Detail.html") or die("error:$!");
	print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>${file}_Detail</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#$body_color">
<a name="top"></a>
<font face="Tahoma" size="6">Detailed Result of $file</font><BR>
EOS
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Detail.txt") or die("error:$!");
	my $iro = 0;
	my $table = 1;
	my $nohit = 0;
	my ($name2, %portal);
	while(<DATA>){
		if($_ =~ /^\n|^\r/){next;}
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ /^Identity/){print OUT "<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";next;}
		elsif($_=~ /^${identity}%/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$iro = 0; $table = 1; next;}
		elsif($_ =~ /^Nohit/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n"; $nohit = 1; next;}
		unless($_ =~ s/^\t//){
			if($nohit){print OUT "<font face=\"Tahoma\" size=\"3\">$_</font><BR>\n";next;}
			$_ =~ /([^\t]+)\t/;
			if($1){$name2 = $1;}
		}else{
			my @temp = split(/\t/, $_);
			if($_ =~ /^ID/){
				if($table == 1){
					print OUT "<table><table border=\"0\" cellspacing=\"1\" bgcolor=\"\#191970\" border=\"1\" align=\"left\">\n<tr>\n";
					print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name&nbsp;&nbsp;</H1></th>\n";
					foreach(@temp){print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>$_&nbsp;&nbsp;</H1></th>\n";}
					print OUT "</tr>\n";
					$table++;
				}else{next;}
			}else{$iro++; $portal{$name2}++; &nakami3($iro, \@temp, $name2, \%portal);}
		}
	}
	close(DATA);
	close(OUT);
	unlink ".\/Results\/4_1_Annotation\/${file}_Detail.txt";
	
#Make neighbor_list html
	open (OUT, ">", ".\/Results\/4_1_Annotation\/${file}_neighbor_list.html") or die("error:$!");
	print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>${file}_neighbor_list</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#$body_color">
<a name="top"></a>
<font face="Tahoma" size="6">Neighbor list of $file</font><BR>
EOS
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_neighbors.txt") or die("error:$!");
	$table = 0;
	while(<DATA>){
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ s/^>//){
			$_ =~ s/_/ /g;
			$table++;
			if($table > 1){print OUT "</table>\n<br clear = \"all\"><br>\n";}
			$iro = 0;
			print OUT "<font face=\"Tahoma\" size=\"4\"><B><I>$_</I>&nbsp;&nbsp;</B></font><BR>\n";
			print OUT "<table><table border=\"0\" cellspacing=\"1\" bgcolor=\"\#191970\" border=\"1\" align=\"left\">\n<tr>\n";
			if($dictionary){
				print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Accession No.&nbsp;&nbsp;</H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name&nbsp;&nbsp;</H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Common name&nbsp;&nbsp;</H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Identity(%)&nbsp;&nbsp;</H1></th>\n</tr>\n";
			}else{
				print OUT "<th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Accession No.&nbsp;&nbsp;</H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Scientific name&nbsp;&nbsp;</H1></th><th bgcolor=\"\#00008B\" align=\"left\" nowrap><H1>Identity(%)&nbsp;&nbsp;</H1></th>\n</tr>\n";
			}
		}else{
			$_ =~ s/^\t//;
			my @temp = split(/\t/, $_);
			$temp[1] =~ s/_/ /g;
			print OUT "<tr>\n";
			if($dictionary){
				print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[0]\" target=Åh_blankÅhrel=\"noopener noreferrer\"><H2>$temp[0]&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[1]</I>&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[3]&nbsp;&nbsp;</H2></td>\n";
			}else{
				print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[0]\" target=Åh_blankÅhrel=\"noopener noreferrer\"><H2>$temp[0]&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[1]</I>&nbsp;&nbsp;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]&nbsp;&nbsp;</H2></td>\n";
			}
			print OUT "</tr>\n";
		}
	}
	print OUT "</table>\n<br clear = \"all\"><br>\n";
	close(DATA);
	close(OUT);
	unlink ".\/Results\/4_1_Annotation\/${file}_neighbors.txt";
}

sub nakami3{
	my ($iro, $temp, $name, $portal) = @_;
	my @temp = @$temp;
	my %portal = %$portal;
	my $synonym;
	if($iro%2 + 1 == 2){
		print OUT "<tr>\n";
		if($portal{$name} and $portal{$name} == 1){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$name</I>&nbsp;&nbsp;</H2></td>\n";}
		else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2></H2></td>\n";}
		for(my $n = 0; $n < @temp; $n++){
			if($n == 8 or $n == 13){
				unless($temp[$n] =~ /N.A./){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[$n]\" target=Åh_blankÅhrel=\"noopener noreferrer\"><H2>$temp[$n]</H2></a></td>\n";}
				else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
			}elsif($n == 9){print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[$n]</I>&nbsp;&nbsp;</H2></td>\n";}
			else{print OUT "<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
		}
		print OUT "</tr>\n";
	}else{
		print OUT "<tr>\n";
		if($portal{$name} and $portal{$name} == 1){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$name</I>&nbsp;&nbsp;</H2></td>\n";}
		else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2></H2></td>\n";}
		for(my $n = 0; $n < @temp; $n++){
			if($n == 8 or $n == 13){
				unless($temp[$n] =~ /N.A./){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/$temp[$n]\" target=Åh_blankÅhrel=\"noopener noreferrer\"><H2>$temp[$n]</H2></a></td>\n";}
				else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
			}elsif($n == 9){print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[$n]</I>&nbsp;&nbsp;</H2></td>\n";}
			else{print OUT "<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[$n]&nbsp;&nbsp;</H2></td>\n";}
		}
		print OUT "</tr>\n";
	}
}
