#!/usr/bin/perl
use strict;
use warnings;

#Setting
my ($algo, $primer, $trim, $denoise, $identity, $identity2, $dictionary, $temporary, $cluster, $type, $lulu, $algo2);
open (SET, "<", "Setting.txt") or die("error:$!");
open (OUT, ">", ".\/Results\/Setting_log.txt") or die("error:$!");
my $logger = 0;
while(<SET>){
	if($logger == 0){print OUT "$_";}
	if($_ =~ /^Algorithm\s*=\s*(\d+)/){$algo = $1;}
	elsif($_ =~ /^Type\s*=\s*(\d+)/){$type = $1;}
	elsif($_ =~ /^Primers\s*=\s*(\S+)/){my $temp = $1;if($temp =~ /^no$/i){$primer = "No";}else{$primer = $temp;}}
	elsif($_ =~ /^Denoise\s*=\s*(\S+)/){$denoise = $1;}
	elsif($_ =~ /^UIdentity\s*=\s*(\S+)/){$identity = $1;}
	elsif($_ =~ /^LIdentity\s*=\s*(\S+)/){$identity2 = $1;}
	elsif($_ =~ /^Common_name\s*=\s*(\S+)/){$dictionary = $1;}
	elsif($_ =~ /^Clustering\s*=\s*(\S+)/){$cluster = $1;}
	elsif($_ =~ /^Algorithm2\s*=\s*(\S+)/){$algo2 = $1;}
	elsif($_ =~ /^Curation\s*=\s*(\S+)/){$lulu = $1;}
}
close(SET);
close(OUT);
unless($algo){$algo = 1;}
unless($algo2){$algo2 = 1;}
unless($type){$type = 1;}
unless($primer){print "Error: Please check the Primer file name in Setting.txt.\n"; exit;}
unless($denoise){$denoise = "yes";}
unless($identity){print "Error: Please check the upper threshold for homology search in Setting.txt.\n"; exit;}
unless($identity2){print "Error: Please check the lower threshold for homology search in Setting.txt.\n"; exit;}
unless(-f ".\/Dictionary\/$dictionary"){
	if($dictionary =~ /^no$/i){undef($dictionary);}
	else{print "Error: Dictionary file for common names not found or name mismatch in Setting.txt.\n"; exit;}
}
unless($cluster){$cluster = "no";}
unless($lulu){$lulu = "no";}

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

my $usv;
if($vv){$usv = "VSEARCH";}
else{$usv = "USEARCH";}

my $cm;
if($algo2 == 1){$cm = "Uclust";}
else{$cm = "SWARM";}

#Usearch options log
open (SET, "<", "Options_usearch.txt") or die("error:$!");
open (OUT, ">", ".\/Results\/Options_usearch_log.txt") or die("error:$!");
while(<SET>){
	print OUT "$_";
}
close(SET);
close(OUT);

#Get data names
opendir (DIR, ".\/Results\/4_1_Annotation") or die ("error:$!");
my @read = readdir DIR;
my %file;
foreach (@read) {
	if ($_ =~ /(.+)_Representative_seq/){$file{$1}++;}
}
closedir DIR;


print "============================================================\n";
print "                        6_1_Portal                          \n";
print "============================================================\n";

my @youbi = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat');
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
if($min < 10){$min = "0$min";}
if($hour < 10){$hour = "0$hour";}
if($mon < 10){$mon = "0$mon";}
if($mday < 10){$mday = "0$mday";}

my $start;
if(-f ".\/Results\/log.txt"){
	open (OUT, "<", ".\/Results\/log.txt") or die("error:$!");
	while(<OUT>){
		chomp($_);
		if($_ =~ /Start/){$start = $_; last;}
	}
	close(OUT);
}else{
	$start = "Start: no data";
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

#Make portal
open (OUT, ">", ".\/Results\/Portal_$year$mon$mday$hour$min.html") or die("error:$!");
print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>Portal_$year$mon$mday$hour$min</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#$body_color">
<a name="top"></a>
<font face="Tahoma" size="3">$start</font><BR>
<font face="Tahoma" size="3">End : $year/$mon/$mday ($youbi[$wday]) $hour:$min</font><BR>
EOS

if(-f ".\/Results\/log.txt"){
	print OUT "<font face=\"Tahoma\" size=\"4\"><a href=\".\/log.html\">Log (Reads after each step)</a></font><BR>";
}
my $title;
if($denoise =~ /^yes$/i){
	if($algo == 1){
		unless($cluster =~ /^yes$/i){
			if($type == 2){
				if($lulu =~ /^yes$/i){$title = "(based on ZOTUs\; Pool Denoising\: $usv\; Curated by LULU)";}
				else{$title = "(based on ZOTUs\; Pool Denoising\: $usv)";}
				#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on ZOTUs\; Pool Denoising\: $usv)</font><BR>\n";
			}else{
				if($lulu =~ /^yes$/i){$title = "(based on ZOTUs\; Individual Denoising\:$usv\; Curated by LULU)";}
				else{$title = "(based on ZOTUs\; Individual Denoising\:$usv)";}
				#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on ZOTUs\; Individual Denoising\:$usv)</font><BR>\n";
			}
		}else{
			if($type == 2){
				if($lulu =~ /^yes$/i){$title = "(based on OTUs\: clusters of ZOTUs by $cm\; Pool Denoising\: $usv\; Curated by LULU)";}
				else{$title = "(based on OTUs\: clusters of ZOTUs by $cm\; Pool Denoising\: $usv)";}
				#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on OTUs\: clusters of ZOTUs\; Pool Denoising\: $usv)</font><BR>\n";
			}else{
				if($lulu =~ /^yes$/i){$title = "(based on OTUs\: clusters of ZOTUs by $cm\; Individual Denoising\: $usv\; Curated by LULU)";}
				else{$title = "(based on OTUs\: clusters of ZOTUs by $cm\; Individual Denoising\: $usv)";}
				#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on OTUs\: clusters of ZOTUs\; Individual Denoising\: $usv)</font><BR>\n";
			}
		}
	}else{
		unless($cluster =~ /^yes$/i){
			if($lulu =~ /^yes$/i){$title = "(based on ASVs\; Denoising\: DADA2\; Curated by LULU)";}
			else{$title = "(based on ASVs\; Denoising\: DADA2)";}
			#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on ASVs\; Denoising\: DADA2)</font><BR>\n";
		}else{
			if($lulu =~ /^yes$/i){$title = "(based on OTUs\: clusters of ASVs by $cm\; Denoising\: DADA2\; Curated by LULU)";}
			else{$title = "(based on OTUs\: clusters of ASVs by $cm\; Denoising\: DADA2)";}
			#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on OTUs\: clusters of ASVs\; Denoising\: DADA2)</font><BR>\n";
		}
	}
}else{
	unless($cluster =~ /^yes$/i){
		if($lulu =~ /^yes$/i){$title = "(Direct Taxonomic Assignment\; Curated by LULU)";}
		else{$title = "(Direct Taxonomic Assignment)";}
		#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (Direct Taxonomic Assignment)</font><BR>\n";
	}else{
		if($lulu =~ /^yes$/i){$title = "(based on OTUs\: clusters of Uniques by $cm\; Curated by LULU)";}
		else{$title = "(based on OTUs\: clusters of Uniques by $cm)";}
		#print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis (based on OTUs\: clusters of Uniques)</font><BR>\n";
	}
}
print OUT "<font face=\"Tahoma\" size=\"6\">Summary of Analysis $title</font><BR>\n";
print OUT <<"EOS";
<font face="Tahoma" size="4"><a href="./5_2_Summary_Table/Summary_Table.tsv">Summary Table</a></font><BR>
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#ff0000" align="left" nowrap><H1>Samples analyzed &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Identity >= ${identity}% species #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>${identity}% > Identity > ${identity2}% species #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Nohit seq. #&nbsp;&nbsp;</H1></th>
</tr>
EOS

my $matome = 0;
foreach(sort keys %file){
	my $file = $_;
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	my $tread = 0;
	my $upper = 0; my $lower = 0; my $nohit = 0;
	my $upper_s = 0; my $lower_s = 0; my $nohit_s = 0;
	my $upper_p = 0; my $lower_p = 0; my $nohit_p = 0;
	$matome++;
	while(<DATA>){
		if($_ =~ /^Identity.+\t(\d+) species\; Total (\d+)\/(\d+) reads \((\S+)\%\)/){
			$upper_s = $1;
			$upper = $2;
			$tread = $3;
			$upper_p = $4;
		}elsif($_ =~ /> Identity > .+\t(\d+) species\; Total (\d+)\/(\d+) reads \((\S+)\%\)/){
			$lower_s = $1;
			$lower = $2;
			$tread = $3;
			$lower_p = $4;
		}elsif($_ =~ /^Nohit\t(\d+) sequences\; Total (\d+)\/(\d+) reads \((\S+)\%\)/){
			$nohit_s = $1;
			$nohit = $2;
			$tread = $3;
			$nohit_p = $4;
			last;
		}
	}
	close(DATA);
	if($tread > 0 and $upper_s > 0){$upper = "$upper_s \($upper reads\/$upper_p\%)";}
	if($tread > 0 and $lower_s > 0){$lower = "$lower_s \($lower reads\/$lower_p\%)";}
	if($tread > 0 and $nohit_s > 0){$nohit = "$nohit_s \($nohit reads\/$nohit_p\%)";}
	if($matome%2 + 1 == 2){print OUT "<tr>\n<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><a href=\#$file><H2>$file\&nbsp\;\&nbsp\;</H2></a></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$tread\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$upper</H2></td></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$lower</H2></td></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$nohit</H2></td>\n</tr>\n";}
	else{print OUT "<tr>\n<td bgcolor=\"\#fceded\" align=\"left\" nowrap><a href=\#$file><H2>$file\&nbsp\;\&nbsp\;</H2></a></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$tread\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$upper</H2></td></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$lower</H2></td></td><td bgcolor=\"\#fceded\" align=\"left\" nowrap><H2>$nohit</H2></td>\n</tr>\n";}
}
print OUT "</table>\n<br clear = \"all\">\n";
print OUT "<a href=\#top><font face=\"Tahoma\">Back to top<font></a><br>\n";
print OUT "<hr style=\"border:0;border-top:thick dotted white;background-color:red;\">\n";
print OUT "<font face=\"Tahoma\" size=\"6\">Each Result of Sample</font><BR><br>\n";


foreach(sort keys %file){
	my $file = $_;
	print OUT "<a name=\"$_\"><font face=\"Tahoma\" size=\"5\"><B>$_</B></font></a><BR>\n<font face=\"Tahoma\" size=\"4\"><a href=\".\/4_1_Annotation\/${file}_Detail.html\" target=\"_blank\" rel=\"noopener noreferrer\">Detailed Results</a></font>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;";
	print OUT "<font face=\"Tahoma\" size=\"4\"><a href=\".\/4_1_Annotation\/${file}_neighbor_list.html\" target=\"_blank\" rel=\"noopener noreferrer\">Neighbor list</a></font>\&nbsp\;\&nbsp\;\&nbsp\;\&nbsp\;";
	print OUT "<font face=\"Tahoma\" size=\"4\"><a href=\".\/4_1_Annotation\/${file}_Representative_seq.fas\">Representative Sequences</a></font><br>\n";
	open (DATA, "<", ".\/Results\/4_1_Annotation\/${file}_Summary.txt") or die("error:$!");
	my $table = 0;
	my $iro = 0;
	while(<DATA>){
		if($_ =~ /^\n/){next;}
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ /^Identity/){print OUT "<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$table = 1;next;}
		if($_ =~ /^${identity}%/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><BR>\n";$iro = 0;$table = 2;next;}
		if($table == 1){
			if($_ =~ /^No applicable sequence/){print OUT "<font face=\"Tahoma\" size=\"2\"><B>$_</B></font><BR>\n"; $table = 2;next;}
			if($_ =~ /^Scientific/){&table;next;}
			my @temp = split(/\t/, $_);
			$iro++;
			if($dictionary){&nakami1($iro, \@temp, $file);}
			else{&nakami2($iro, \@temp, $file);}
		}
		if($table == 2){
			if($_ =~ /^No applicable sequence/){print OUT "<font face=\"Tahoma\" size=\"2\"><B>$_</B></font><BR>\n"; next;}
			if($_ =~ /^Nohit/){print OUT "</table>\n<br clear = \"all\"><br>\n<font face=\"Tahoma\" size=\"3\"><B>$_</B></font><br><br>\n<a href=\#top><font face=\"Tahoma\">Back to list<font></a><br><br><br>\n\n\n\n"; last;}
			if($_ =~ /^Scientific/){&table;next;}
			my @temp = split(/\t/, $_);
			$iro++;
			if($dictionary){&nakami1($iro, \@temp, $file);}
			else{&nakami2($iro, \@temp, $file);}
		}
	}
	close(DATA);
}
close(OUT);
print "Portal_$year$mon$mday$hour$min.html was created.\n";

#Make log.html
unless(-f ".\/Results\/log.txt"){exit;}
open (OUT, ">", ".\/Results\/log.html") or die("error:$!");
print OUT <<"EOS";
<html>
<head>
<meta http-equiv="Content-type" content="text/html" charset="Shift_JIS">
<title>Log</title>
<style type="text/css">
H1{color: #ffffff; text-align: left; font: 100% Tahoma;}
H2{color: #000000; text-align: left; font: 100% Tahoma;}
H3{text-align: left; font: 100% Tahoma; line-height: 4px;}
</style>
</head>
<body bgcolor="#$body_color">
<a name="top"></a>
<font face="Tahoma" size="6">Reads after each step</font><BR>
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
EOS

if(-f ".\/Results\/read_tracking.csv"){
	print OUT <<"EOS";
<th bgcolor="#ff0000" align="left" nowrap><H1>Samples analyzed &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Raw reads #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Filtered &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>DenoisedF &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>DenoisedR &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Merged &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>NonChimeric &nbsp;&nbsp;</H1></th>
EOS
	open (LOG, "<", ".\/Results\/read_tracking.csv") or die("error:$!");
	my $count = 0;
	while(<LOG>){
		$_ =~ s/\r//g;
		$_ =~ s/\"//g;
		chomp($_);
		if($_ =~ /Sample/){next;}
		my @log = split(/,/, $_);
		unless($count){
			print OUT "</tr>\n";
			foreach(@log){print OUT "<td bgcolor=\"#ffffff\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
			print OUT "<tr>\n";
			$count = 1;
		}else{
			print OUT "</tr>\n";
			foreach(@log){print OUT "<td bgcolor=\"#fceded\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
			print OUT "<tr>\n";
			$count = 0;
		}
	}
	close(LOG);
}else{
	print OUT <<"EOS";
<th bgcolor="#ff0000" align="left" nowrap><H1>Samples analyzed &nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Raw reads #&nbsp;&nbsp;</H1></th>
<th bgcolor="#ff0000" align="left" nowrap><H1>Merged &nbsp;&nbsp;</H1></th>
EOS
	if($trim){$primer = "no";}
	unless($primer =~ /^no$/i){
		print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Strip primer &nbsp;&nbsp;</H1></th>\n";
	}
	print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Quality filter &nbsp;&nbsp;</H1></th>\n";
	if($denoise =~ /^yes$/i){print OUT "<th bgcolor=\"#ff0000\" align=\"left\" nowrap><H1>Denoise &nbsp;&nbsp;</H1></th>\n</tr>\n";}
	open (LOG, "<", ".\/Results\/log.txt") or die("error:$!");
	my $count = 0;
	while(<LOG>){
		chomp($_);
		$_ =~ s/\r//g;
		if($_ =~ /Start/){next;}
		my @log = split(/\t/, $_);
		unless($count){
			print OUT "</tr>\n";
			foreach(@log){print OUT "<td bgcolor=\"#ffffff\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
			print OUT "<tr>\n";
			$count = 1;
		}else{
			print OUT "</tr>\n";
			foreach(@log){print OUT "<td bgcolor=\"#fceded\" align=\"left\" nowrap><H2>$_&nbsp;&nbsp;</H2></td>"}
			print OUT "<tr>\n";
			$count = 0;
		}
	}
	close(LOG);
}

print OUT "</table>\n<br clear = \"all\">\n";
unlink ".\/Results\/log.txt";

#sub
sub table{
	if($dictionary){
		print OUT <<'EOS';
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#00008B" align="left" nowrap><H1>Scientific name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Common name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Identity(%) &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Confidence &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>LCA Rank &nbsp;&nbsp;</H1></th>
</tr>
EOS
	}else{
		print OUT <<'EOS';
<table><table border="0" cellspacing="1" bgcolor="#191970" border="1" align="left">
<tr>
<th bgcolor="#00008B" align="left" nowrap><H1>Scientific name &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Total read #&nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Identity(%) &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>Confidence &nbsp;&nbsp;</H1></th>
<th bgcolor="#00008B" align="left" nowrap><H1>LCA Rank &nbsp;&nbsp;</H1></th>
</tr>
EOS
	}
}

sub nakami1{
	my ($iro, $temp, $file) = @_;
	my @temp = @$temp;
	my $synonym;
	if($iro%2 + 1 == 2){
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[4]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[5]</H2></td>
</tr>
EOS
	}else{
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[4]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[5]</H2></td>
</tr>
EOS
	}
}

sub nakami2{
	my ($iro, $temp, $file) = @_;
	my @temp = @$temp;
	my $synonym;
	if($iro%2 + 1 == 2){
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#ffffff\" align=\"left\" nowrap><H2>$temp[4]</H2></td>
</tr>
EOS
	}else{
		print OUT <<"EOS";
<tr>
<td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2><I>$temp[0]\&nbsp\;\&nbsp\;</I></H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[1]\&nbsp\;\&nbsp\;</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[2]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[3]</H2></td><td bgcolor=\"\#E6E6FA\" align=\"left\" nowrap><H2>$temp[4]</H2></td>
</tr>
EOS
	}
}
