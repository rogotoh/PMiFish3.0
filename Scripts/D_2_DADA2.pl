#!/usr/bin/perl
use strict;
use warnings;

#DADA2 options
my ($trunclenf, $trunclenr, $maxn, $maxeef, $maxeer, $truncq, $minoverlap, $maxmismatch, $quality);
open (SET, "<", "Options_DADA2.txt") or die("error:$!");
open (OUT, ">", ".\/Results\/Options_DADA2_log.txt") or die("error:$!");
while(<SET>){
	if($_ =~ /^\t\*truncLen\s*=\s*c\((\d+),(\d+)\)/){$trunclenf = $1; $trunclenr = $2;}
	elsif($_ =~ /^\t\*maxN\s*=\s*(\d+)/){$maxn = $1;}
	elsif($_ =~ /^\t\*maxEE\s*=\s*c\((\S+),(\S+)\)/){$maxeef = $1; $maxeer = $2;}
	elsif($_ =~ /^\t\*truncQ\s*=\s*(\d+)/){$truncq = $1;}
	elsif($_ =~ /^\t\*minOverlap\s*=\s*(\d+)/){$minoverlap = $1;}
	elsif($_ =~ /^\t\*maxMismatch\s*=\s*(\d+)/){$maxmismatch = $1;}
	elsif($_ =~ /^\tQuality_profile\s*=\s*(\S+)/){$quality = $1;}
	print OUT "$_";
}
close(SET);
close(OUT);
unless($trunclenf){$trunclenf = 0;}
unless($trunclenr){$trunclenr = 0;}
unless($maxn){$maxn = 0;}
unless($maxeef){$maxeef = 2;}
unless($maxeer){$maxeer = 2;}
unless($truncq){$truncq = 2;}
unless($minoverlap){$minoverlap = 12;}
unless($maxmismatch){$maxmismatch = 0;}
unless($quality){$quality = "yes";}

unless(-d ".\/Results"){mkdir ".\/Results";}
#DADA2.R
open (OUT, ">", ".\/Results\/DADA2.R") or die("error:$!");

print "============================================================\n";
print "                         2_DADA2                            \n";
print "============================================================\n";

print OUT <<"EOS";
library(dada2)

\# Specify the data directory
path <- "./Results/1_1_Primer_Trimmed_fastq"
path2 <- "./Results"

\# Get only FASTQ files (excluding _log.txt and others)
all_files <- list.files(path, pattern="\\\\.fq\$", full.names = TRUE)
fnFs <- sort(all_files[grepl("_R1_", all_files)])
fnRs <- sort(all_files[grepl("_R2_", all_files)])

\# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_stripped"), `[`, 1)

\# Set output directory
filt_path <- file.path(path2, "1_2_Filtered")
dir.create(filt_path, showWarnings = FALSE)

EOS

if($quality =~ /yes/i){
	print OUT <<"EOS";
\# Output quality profile
dir.create("./Results/Quality_profile", showWarnings = FALSE)
for (file in fnFs) {
  \# Create output file name
  output_file <- paste0("./Results/Quality_profile/Quality_profile_", basename(file), ".png")
  
  \# Start PNG output
  png(output_file, width=800, height=600)

  \# Plot quality profile
  print(plotQualityProfile(file))

  \# End PNG output
  dev.off()
}

for (file in fnRs) {
  output_file <- paste0("./Results/Quality_profile/Quality_profile_", basename(file), ".png")
  png(output_file, width=800, height=600)
  print(plotQualityProfile(file))
  dev.off()
}
EOS
}

print OUT <<"EOS";
\# Filtering and trimming
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

\# Execute filtering
filter_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c($trunclenf,$trunclenr),
              maxN=$maxn, maxEE=c($maxeef,$maxeer), truncQ=$truncq, rm.phix=TRUE)

\# Exclude samples with zero reads after filtering
valid_samples <- filter_out[,1] > 0 
filtFs <- filtFs[valid_samples]
filtRs <- filtRs[valid_samples]
sample.names <- sample.names[valid_samples]

\# Learn error models
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

\# Plot error rates for forward reads
png("./Results/error_plot_forward.png", width=800, height=600)
plotErrors(errF, nominalQ=TRUE)
dev.off()

\# Plot error rates for reverse reads
png("./Results/error_plot_reverse.png", width=800, height=600)
plotErrors(errR, nominalQ=TRUE)
dev.off()

\# Denoising (excluding samples with zero reads)
dadaFs <- lapply(filtFs, function(f) if (file.exists(f)) dada(f, err=errF, multithread=FALSE) else NULL)
dadaRs <- lapply(filtRs, function(f) if (file.exists(f)) dada(f, err=errR, multithread=FALSE) else NULL)

\# Merging (skip NULL)
mergers <- mapply(function(dF, fF, dR, fR) {
  if (!is.null(dF) & !is.null(dR)) mergePairs(dF, fF, dR, fR, minOverlap = $minoverlap, maxMismatch = $maxmismatch, verbose=TRUE) else NULL
}, dadaFs, filtFs, dadaRs, filtRs, SIMPLIFY=FALSE)

\# Create sequence table
seqtab <- makeSequenceTable(mergers)

\# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

\# Check removal rate
sum(seqtab.nochim)/sum(seqtab)

\# Set sample names as row names
rownames(seqtab.nochim) <- sample.names

\# Output as CSV file
write.csv(seqtab.nochim, "./Results/ASV_table_original.csv", row.names=TRUE)

\# Record the number of reads at each filtering and processing stage
reads_filtered <- filter_out[valid_samples, 1]  
reads_denoisedF <- sapply(dadaFs, function(x) ifelse(is.null(x), 0, sum(x\$denoised)))
reads_denoisedR <- sapply(dadaRs, function(x) ifelse(is.null(x), 0, sum(x\$denoised)))
reads_merged <- sapply(mergers, function(x) ifelse(is.null(x), 0, sum(x\$abundance)))
reads_nochim <- rowSums(seqtab.nochim)

\# Convert statistics into a data frame
tracking <- data.frame(
  Sample = sample.names,
  Input = filter_out[valid_samples, 1],
  Filtered = reads_filtered,
  DenoisedF = reads_denoisedF,
  DenoisedR = reads_denoisedR,
  Merged = reads_merged,
  NonChimeric = reads_nochim
)

\# Output as CSV file
write.csv(tracking, file="./Results/read_tracking.csv", row.names=FALSE)

EOS
close(OUT);

my $command = "Rscript .\/Results\/DADA2.R";
system $command;
