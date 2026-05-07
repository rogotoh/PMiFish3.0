# PMiFish (A Pipeline for MiFish Data Analysis) ver. 3.0

PMiFish is a streamlined and user-friendly pipeline designed for the efficient analysis of environmental DNA (eDNA) metabarcoding data. It takes FASTQ files generated from Illumina next-generation sequencers as input and performs a comprehensive workflow through to taxonomic assignment.

## 🌟 Key Features
- **Cross-platform Support**: Runs on Windows, macOS, and Linux.
- **Versatile Tool Integration**: Uses **USEARCH (Edgar, 2010)** or **VSEARCH (Rognes et al., 2016)** for core analysis, with an option to use **DADA2 (Callahan et al., 2016)** for denoising.
- **High Flexibility**:
  - Supports samples containing reads amplified by multiple primer sets, generating individual results for each primer.
  - Applicable to any taxonomic group (e.g., fish, crustaceans, fungi), not just fish.
  - Includes advanced options such as post-clustering curation with **LULU (Frøslev et al., 2017)** and clustering with **Swarm (Mahé et al., 2022)**.
- **Comprehensive Reporting**: Generates a browser-viewable **Portal (HTML)** summarizing read counts, detected species lists, and detailed taxonomic assignments.

## 📥 Getting Started

### Download the Pipeline
To use PMiFish, first download the entire script repository to your local machine using one of the following methods:

#### Method A: Manual Download (ZIP)
1. Click the green "<> Code" button at the top of this GitHub page.

2. Select "Download ZIP".

3. Extract the downloaded ZIP file to your desired location.

#### Method B: Using Git 
Open your terminal or PowerShell and run:
```bash
git clone https://github.com/rogotoh/PMiFish3.0.git
cd PMiFish
```

## 📂 Directory Structure
```text
PMiFish/
├── DataBase/    # Reference databases (FASTA) and primer files
├── Dictionary/  # Species/common name dictionaries
├── Results/     # All output files
├── Run/         # Input FASTQ or FASTQ.gz files
├── Scripts/     # Perl scripts for each processing step
├── Tools/       # USEARCH / VSEARCH / Swarm executables
├── PMiFish.pl   # Main pipeline script
├── Demultiplex.pl   # Script for demultiplex
├── PA_with_DB.pl    # Script for phylogenetic analysis
├── Setting.txt  # User configuration file
└── Options_usearch.txt / Options_DADA2.txt # Detailed tool options
```

## 🛠 Prerequisites & Installation

### 【Windows Only】
1. **Install Strawberry Perl**: [https://strawberryperl.com/](https://strawberryperl.com/).
2. **Install gzip**: [http://gnuwin32.sourceforge.net/packages/gzip.htm](http://gnuwin32.sourceforge.net/packages/gzip.htm).
   - **Important**: You must manually add the directory containing `gzip.exe` to your system **PATH**.

### 【Common Settings (Windows, macOS, Linux)】
1. **USEARCH or VSEARCH**: Place the executable in the `Tools/` directory.
2. **R & Libraries**: If using **DADA2** or **LULU**, install R and the respective libraries.
3. **MEGACC (MEGA12) (Kumer et al., 2024)**: Required for phylogenetic analysis.


## 🚀 Quick Start

1. **Configure Settings**: Edit `Setting.txt` to specify your reference database, primers, and quality thresholds. If you need to customize tool-specific parameters, edit `Options_usearch.txt` or `Options_DADA2.txt`.
2. **Prepare Input**: Place your FASTQ or FASTQ.gz files into the `Run/` directory.
3. **Run the Pipeline**: Open a terminal (macOS/Linux) or PowerShell (Windows) in the PMiFish directory and execute:
   ```bash
   perl ./PMiFish.pl
   ```
  
## 📊 Output Results
All results are automatically organized and saved in the `Results/` directory. The pipeline generates structured outputs, from raw sequence processing to final taxonomic reports.

### Key Deliverables
- **Portal_[Date_Time].html**: A comprehensive, interactive HTML report. It serves as the main dashboard to review read counts, detection lists, and summary statistics in your web browser.
- **5-2 Summary_Table**: Tables summarizing results by taxonomic group and representative sequences. This is typically the primary file used for further statistical analysis.

### Logs and Records
- **log.html**: A visual summary of read count transitions (filtering steps) for each sample.
- **Setting_log.txt**: A complete record of the `Setting.txt` parameters used for the run.
- **Options_log files**: Copies of `Options_usearch.txt` or `Options_DADA2.txt` to ensure reproducibility of your analysis.

## 🤝 Acknowledgment
This work was supported by **JSPS KAKENHI Grant Number 23K05282**. This script was refined and released as part of the research project funded by this grant.

## License
[MIT] https://github.com/rogotoh/PMiFish3.0/blob/master/LICENSE
