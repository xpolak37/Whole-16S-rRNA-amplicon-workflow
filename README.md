# Whole 16S rRNA amplicon workflow

This repository includes the bioinformatic workflow for *long-read* 16S rRNA amplicon processing, sequenced using pacbio technology. All parameters are set to work properly on specific data from our lab. Before any usage, please check if it fits your data as well. 

## 🚀 Quick usage

Whole workflow can be executed by main_whole16s.sh, where you have to change paths to conda environments and to your data.

```bash
bash main_whole16s.sh
```

### 📥 Inputs

<u>REQUIRED</u>

**PATHS**

- input_dir="/path/to/raw/sequencing/data" - preferably in fastq.gz format	
- project_dir="/path/to/pre-created/output/directory" - in this directory, all outputs will be saved in pre-programmed structure	
- TMPDIR="/path/to/tmp" - some tools save intermediate files as temporary, however, it is important to have pre-alocated path and storage for this, in order to not get error message because of the storage	
- decipher_classifier="/path/to/decipher/classifier" - RData file which includes pre-trained IDTAXA classifier based on your used primers. See how to build such classifier on another tutorial	
- pacbio_16s_env="/path/to/conda/environment/" - conda environment, which have pre-installed all required software	
- blast_env="/path/to/conda/blast_env" - conda environment with pre-installed BLAST (it is better to keep it in different environment, as conda struggled with some dependencies)

<u>OPTIONAL</u>

**DADA2 (default) PARAMETERS**
- minQ=3		
- minLen=1000	
- maxLen=1600	
- maxN=0	
- maxEE=2	

**IDTAXA (default) PARAMETERS**	
- threshold=60	
- strand="both"

### 📤 Outputs
Multiple directories in *project_dir*:
- run_info - all important information about the reads, their statistics, bioinformatics processing and final tables that are ready to be analyzed  
- quality_raw - fastqc and multiqc quality reports of raw reads 
- preprocessed - preprocessed fastqs (with cutadapt) with quality reports
- decontaminated - decontaminated fastqs (after human and phix reads removal)
- denoised - denoised reads after DADA2 pipeline processing in form of ASV and TAXA tables

### 📦 Requirements
- phix genomic fasta
- fastqc
- multiqc
- hostile
- cutadapt
- R with preinstalled packages (ShortRead, dada2, Biostrings, reshape2, DECIPHER, tibble)
- DADA2
- pre-trained taxonomic classifier

## ℹ️ About
The workflow consists of multiple modules, all will be executed if not set otherwise:

---
### 🧹 **<u>1. Quality control</u>**

Runs **FastQC & MultiQC** to provide quality control report of all samples. Additionally, it creates a **QUICK SUMMARY** 
of adapter content, sequencing depth, overrepresented sequences and their BLAST hits, for overall summary of the run. Lastly, it runs a primer check, generating a report on the number of occurrences of specific primer sequences in each sample. This is particularly useful when using phased staggered primers and you want to verify whether their distribution is relatively uniform.

**Input**: raw fastqs  
**Output**:
- *quality_raw* folder with fastqc and multiqc data
- *primers_check* folder with forward and reverse primers summary
- *run_info* folder with:
     - raw_custom_summary.txt - adapter content, sequencing depth, overrepresented sequences
     - raw_multiqc_report.html - copy of the multiqc_report.html
     - reads_count_summary.txt - read counts
---    


### ⚙️ **<u>2. Preprocessing</u>**

Runs **cutadapt** for preprocessing the raw fastq files. It trims primers and Nextera Transposase Adapters. Generally, it is executed with the default parameters with the minor changes. See the command below.

```bash
cutadapt \
    -g ^<forward_primer_sequence> \
    -a <reverse_primer_sequence>$ \
    -a 'A{10}' -a 'G{10}' \
    --cores 20 \
    --discard-untrimmed \
    --revcomp \
    -o /path/to/output/sample_trimmed.fastq.gz 
    </path/to/input/fastq.gz>

```

Next, it runs **hostile** to remove human and phix contamination. As for human decontamination, it uses hostile's prepared indexes (human-t2t-hla-argos985-mycob140) that consist of T2T-CHM13v2.0 + IPD-IMGT/HLA v3.51 masked with 150mers for 985 FDA-ARGOS bacterial & 140 mycobacterial genomes. As for phix decontamination, it uses phiX174 genomic fasta file, which can be downloaded here: [phiX174](https://www.ncbi.nlm.nih.gov/nuccore/9626372)

See the executed commands below. 

```bash
# human decontamination
hostile clean \
--fastq1 sample_trimmed.fastq.gz \
--index human-t2t-hla-argos985-mycob140 \
--output /path/to/project_dir/decontaminated/human/

# phix decontamination
hostile clean \
--fastq1 sample_trimmed.clean.fastq.gz \
--index path/to/phix/fasta \
--output /path/to/project_dir/decontaminated/human_phix/
```

**Input**: raw fastqs  
**Output**: 
- *trimmed* folder - results of trimmomatic
- *decontaminated* folder:
    - *decontaminated/human/* folder - results of hostile's human removal
    - *decontmainated/human_phix/* folder - results of hostile's phix removal
- *run_info* folder with:
     - read_counts_summary.txt - updated read counts track with <u>trimmed</u> and <u>decontaminated</u> columns. 
     - trimmed_multiqc_report.html - copy of the multiqc_report.html containing fastqc reports of trimmed reads
---

### 🧬 **<u> 3. Denoising</u>**

Performs denosing using **DADA2** following this particular tutorial:  https://github.com/benjjneb/LRASManuscript/blob/master/LRASms_fecal.Rmd. 

- **DADA2** uses a model-based approach to correct sequencing errors by learning error rates from the data itself and distinguishing true biological sequences (amplicon sequence variants, ASVs) from errors. It performs quality filtering, error modeling, dereplication, and chimera removal, providing single-nucleotide resolution of microbial variants.

!!! With PacBio reads, DADA2 processes full-length amplicons in single reads, so paired-end merging is not needed. The pipeline uses PacBio-specific filtering and error modeling adapted to long, high-accuracy HiFi sequences!!!

---
### 🦠 **<u> 4. Taxonomic assignment</u>**

This part is also implemented in the denoising module, specifically, in the DADA2 R script.  **idtaxa()** from DECIPHER R package is implemented here.

**IdTaxa()** is a function for taxonomic classification of DNA sequences (typically 16S, 18S, or ITS) using a machine-learning approach based on discriminant analysis of sequence k-mer composition. It assigns taxonomy by comparing query sequences to a trained reference database, estimating confidence scores for each taxonomic rank.

It is important, that proper pre-trained classifier is used as input here (see tutorial at > TO BE ADDED)

