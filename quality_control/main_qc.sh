#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
amplicon_16s_env=$1
conda_env_dir_blast=$2

# paths for data storage
path_input=$3
path_project_dir=$4
path_output=${path_project_dir}/quality_raw

# activate conda environment 
conda activate ${amplicon_16s_env}

# info for tools and versions txt
echo "main_qc.sh:" >> ${path_project_dir}/run_info/tools.txt

# making dirs
mkdir ${path_output}

# run fastqc
find ${path_input} -type f -name "*.fastq.gz" | parallel -j 10 "fastqc -o "${path_output}" -quiet -t 5"

## track version
fastqc --version >> ${path_project_dir}/run_info/tools.txt

# run multiqc
cd ${path_output}
multiqc . --interactive
## track version
multiqc --version >> ${path_project_dir}/run_info/tools.txt

cd multiqc_data

# QUICK SUMMARY
## Adapter content
echo "ADAPTER CONTENT:" >  custom_summary.txt
awk 'NR>1 {print $3}' fastqc_adapter_content_plot.txt | sort -u >> custom_summary.txt

## sequencing septh
conda activate ${conda_env_dir_R}
echo "SEQUENCING DEPTH SUMMARY" >> custom_summary.txt
awk 'NR>1 {print $NF*1000000}' multiqc_general_stats.txt | Rscript -e 'x <- scan("stdin", quiet=TRUE); summary(x)' >> custom_summary.txt
### track version
R --version | head -n 1 >> ${path_project_dir}/run_info/tools.txt

## overrepresented sequences
conda activate ${conda_env_dir_blast}

### get overrepresented sequences to fasta
awk 'NR>1 {print ">seq" NR-1 "\n" $1}' fastqc_top_overrepresented_sequences_table.txt > overrepresented.fasta

### run blastn
blastn -query overrepresented.fasta -db nt -remote -out overrepresented_results.txt -outfmt 6 
### track version
blastn -version | head -n 1 >> ${path_project_dir}/run_info/tools.txt

### filter results
awk '{
    if (!($1 in best) || $12 > best[$1]) {
        best[$1] = $12
        line[$1] = $0
    }
}
END {
    for (q in line) {
        print line[q]
    }
}' overrepresented_results.txt > best_hits.txt

### change the headers to identificants
# First: build a mapping from BLAST results (seqID -> accession)
awk '{print $1, $2}' best_hits.txt | sort -u > mapping.txt

# Now rewrite the fasta
awk '
BEGIN {
    # read mapping into array
    while ((getline < "mapping.txt") > 0) {
        map[$1] = $2
    }
}
{
    if ($0 ~ /^>/) {
        seqid = substr($0,2)   # remove ">"
        if (seqid in map) {
            # fetch NCBI header
            cmd = "esearch -db nucleotide -query \"" map[seqid] "\" | efetch -format fasta | head -n1"
            cmd | getline newheader
            close(cmd)
            print newheader
        } else {
            # no BLAST hit → keep original header
            print $0
        }
    } else {
        # sequence line
        print $0
    }
}' overrepresented.fasta > overrepresented_named.fasta

echo "OVERREPRESENTED SEQUENCES:" >> custom_summary.txt
cat overrepresented_named.fasta >> custom_summary.txt

### remove all redundant files
rm overrepresented.fasta overrepresented_named.fasta mapping.txt best_hits.txt overrepresented_results.txt

# The final counts statistics
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt

# Header
echo -e "Sample\tRaw\tTrimmed\tDecontaminated" > ${outfile}
awk 'NR>1 {printf "%s\t%.0f\n", $1, $7*1000000}' ${path_output}/multiqc_data/multiqc_general_stats.txt >> ${outfile}

## LOW SEQUENCING DEPTH SAMPLES
echo "SAMPLES BELOW 10k:" >> custom_summary.txt
awk 'NR > 1 && $2 < 10000' ${path_project_dir}/run_info/read_counts_summary.txt >> custom_summary.txt

# copy the summaries to the run_info file
cp custom_summary.txt ${path_project_dir}/run_info/raw_custom_summary.txt
cp ../multiqc_report.html ${path_project_dir}/run_info/raw_multiqc_report.html
