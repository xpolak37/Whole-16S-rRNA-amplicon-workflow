#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
amplicon_16s_env=$1

# paths for data storage
path_input=$2
path_project_dir=$3
path_output=${path_project_dir}/preprocessed

# activate environment
conda activate ${amplicon_16s_env}

# info for tools and versions txt
echo -e "\ncutadapt.sh:" >> ${path_project_dir}/run_info/tools.txt

# create output dirs
mkdir ${path_output}
mkdir ${path_output}/reports/

# adapters + primers trimming 
for read in ${path_input}/*.fastq.gz
do
    sample_name=$(basename "${read}" | sed 's/.fastq.gz//')
	cutadapt \
    -g ^AGRGTTYGATYMTGGCTCAG \
    -a RGYTACCTTGTTACGACTT$ \
    -a 'A{10}' -a 'G{10}' \
	--cores 20 \
    --discard-untrimmed \
    --revcomp \
    -o ${path_output}/${sample_name}_trimmed.fastq.gz \
    ${read}
	
done

## track version
echo "cutadapt " $(cutadapt --version) >> ${path_project_dir}/run_info/tools.txt

# FASTQC+MULTIQC for trimmed report
find ${path_output} -type f -name "*.fastq.gz" | parallel -j 10 "fastqc -o "${path_output}/reports/" -quiet -t 5"
cd ${path_output}/reports/
multiqc . --interactive

## track version
fastqc --version >> ${path_project_dir}/run_info/tools.txt
multiqc --version >> ${path_project_dir}/run_info/tools.txt

# The final counts statistics
## Output file
outfile=${path_project_dir}/run_info/read_counts_summary.txt
infile=${path_output}/reports/multiqc_data/multiqc_general_stats.txt

cd ${path_scripts}/preprocessing
python3 counts_summary.py ${infile} ${outfile} 2

## track version
python --version >> ${path_project_dir}/run_info/tools.txt

# copy the summaries to the run_info file
cp ${path_output}/reports/multiqc_report.html ${path_project_dir}/run_info/trimmed_multiqc_report.html
