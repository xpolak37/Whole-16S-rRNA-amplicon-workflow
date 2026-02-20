#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
pacbio_16s_env=$1

# paths for data storage
path_input=$2
path_project_dir=$3
export path_output=${path_project_dir}/oriented

conda activate ${pacbio_16s_env}

mkdir ${path_output}

# info for tools and versions txt
echo -e "\ndenoise_preprocessing.sh:" >> ${path_project_dir}/run_info/tools.txt

# ORIENTING READS
run_vsearch(){
        reads=$1
        sample=$(echo "$reads" | sed 's/.fastq.gz//')
        sample=$(basename ${sample})
        vsearch --orient ${reads} --db ${orienting_db} --fastqout ${path_output}/${sample}-oriented.fq 
}

export -f run_vsearch

find ${path_input} -type f -name "*.fastq.gz" | parallel -j 10 run_vsearch

## track version
vsearch --version >> ${path_project_dir}/run_info/tools.txt
