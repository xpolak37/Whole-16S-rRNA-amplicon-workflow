#!/bin/bash

eval "$(conda shell.bash hook)"

####### USER INPUTS #########
# DIRECTORIES
export input_dir="/home/povp/seq_data/16S/pacbio_novogene_jan26/data/pool1/"
export project_dir="/home/povp/Projects/pacbio_novogene_jan26/pool1/"
export TMPDIR="/home/povp/tmp"
export decipher_classifier="/home/povp/taxonomic_classifiers/decipher_classifier/idtaxa_trainingSet_V3V4_silva_138_2.RData"

# CONDA ENVIRONMENTS
export pacbio_16s_env="/home/povp/conda_envs/pacbio_16/"
export blast_env="/home/povp/conda_envs/blast_env"

# PARAMETERS
## DADA
minQ=3
minLen=1000
maxLen=1600
maxN=0
maxEE=2

## IDTAXA
export threshold=60
export strand="both"

##############################

conda activate ${pacbio_16s_env}

perform_quality_control=TRUE
perform_preprocessing=TRUE
perform_denoising=TRUE
perform_taxassignment=TRUE
keep_intermediate=FALSE

# set all variables
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir=$2; shift 2 ;;
        --project_dir) project_dir=$2; shift 2 ;;
        --keep_intermediate) keep_intermediate=TRUE; shift ;;
        --skip_quality) perform_quality_control=FALSE; shift ;;
        --skip_preprocessing) perform_preprocessing=FALSE; shift ;;
        --skip_denoising) perform_denoising=FALSE; shift ;;
        --skip_taxassignment) perform_taxassignment=FALSE; shift ;;
        --blast_env) blast_env=$2; shift 2;;
        --main_env) pacbio_16s_env=$2; shift 2;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

export path_scripts="$(pwd)"

# define the txt file for version tracking
mkdir -p ${project_dir}/run_info/
echo "# Tools and their versions throughout the run" > ${project_dir}/run_info/tools.txt

# RUN the pipeline

## Quality control
if [[ "${perform_quality_control}" == "TRUE" ]]; then
    echo "Performing quality control..."
    cd quality_control
    bash main_qc.sh ${pacbio_16s_env} ${blast_env} ${input_dir} ${project_dir}
    conda activate ${pacbio_16s_env}
    cd ${path_scripts}
fi

# Preprocessing
if [[ "${perform_preprocessing}" == "TRUE" ]]; then
    cd preprocessing
    echo "Performing trimming with cutadapt"
    bash cutadapt.sh ${pacbio_16s_env} ${input_dir} ${project_dir} 
    echo "Performing host decontamination with hostile"
    bash host_decontamination.sh ${pacbio_16s_env} ${project_dir}/preprocessed/ ${project_dir} 
    echo "Done"
    cd ${path_scripts}
fi

# Denoising
if [[ "${perform_denoising}" == "TRUE" ]]; then
    cd denoising

    mkdir ${project_dir}/denoised
    echo "Running Dada2 for denoising ..."
    Rscript main_dada2.R ${project_dir}/decontaminated/human_phix ${project_dir}/denoised ${minQ} ${minLen} ${maxLen} ${maxN} ${maxEE} ${decipher_classifier}
    echo Done
    
    # cp the result in run_info
    cp ${project_dir}/denoised/ASV_table.csv ${project_dir}/run_info/ASV_table.tsv
    cp ${project_dir}/denoised/taxa_table.csv ${project_dir}/run_info/taxa_table.tsv
    cp ${project_dir}/denoised/taxa_table_conf.csv ${project_dir}/run_info/taxa_table_conf.csv

    cd ${path_scripts}
fi


