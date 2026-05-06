#!/bin/bash

eval "$(conda shell.bash hook)"

####### USER INPUTS #########
# DIRECTORIES
export input_dir="/home/povp/Projects/g363_longreads/plzen_library/multiplexed_analysis/trimmed/"
export project_dir="/home/povp/Projects/g363_longreads/plzen_library/multiplexed_analysis/"
export TMPDIR="/home/povp/tmp"
export decipher_classifier="/home/povp/taxonomic_classifiers/decipher_classifier/idtaxa_trainingSet_W16S_silva_138_2.RData"
export orienting_db="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_W16S/silva-138.2-ssu-nr99-seqs-27F-1492R-orienting.fasta"
export vsearch_db="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_W16S/silva-138.2-ssu-nr99-seqs-27F-1492R-uniq.qza"
export vsearch_tax="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_W16S/silva-138.2-ssu-nr99-tax-27F-1492R-derep-uniq.qza"

# CONDA ENVIRONMENTS
export pacbio_16s_env="/home/povp/conda_envs/pacbio_16s/"
export blast_env="/home/povp/conda_envs/blast_env"
export qiime_env="/home/povp/conda_envs/qiime2-amplicon-2026.1/"

# PARAMETERS
## DADA
minQ=0
minLen=1000
maxLen=1600
maxN=0
maxEE=Inf
minAbundance=0.001

## VSEARCH - copying default from https://github.com/PacificBiosciences/HiFi-16S-workflow/blob/main/nextflow.config
vsearch_identity=0.95
vsearch_maxreject=100
vsearch_maxaccept=100
vsearch_pminconsensus=0.51

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
perform_trimming=TRUE
perform_hostremoval=TRUE
perform_dada=TRUE
perform_taxonomy=TRUE
perform_vsearch=FALSE


# HELP OPTIONS
show_help() {
    cat << EOF
Usage: bash main.sh [OPTIONS]

Description:
    Bioinformatic workflow for long-read 16S rRNA amplicon processing, sequenced using PacBio technology.

Options:
    Required:
    -i,  --input_dir            Path to input folder containing raw fastq.gz files
    -o,  --project_dir          Path to project directory where all results will be saved

    Environments:
    -m,  --main_env             Name of the main conda environment used to run the pipeline
    -b,  --blast_env            Name of the conda environment used for BLAST taxonomy classification

    Pipeline control:
         --skip_quality         Skip the quality control step
         --skip_preprocessing   Skip the preprocessing step
         --skip_trimming        Skip the primer and adapter trimming step
         --skip_hostremoval     Skip the host read removal step
         --skip_denoising       Skip the denoising step
         --skip_dada            Skip the DADA2 ASV inference step
         --skip_taxonomy        Skip the taxonomy classification step

    Additional:
         --keep_intermediate    Keep intermediate files generated during the pipeline run
         --vsearch              Use CONSENSUS-VSEARCH taxonomy assignment instead of IDTAXA
    -h,  --help                 Show this help message and exit

Example:
    bash main.sh --input_dir /data/raw --project_dir /data/results --main_env pacbio_env --blast_env blast_env
EOF
}


# If no arguments provided
if [[ $# -eq 0 ]]; then
    show_help
    exit 1
fi

# set all variables
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir) input_dir=$2; shift 2 ;;
        --project_dir) project_dir=$2; shift 2 ;;
        --keep_intermediate) keep_intermediate=TRUE; shift ;;
        --skip_quality) perform_quality_control=FALSE; shift ;;
        --skip_preprocessing) perform_preprocessing=FALSE; shift ;;
        --skip_trimming) perform_trimming=FALSE; shift ;;
        --skip_hostremoval) perform_hostremoval=FALSE; shift ;;
        --skip_denoising) perform_denoising=FALSE; shift ;;
        --skip_dada) perform_dada=FALSE; shift ;;
        --skip_taxonomy) perform_taxonomy=FALSE; shift ;;
        --blast_env) blast_env=$2; shift 2;;
        --main_env) pacbio_16s_env=$2; shift 2;;
        --vsearch) perform_vsearch=TRUE; shift ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
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
    cd ${path_scripts}
fi


# Preprocessing
if [[ "${perform_preprocessing}" == "TRUE" ]]; then
    cd preprocessing
    if [[ "${perform_trimming}" == "TRUE" ]]; then
        echo "Performing trimming with cutadapt..."
        bash cutadapt.sh ${pacbio_16s_env} ${input_dir} ${project_dir} 
    fi

    if [[ "${perform_hostremoval}" == "TRUE" ]]; then
        if [[ "${perform_trimming}" == "FALSE" ]]; then
            mkdir ${project_dir}/preprocessed/
            cp ${input_dir}/*fastq.gz ${project_dir}/preprocessed/
        fi

        echo "Performing host decontamination with hostile"
        bash host_decontamination.sh ${pacbio_16s_env} ${project_dir}/preprocessed/ ${project_dir} 
    fi

    echo "Done"
    cd ${path_scripts}
fi

# Denoising
if [[ "${perform_dada}" == "TRUE" ]]; then
    cd denoising
    mkdir ${project_dir}/denoised

    bash denoise_preprocessing.sh ${pacbio_16s_env} ${project_dir}/decontaminated/human_phix ${project_dir}
    
    conda activate ${pacbio_16s_env}

    if [[ "${perform_denoising}" == "TRUE" ]]; then
        echo "Running Dada2 for denoising ..."
        Rscript main_dada2.R ${project_dir}/oriented ${project_dir}/denoised ${minQ} ${minLen} ${maxLen} ${maxN} ${maxEE}
        echo Done
    fi

     if [[ "${perform_denoising}" == "FALSE" ]]; then
        echo "Running Dada2 without denoising ..."
        Rscript main_dada2.R ${project_dir}/oriented ${project_dir}/denoised ${minQ} ${minLen} ${maxLen} ${maxN} ${maxEE} ${minAbundance}
        echo Done
    fi

    # cp the result in run_info
    cp ${project_dir}/denoised/ASV_table.tsv ${project_dir}/run_info/ASV_table.tsv

    cd ${path_scripts}
fi

if [[ "${perform_taxonomy}" == "TRUE" ]]; then
    cd taxonomy
    mkdir ${project_dir}/taxonomy

    if [[ "${perform_vsearch}" == "TRUE" ]]; then
        echo "Running qiime feature-classifier-vsearch for taxonomic classification ..."
        bash qiime_vsearch.sh ${qiime_env} ${project_dir}/denoised ${project_dir} ${vsearch_db} ${vsearch_tax} ${vsearch_identity} ${vsearch_maxreject} ${vsearch_maxaccept} ${vsearch_pminconsensus}
        echo Done
    fi

    else
        echo "Running DECIPHER::IDTAXA for taxonomic classification ..."
        Rscript idtaxa.R ${project_dir}/denoised ${project_dir}/denoised ${decipher_classifier}
        echo Done
    fi

    # cp the result in run_info
    cp ${project_dir}/denoised/taxa_table.tsv ${project_dir}/run_info/taxa_table.tsv
    cp ${project_dir}/denoised/taxa_table_conf.tsv ${project_dir}/run_info/taxa_table_conf.csv


    cd ${path_scripts}
fi

