#!/bin/bash
eval "$(conda shell.bash hook)"

# variables for conda environment
qiime_env=$1

# paths for data storage
path_input=$2
path_project_dir=$3
path_output=${path_project_dir}/taxonomy
vsearch_identity=$4
vsearch_maxreject=$5
vsearch_maxaccept=$6

# info for tools and versions txt
echo -e "\nqiime_vsearch.sh:" >> ${path_project_dir}/run_info/tools.txt

conda activate ${qiime_env}
mkdir ${path_output}

# IMPORTING
qiime tools import \
  --input-path ${path_input}/ASV_sequences.fasta \
  --output-path ${path_output}/rep-seqs.qza \
  --type 'FeatureData[Sequence]'

cd ${path_output}

# TAXONOMY
qiime feature-classifier classify-consensus-vsearch --i-query rep-seqs.qza \
        --o-classification taxonomy.vsearch.qza \
        --i-reference-reads ${vsearch_db} \
        --i-reference-taxonomy ${vsearch_tax} \
        --p-threads 20 \
        --p-maxrejects ${vsearch_maxreject} \
        --p-maxaccepts ${vsearch_maxaccept} \
        --p-perc-identity ${vsearch_identity} \
        --p-strand both \
        --p-top-hits-only \
        --o-search-results taxonomy.vsearch.blast_results.qza \
        --p-min-consensus ${vsearch_pminconsensus}

qiime tools export --input-path taxonomy.vsearch.qza --output-path tax_export

qiime tools export \
--input-path taxonomy.vsearch.qza \
--output-path qiime_results
