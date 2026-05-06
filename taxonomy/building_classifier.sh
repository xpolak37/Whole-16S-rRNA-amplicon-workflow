#!/bin/bash

conda activate /home/povp/conda_envs/qiime2-amplicon-2026.1

export TMPDIR="/home/povp/tmp/"

qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-include-species-labels

qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

# the quality control filter
qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

# min len FILTER
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1000 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza 

# Qiime rescript dereplicate, the --p-mode determines how the plugin handles taxonomy when 
# it finds multiple sequences that are 100% identical but have different taxonomic labels.
# When you "dereplicate," you are collapsing all identical sequences into a 
# single representative sequence to make your database smaller and faster. However, if Sequence A and Sequence B are identical, 
# but Sequence A is labeled "Genus Escherichia" and Sequence B is labeled "Genus Shigella", the mode tells RESCRIPt how to resolve that conflict.

# choices('uniq', 'lca', 'majority', 'super')
# "uniq" will retain all sequences with unique taxonomic affiliations. 
# "lca" will find the least common ancestor among all taxa sharing a sequence. 
# "majority" will find the most common taxonomic label associated with that sequence; note that in the event of a tie, "majority" will pick the winner arbitrarily. 
# "super" finds the LCA consensus while giving preference to majority labels and collapsing substrings into superstrings.

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.2-ssu-nr99-tax.qza \
    --p-rank-handles 'domain' 'phylum' 'class' 'order' 'family' 'genus' 'species' \
    --p-mode 'super' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-super.qza

f_primer = "CCTACGGGNGGCWGCAG"
    r_primer = "GACTACHVGGGTATCTAATCC"
    
qiime feature-classifier extract-reads \
    --i-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
    --p-f-primer AGRGTTYGATYMTGGCTCAG \
    --p-r-primer RGYTACCTTGTTACGACTT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.2-ssu-nr99-seqs-27F-1492R.qza

qiime feature-classifier extract-reads \
    --i-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
    --p-f-primer CCTACGGGAGGCAGCAG \
    --p-r-primer GGACTACHVGGGTWTCTAAT \
    --p-n-jobs 10 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.2-ssu-nr99-seqs-341F-805R.qza

# second dereplication 
# uniq mode: 
# Collapses identical sequences only
# Keeps taxonomy exactly as-is
# Does NOT attempt taxonomy merging or conflict resolution

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-27F-1492R.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-27F-1492R-uniq.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-27F-1492R-derep-uniq.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-341F-805R.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-341F-805R-uniq.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-341F-805R-derep-uniq.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-27F-1492R.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
    --p-mode 'lca' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-27F-1492R-lca.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-27F-1492R-derep-lca.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-341F-805R.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
    --p-mode 'lca' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-341F-805R-lca.qza \
    --o-dereplicated-taxa  silva-138.2-ssu-nr99-tax-341F-805R-derep-lca.qza

# Now: 
# vsearch_db="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_W16S/silva-138.2-ssu-nr99-seqs-27F-1492R-uniq.qza"
# vsearch_tax="/home/povp/taxonomic_classifiers/Rescript_classifier/classifier_W16S/silva-138.2-ssu-nr99-tax-27F-1492R-derep-uniq.qza"

# Naive bayes
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.2-ssu-nr99-seqs-27F-1492R-uniq.qza \
    --i-reference-taxonomy silva-138.2-ssu-nr99-tax-27F-1492R-derep-uniq.qza \
    --o-classifier silva-138.2-ssu-nr99-27F-1492R-classifier.qza