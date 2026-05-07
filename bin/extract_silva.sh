#!/bin/bash
#
# Build the SILVA-derived QIIME2 artefacts feeding the four classifiers
# (qnb, qblast, idtaxa, assigntaxonomy) for full-length 16S (PacBio HiFi).
#
# Run this manually inside a QIIME2 amplicon environment when provisioning
# a fresh --classifiers_dir. The Nextflow pipeline does not invoke this; it
# expects the resulting artefacts to already live in --classifiers_dir.
#
# Primers anchored to 27F-1492R for full-length 16S (matches main.nf
# default f_primer/r_primer).
#
# Usage:
#   conda activate /path/to/qiime2-amplicon-2026.1
#   bash extract_silva.sh

set -euo pipefail

qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-include-species-labels

qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1000 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza \
    --i-taxa silva-138.2-ssu-nr99-tax.qza \
    --p-rank-handles 'domain' 'phylum' 'class' 'order' 'family' 'genus' 'species' \
    --p-mode 'super' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-super.qza

# Full-length 16S extraction (27F / 1492R)
qiime feature-classifier extract-reads \
    --i-sequences silva-138.2-ssu-nr99-seqs-derep-super.qza \
    --p-f-primer AGRGTTYGATYMTGGCTCAG \
    --p-r-primer RGYTACCTTGTTACGACTT \
    --p-n-jobs 10 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.2-ssu-nr99-seqs-27F-1492R.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-27F-1492R.qza \
    --i-taxa silva-138.2-ssu-nr99-tax-derep-super.qza \
    --p-mode 'lca' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-27F-1492R-lca.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-27F-1492R-derep-lca.qza

# Naive Bayes classifier (qnb_classifier.qza)
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.2-ssu-nr99-seqs-27F-1492R-lca.qza \
    --i-reference-taxonomy silva-138.2-ssu-nr99-tax-27F-1492R-derep-lca.qza \
    --o-classifier qnb_classifier.qza

# QIIME consensus-BLAST inputs (qblast_seqs.qza / qblast_tax.qza)
cp silva-138.2-ssu-nr99-seqs-27F-1492R-lca.qza qblast_seqs.qza
cp silva-138.2-ssu-nr99-tax-27F-1492R-derep-lca.qza qblast_tax.qza

# IDTAXA / assignTaxonomy artefacts are built from the exported sequences +
# taxonomy via bin/train_classifiers.py (run after this script).
qiime tools export --input-path silva-138.2-ssu-nr99-seqs-27F-1492R-lca.qza --output-path .
qiime tools export --input-path silva-138.2-ssu-nr99-tax-27F-1492R-derep-lca.qza --output-path .
