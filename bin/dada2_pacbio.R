#!/usr/bin/env Rscript
# Stub — real implementation arrives in the DADA2 stage spec.
# Positional args (matches calling module DADA2_PACBIO):
#   1: input_dir (oriented fastqs)
#   2: output_dir
#   3: minQ
#   4: minLen
#   5: maxLen
#   6: maxN
#   7: maxEE  (string 'Inf' allowed)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("dada2_pacbio.R: expected 7 positional args, got ", length(args))
input_dir <- args[1]; output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
file.create(file.path(output_dir, "ASV_table.tsv"))
file.create(file.path(output_dir, "ASV_sequences.fasta"))
file.create(file.path(output_dir, "track_control.tsv"))
cat("[stub] dada2_pacbio.R touched outputs in", output_dir, "\n")
