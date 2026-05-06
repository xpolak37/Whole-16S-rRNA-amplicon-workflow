#!/usr/bin/env Rscript
# Stub — real implementation arrives in the DADA2 stage spec.
# Positional args (matches DADA2_PACBIO_NODENOISE):
#   1: input_dir, 2: output_dir, 3: minQ, 4: minLen, 5: maxLen,
#   6: maxN, 7: maxEE, 8: minAbundance
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) stop("dada2_pacbio_nodenoise.R: expected 8 positional args, got ", length(args))
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
file.create(file.path(output_dir, "ASV_table.tsv"))
file.create(file.path(output_dir, "ASV_sequences.fasta"))
file.create(file.path(output_dir, "track_control.tsv"))
cat("[stub] dada2_pacbio_nodenoise.R touched outputs in", output_dir, "\n")
