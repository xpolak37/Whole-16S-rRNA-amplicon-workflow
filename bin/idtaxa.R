#!/usr/bin/env Rscript
# Stub — real IDTAXA implementation arrives in the IDTAXA stage spec
# (which will reconcile the existing taxonomy/idtaxa.R against the sibling version).
# Positional args (matches IDTAXA module): asv_fasta, training_set, output_tsv [, conf_output_tsv]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("idtaxa.R: expected at least 3 positional args, got ", length(args))
out <- args[3]
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
file.create(out)
if (length(args) >= 4) file.create(args[4])
cat("[stub] idtaxa.R touched", out, "\n")
