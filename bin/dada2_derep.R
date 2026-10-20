#!/usr/bin/env Rscript
# DADA2 dereplication for a SINGLE PacBio HiFi 16S sample (no-denoise path).
#
# Split-mode stage 2 of 3 for dada2_nodenoise (filter -> derep -> merge).
#
# derepFastq is per-file and independent, so running it one sample per task is
# equivalent to the cohort-level call in dada2_pacbio_nodenoise.R.
#
# Unlike the denoised path there is no dada() step to collapse the output, so
# these files hold every raw unique sequence and are large. They are gzipped for
# that reason, and dada2_nodenoise_merge.R streams them rather than loading the
# whole cohort at once.
#
# Args (sibling-style long flags):
#   --input <fq>        one *_filt.fastq.gz file
#   --sample <name>     sample id used for output naming
#
# Outputs (cwd):
#   <sample>.uniques.tsv.gz    Sequence, Abundance (raw uniques)
#   <sample>.derep_stats.tsv   SampleID, reads

suppressMessages(suppressWarnings({
    library(dada2)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) > 0) args[idx + 1] else default
}

input_fq <- get_arg(args, "--input")
sample   <- get_arg(args, "--sample")

if (is.null(input_fq)) stop("dada2_derep.R: --input is required")
if (is.null(sample))   stop("dada2_derep.R: --sample is required")

drp     <- derepFastq(input_fq, verbose = FALSE)
uniques <- getUniques(drp)

out <- data.frame(
    Sequence  = names(uniques),
    Abundance = as.integer(uniques),
    row.names = NULL,
    check.names = FALSE
)
gz <- gzfile(paste0(sample, ".uniques.tsv.gz"), "w")
write.table(out, file = gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

# Matches the cohort script's getN(): sum of the uniques vector, i.e. reads
# surviving to dereplication, not the number of distinct sequences.
stats <- data.frame(
    SampleID = sample,
    reads    = as.integer(sum(uniques)),
    row.names = NULL,
    check.names = FALSE
)
write.table(stats, file = paste0(sample, ".derep_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_derep.R:", sample, "-", nrow(out), "unique sequences,",
    sum(uniques), "reads\n")
