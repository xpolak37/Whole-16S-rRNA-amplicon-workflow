#!/usr/bin/env Rscript
# DADA2 dereplication + denoising for a SINGLE PacBio HiFi 16S sample.
#
# Split-mode stage 3 of 4 (filter -> learn_errors -> denoise -> merge).
#
# Correctness note: dada() defaults to pool = FALSE, which denoises every
# sample independently of the others. The cohort-level dada2_pacbio.R relied on
# that same default, so denoising one sample per task yields the same partition
# as the cohort-level call -- this stage is a scheduling change, not a
# statistical one.
#
# This is also the stage that made the cohort-level script blow up on memory:
# it held a derep object (unique sequences x quality matrix) for every sample at
# once. Here peak memory is one sample's worth regardless of cohort size.
#
# Args (sibling-style long flags):
#   --input <fq>        one *_filt.fastq.gz file
#   --sample <name>     sample id used for output naming
#   --err <rds>         error model from dada2_learn_errors.R
#   --nproc <n>         worker threads
#   --band_size <int>   BAND_SIZE passed to dada
#
# Outputs (cwd):
#   <sample>.uniques.tsv        Sequence, Abundance (the denoised ASVs)
#   <sample>.denoise_stats.tsv  SampleID, denoised

suppressMessages(suppressWarnings({
    library(dada2)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) > 0) args[idx + 1] else default
}

input_fq  <- get_arg(args, "--input")
sample    <- get_arg(args, "--sample")
err_rds   <- get_arg(args, "--err")
nproc     <- as.integer(get_arg(args, "--nproc", "1"))
band_size <- as.integer(get_arg(args, "--band_size", "32"))

if (is.null(input_fq)) stop("dada2_denoise.R: --input is required")
if (is.null(sample))   stop("dada2_denoise.R: --sample is required")
if (is.null(err_rds))  stop("dada2_denoise.R: --err is required")

err <- readRDS(err_rds)

drp <- derepFastq(input_fq, verbose = FALSE)
dd  <- dada(drp, err = err, BAND_SIZE = band_size, multithread = nproc)

uniques <- getUniques(dd)

out <- data.frame(
    Sequence  = names(uniques),
    Abundance = as.integer(uniques),
    row.names = NULL,
    check.names = FALSE
)
write.table(out, file = paste0(sample, ".uniques.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

stats <- data.frame(
    SampleID = sample,
    denoised = as.integer(sum(uniques)),
    row.names = NULL,
    check.names = FALSE
)
write.table(stats, file = paste0(sample, ".denoise_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_denoise.R:", sample, "-", nrow(out), "ASVs,",
    sum(uniques), "reads\n")
