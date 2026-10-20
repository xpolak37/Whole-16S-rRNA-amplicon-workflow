#!/usr/bin/env Rscript
# DADA2 filterAndTrim for a SINGLE PacBio HiFi 16S sample.
#
# Split-mode stage 1 of 4 (filter -> learn_errors -> denoise -> merge).
# filterAndTrim treats every input file independently, so running it one
# sample per task produces byte-identical filtered reads to the cohort-level
# call in dada2_pacbio.R.
#
# Args (sibling-style long flags):
#   --input <fq>        one *-oriented.fq file
#   --sample <name>     sample id used for output naming
#   --nproc <n>         worker threads
#   --minQ <int>        filterAndTrim minQ
#   --minLen <int>      filterAndTrim minLen
#   --maxLen <int>      filterAndTrim maxLen
#   --maxN <int>        filterAndTrim maxN
#   --maxEE <num|Inf>   filterAndTrim maxEE
#
# Outputs (cwd):
#   <sample>_filt.fastq.gz     filtered reads (absent if nothing survived)
#   <sample>.filter_stats.tsv  SampleID, input, filtered

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
nproc     <- as.integer(get_arg(args, "--nproc", "1"))
minQ      <- as.integer(get_arg(args, "--minQ",   "0"))
minLen    <- as.integer(get_arg(args, "--minLen", "1000"))
maxLen    <- as.integer(get_arg(args, "--maxLen", "1600"))
maxN      <- as.integer(get_arg(args, "--maxN",   "0"))
maxEE_raw <- get_arg(args, "--maxEE", "Inf")
maxEE     <- if (identical(maxEE_raw, "Inf")) Inf else as.numeric(maxEE_raw)

if (is.null(input_fq)) stop("dada2_filter.R: --input is required")
if (is.null(sample))   stop("dada2_filter.R: --sample is required")

filt <- paste0(sample, "_filt.fastq.gz")

out <- filterAndTrim(input_fq, filt,
                     minQ = minQ, minLen = minLen, maxLen = maxLen,
                     maxN = maxN, maxEE = maxEE, rm.phix = FALSE,
                     multithread = nproc, verbose = FALSE)

stats <- data.frame(
    SampleID = sample,
    input    = as.integer(out[1, 1]),
    filtered = as.integer(out[1, 2]),
    row.names = NULL,
    check.names = FALSE
)
write.table(stats, file = paste0(sample, ".filter_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# filterAndTrim does not create the output file when no read survives. The
# cohort-level script crashed later in derepFastq on the missing path; here the
# file is simply absent, the Nextflow output is optional, and dada2_merge.R
# backfills the sample with zero counts.
if (!file.exists(filt)) {
    cat("dada2_filter.R:", sample, "- no reads survived filtering\n")
} else {
    cat("dada2_filter.R:", sample, "-", out[1, 2], "/", out[1, 1], "reads kept\n")
}
