#!/usr/bin/env Rscript
# DADA2 error-model training for PacBio HiFi 16S (cohort-level, runs once).
#
# Split-mode stage 2 of 4 (filter -> learn_errors -> denoise -> merge).
#
# This stage is deliberately kept whole-cohort: learnErrors consumes filtered
# files in order until --nbases is reached, so the resulting model depends on
# the file ordering. Sorting by sample name here reproduces the ordering the
# cohort-level dada2_pacbio.R got from list.files().
#
# Note on cost: with the default --nbases 1e8 and ~200 Mbases per PacBio HiFi
# sample, this stage reads well under one sample regardless of cohort size, so
# its runtime is effectively constant in the number of samples.
#
# Args (sibling-style long flags):
#   --input <dir>       directory containing *_filt.fastq.gz files
#   --nproc <n>         worker threads
#   --nbases <num>      learnErrors nbases cap (default 1e8, the dada2 default)
#   --randomize <bool>  learnErrors randomize (default FALSE, the dada2 default)
#   --band_size <int>   BAND_SIZE passed to learnErrors
#
# Outputs (cwd):
#   err.rds             the fitted error model, consumed by dada2_denoise.R
#   qualBins.txt        the quality bins extracted from the data (provenance)

suppressMessages(suppressWarnings({
    library(dada2)
    library(ShortRead)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) > 0) args[idx + 1] else default
}

input_dir <- get_arg(args, "--input")
nproc     <- as.integer(get_arg(args, "--nproc", "1"))
nbases    <- as.numeric(get_arg(args, "--nbases", "1e8"))
randomize <- as.logical(get_arg(args, "--randomize", "FALSE"))
band_size <- as.integer(get_arg(args, "--band_size", "32"))

if (is.null(input_dir)) stop("dada2_learn_errors.R: --input is required")
if (is.na(randomize))   stop("dada2_learn_errors.R: --randomize must be TRUE or FALSE")

filts <- list.files(input_dir, pattern = "_filt\\.fastq\\.gz$", full.names = TRUE)
if (length(filts) == 0) stop("dada2_learn_errors.R: no *_filt.fastq.gz in ", input_dir)

# Deterministic, locale-independent ordering by sample name. learnErrors walks
# the files in this order until the nbases cap is hit, so the order is part of
# the reproducibility contract shared with dada2_merge.R.
sample.names <- sub("_filt\\.fastq\\.gz$", "", basename(filts))
ord          <- order(sample.names, method = "radix")
filts        <- filts[ord]
sample.names <- sample.names[ord]
names(filts) <- sample.names

cat("Found", length(filts), "filtered samples\n")

# ---- qualBins (auto-extract from a sample of the first filtered file) ----
# Per benjjneb/dada2#2107: qualBins should reflect the binned quality scores
# actually present in the data. PacBio CCS binning is consistent within a run,
# so sampling one file is sufficient.
cat("Extracting qualBins from", filts[[1]], "\n")
qual_sample <- yield(FastqStreamer(filts[[1]], n = 10000))
qualBins    <- sort(unique(as.vector(as(quality(qual_sample), "matrix"))))
qualBins    <- qualBins[!is.na(qualBins)]
cat("qualBins:", paste(qualBins, collapse = ","), "\n")
writeLines(as.character(qualBins), "qualBins.txt")

# ---- learn errors ----
err <- learnErrors(filts,
                   nbases = nbases,
                   randomize = randomize,
                   errorEstimationFunction = makeBinnedQualErrfun(qualBins),
                   BAND_SIZE = band_size, multithread = nproc)

saveRDS(err, "err.rds")
cat("dada2_learn_errors.R done\n")
