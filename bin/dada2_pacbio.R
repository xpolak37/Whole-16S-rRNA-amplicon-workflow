#!/usr/bin/env Rscript
# DADA2 denoising for PacBio HiFi 16S (single-end, cohort-level).
#
# Args (sibling-style long flags):
#   --input <dir>       directory containing *-oriented.fq files
#   --nproc <n>         worker threads
#   --minQ <int>        filterAndTrim minQ
#   --minLen <int>      filterAndTrim minLen
#   --maxLen <int>      filterAndTrim maxLen
#   --maxN <int>        filterAndTrim maxN
#   --maxEE <num|Inf>   filterAndTrim maxEE
#
# Outputs (cwd):
#   ASV_table.tsv          SeqID + per-sample abundance columns
#   ASV_sequences.fasta    >ASV_<i> headers
#   track_control.tsv      SampleID, input, filtered, denoised, nonchim
#   dada2.counts.tsv       sample\tcount  (count = nonchimeric reads, one row per sample)

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
minQ      <- as.integer(get_arg(args, "--minQ",   "0"))
minLen    <- as.integer(get_arg(args, "--minLen", "1000"))
maxLen    <- as.integer(get_arg(args, "--maxLen", "1600"))
maxN      <- as.integer(get_arg(args, "--maxN",   "0"))
maxEE_raw <- get_arg(args, "--maxEE", "Inf")
maxEE     <- if (identical(maxEE_raw, "Inf")) Inf else as.numeric(maxEE_raw)

if (is.null(input_dir)) stop("dada2_pacbio.R: --input is required")

reads_path   <- list.files(input_dir, pattern = "-oriented\\.fq$", full.names = TRUE)
if (length(reads_path) == 0) stop("dada2_pacbio.R: no *-oriented.fq files in ", input_dir)
sample.names <- sub("-oriented\\.fq$", "", basename(reads_path))
sample.names <- sub("_trimmed_cleaned$", "", sample.names)

cat("Found", length(reads_path), "samples\n")

# ---- filter ----
filts <- file.path("filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filts) <- sample.names
out <- filterAndTrim(reads_path, filts,
                     minQ = minQ, minLen = minLen, maxLen = maxLen,
                     maxN = maxN, maxEE = maxEE, rm.phix = FALSE,
                     multithread = nproc, verbose = FALSE)
rownames(out) <- sample.names

# ---- qualBins (auto-extract from a sample of the first filtered file) ----
# Per benjjneb/dada2#2107: qualBins should reflect the binned quality scores
# actually present in the data. PacBio CCS binning is consistent within a run,
# so sampling one file is sufficient.
cat("Extracting qualBins from", filts[[1]], "\n")
qual_sample <- yield(FastqStreamer(filts[[1]], n = 10000))
qualBins    <- sort(unique(as.vector(as(quality(qual_sample), "matrix"))))
qualBins    <- qualBins[!is.na(qualBins)]
cat("qualBins:", paste(qualBins, collapse = ","), "\n")

# ---- learn errors + denoise ----
err <- learnErrors(filts,
                   errorEstimationFunction = makeBinnedQualErrfun(qualBins),
                   BAND_SIZE = 32, multithread = nproc)
drp <- derepFastq(filts, verbose = FALSE)
dd  <- dada(drp, err = err, BAND_SIZE = 32, multithread = nproc)
if (inherits(dd, "dada")) dd <- setNames(list(dd), sample.names)

# ---- ASV table + chimeras ----
seqtab        <- makeSequenceTable(dd)
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    minFoldParentOverAbundance = 3.5,
                                    multithread = nproc)

# ---- ASV table (SeqID-keyed) ----
asv_seqs <- colnames(seqtab.nochim)
final_seqtab <- as.data.frame(t(seqtab.nochim))
rownames(final_seqtab) <- NULL
final_seqtab <- cbind(SeqID = asv_seqs, final_seqtab)
write.table(final_seqtab, file = "ASV_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- ASV fasta (>ASV_<i> headers) ----
fasta_lines <- unlist(lapply(seq_along(asv_seqs), function(i) {
    c(paste0(">ASV_", i), asv_seqs[i])
}))
writeLines(fasta_lines, "ASV_sequences.fasta")

# ---- track_control: SampleID, input, filtered, denoised, nonchim ----
getN <- function(x) sum(getUniques(x))
denoised <- sapply(dd, getN)
nonchim  <- rowSums(seqtab.nochim)
track <- data.frame(
    SampleID = sample.names,
    input    = out[sample.names, 1],
    filtered = out[sample.names, 2],
    denoised = denoised[sample.names],
    nonchim  = nonchim[sample.names],
    row.names = NULL,
    check.names = FALSE
)
write.table(track, file = "track_control.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- counts file (sample\tcount, one row per sample) ----
counts <- data.frame(sample = sample.names, count = nonchim[sample.names])
write.table(counts, file = "dada2.counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_pacbio.R done\n")
