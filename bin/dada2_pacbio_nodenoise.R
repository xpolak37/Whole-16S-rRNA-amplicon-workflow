#!/usr/bin/env Rscript
# DADA2 no-denoise path for PacBio HiFi 16S (single-end, cohort-level).
# Skips learnErrors/dada — dereplicate, drop low-abundance ASVs, remove chimeras.
#
# Args (sibling-style long flags):
#   --input <dir>           directory containing *-oriented.fq files
#   --nproc <n>             worker threads
#   --minQ <int>            filterAndTrim minQ
#   --minLen <int>          filterAndTrim minLen
#   --maxLen <int>          filterAndTrim maxLen
#   --maxN <int>            filterAndTrim maxN
#   --maxEE <num|Inf>       filterAndTrim maxEE
#   --minAbundance <num>    drop ASVs whose total proportion of reads is below
#                           this threshold (default 0.001 = 0.1%)
#
# Outputs (cwd):
#   ASV_table.tsv               SeqID (ASV_<i>) + Sequence + per-sample abundance cols
#   ASV_sequences.fasta         >ASV_<i> headers (same order as ASV_table.tsv)
#   track_control.tsv           SampleID, input, filtered, denoised, nonchim
#                               (denoised here = unique-after-derep, per-sample)
#   dada2_nodenoise.counts.tsv  sample\tcount  (count = nonchimeric reads)

suppressMessages(suppressWarnings({
    library(dada2)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) > 0) args[idx + 1] else default
}

input_dir    <- get_arg(args, "--input")
nproc        <- as.integer(get_arg(args, "--nproc", "1"))
minQ         <- as.integer(get_arg(args, "--minQ",   "0"))
minLen       <- as.integer(get_arg(args, "--minLen", "1000"))
maxLen       <- as.integer(get_arg(args, "--maxLen", "1600"))
maxN         <- as.integer(get_arg(args, "--maxN",   "0"))
maxEE_raw    <- get_arg(args, "--maxEE", "Inf")
maxEE        <- if (identical(maxEE_raw, "Inf")) Inf else as.numeric(maxEE_raw)
minAbundance <- as.numeric(get_arg(args, "--minAbundance", "0.001"))

if (is.null(input_dir)) stop("dada2_pacbio_nodenoise.R: --input is required")

reads_path   <- list.files(input_dir, pattern = "-oriented\\.fq$", full.names = TRUE)
if (length(reads_path) == 0) stop("dada2_pacbio_nodenoise.R: no *-oriented.fq files in ", input_dir)
sample.names <- sub("-oriented\\.fq$", "", basename(reads_path))

cat("Found", length(reads_path), "samples (no-denoise mode)\n")

# ---- filter ----
filts <- file.path("filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filts) <- sample.names
out <- filterAndTrim(reads_path, filts,
                     minQ = minQ, minLen = minLen, maxLen = maxLen,
                     maxN = maxN, maxEE = maxEE, rm.phix = FALSE,
                     multithread = nproc, verbose = FALSE)
rownames(out) <- sample.names

# ---- dereplicate only ----
drp <- derepFastq(filts, verbose = FALSE)
if (inherits(drp, "derep")) drp <- setNames(list(drp), sample.names)

# ---- ASV table + low-abundance filter + chimeras ----
seqtab <- makeSequenceTable(drp)

total_reads   <- sum(seqtab)
asv_total     <- colSums(seqtab)
asv_proportion <- asv_total / total_reads
keep <- asv_proportion >= minAbundance
cat("Dropping", sum(!keep), "/", length(keep),
    "ASVs below minAbundance =", minAbundance, "\n")
seqtab <- seqtab[, keep, drop = FALSE]

seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    minFoldParentOverAbundance = 3.5,
                                    multithread = nproc)

# ---- ASV table (ASV_<i>-keyed, sequence preserved in Sequence column) ----
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV_", seq_along(asv_seqs))
final_seqtab <- as.data.frame(t(seqtab.nochim))
rownames(final_seqtab) <- NULL
final_seqtab <- cbind(SeqID = asv_ids, Sequence = asv_seqs, final_seqtab)
write.table(final_seqtab, file = "ASV_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- ASV fasta (>ASV_<i> headers, same order as ASV_table.tsv) ----
fasta_lines <- unlist(lapply(seq_along(asv_seqs), function(i) {
    c(paste0(">", asv_ids[i]), asv_seqs[i])
}))
writeLines(fasta_lines, "ASV_sequences.fasta")

# ---- track_control: per-sample (fixes the bash whole-batch shape bug) ----
getN <- function(x) sum(getUniques(x))
denoised <- sapply(drp, getN)
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

# ---- counts file ----
counts <- data.frame(sample = sample.names, count = nonchim[sample.names])
write.table(counts, file = "dada2_nodenoise.counts.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_pacbio_nodenoise.R done\n")
