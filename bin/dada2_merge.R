#!/usr/bin/env Rscript
# DADA2 sequence-table assembly + chimera removal for PacBio HiFi 16S.
#
# Split-mode stage 4 of 4 (filter -> learn_errors -> denoise -> merge).
#
# Chimera removal is a cohort-level decision (removeBimeraDenovo's consensus
# method votes across samples), so it stays a single task. It is cheap on
# memory here because it operates on denoised ASVs -- roughly 700 per sample
# rather than the ~80k raw uniques a derep object carries.
#
# The table-construction calls below deliberately mirror dada2_pacbio.R
# line-for-line so that split mode and cohort mode emit identical files.
#
# Args (sibling-style long flags):
#   --uniques <dir>         directory of <sample>.uniques.tsv
#   --filter_stats <dir>    directory of <sample>.filter_stats.tsv
#   --denoise_stats <dir>   directory of <sample>.denoise_stats.tsv
#   --nproc <n>             worker threads
#   --min_fold <num>        removeBimeraDenovo minFoldParentOverAbundance
#   --counts_name <file>    name of the counts output file
#
# Outputs (cwd):
#   ASV_table.tsv          SeqID (ASV_<i>) + Sequence + per-sample abundance cols
#   ASV_sequences.fasta    >ASV_<i> headers (in the same order as ASV_table.tsv)
#   track_control.tsv      SampleID, input, filtered, denoised, nonchim
#   dada2.counts.tsv       sample\tcount  (count = nonchimeric reads)

suppressMessages(suppressWarnings({
    library(dada2)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) > 0) args[idx + 1] else default
}

uniques_dir  <- get_arg(args, "--uniques")
filt_dir     <- get_arg(args, "--filter_stats")
den_dir      <- get_arg(args, "--denoise_stats")
nproc        <- as.integer(get_arg(args, "--nproc", "1"))
min_fold     <- as.numeric(get_arg(args, "--min_fold", "3.5"))
counts_name  <- get_arg(args, "--counts_name", "dada2.counts.tsv")

if (is.null(uniques_dir)) stop("dada2_merge.R: --uniques is required")
if (is.null(filt_dir))    stop("dada2_merge.R: --filter_stats is required")
if (is.null(den_dir))     stop("dada2_merge.R: --denoise_stats is required")

read_stats <- function(dir, pattern) {
    files <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop("dada2_merge.R: no ", pattern, " files in ", dir)
    do.call(rbind, lapply(files, read.delim))
}

# The filter stats are the authoritative sample list: every sample reaches this
# stage even if it lost all of its reads during filtering.
filter_stats <- read_stats(filt_dir, "\\.filter_stats\\.tsv$")

# Deterministic, locale-independent ordering shared with dada2_learn_errors.R.
sample.names <- sort(filter_stats$SampleID, method = "radix")
filter_stats <- filter_stats[match(sample.names, filter_stats$SampleID), ]

cat("Merging", length(sample.names), "samples\n")

# ---- load per-sample denoised uniques ----
uniques_for <- function(s) {
    f <- file.path(uniques_dir, paste0(s, ".uniques.tsv"))
    if (!file.exists(f)) return(NULL)
    tab <- read.delim(f, colClasses = c("character", "integer"))
    setNames(tab$Abundance, tab$Sequence)
}
drp_list <- lapply(sample.names, uniques_for)
names(drp_list) <- sample.names

empty <- vapply(drp_list, is.null, logical(1))
if (all(empty)) stop("dada2_merge.R: every sample is empty after denoising")
if (any(empty)) {
    cat("Samples with no surviving reads (backfilled as zero):",
        paste(sample.names[empty], collapse = ", "), "\n")
}

# ---- ASV table + chimeras ----
# Empty samples are held out of removeBimeraDenovo so they cannot perturb the
# cross-sample consensus vote, then re-inserted as zero columns below.
seqtab        <- makeSequenceTable(drp_list[!empty])
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    minFoldParentOverAbundance = min_fold,
                                    multithread = nproc)

# ---- ASV table (ASV_<i>-keyed, sequence preserved in Sequence column) ----
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV_", seq_along(asv_seqs))
final_seqtab <- as.data.frame(t(seqtab.nochim))
rownames(final_seqtab) <- NULL

if (any(empty)) {
    # Re-insert dropped samples as zero columns, restoring full sample order.
    zero <- as.data.frame(matrix(0L, nrow = length(asv_seqs), ncol = sum(empty)))
    colnames(zero) <- make.names(sample.names[empty])
    final_seqtab <- cbind(final_seqtab, zero)
    final_seqtab <- final_seqtab[, make.names(sample.names), drop = FALSE]
}

final_seqtab <- cbind(SeqID = asv_ids, Sequence = asv_seqs, final_seqtab)
write.table(final_seqtab, file = "ASV_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- ASV fasta (>ASV_<i> headers, same order as ASV_table.tsv) ----
fasta_lines <- unlist(lapply(seq_along(asv_seqs), function(i) {
    c(paste0(">", asv_ids[i]), asv_seqs[i])
}))
writeLines(fasta_lines, "ASV_sequences.fasta")

# ---- track_control: SampleID, input, filtered, denoised, nonchim ----
denoise_stats <- do.call(rbind, lapply(
    list.files(den_dir, pattern = "\\.denoise_stats\\.tsv$", full.names = TRUE),
    read.delim))

denoised <- setNames(rep(0L, length(sample.names)), sample.names)
denoised[denoise_stats$SampleID] <- as.integer(denoise_stats$denoised)

nonchim <- setNames(rep(0L, length(sample.names)), sample.names)
row_sums <- rowSums(seqtab.nochim)
nonchim[names(row_sums)] <- as.integer(row_sums)

track <- data.frame(
    SampleID = sample.names,
    input    = filter_stats[, "input"],
    filtered = filter_stats[, "filtered"],
    denoised = denoised[sample.names],
    nonchim  = nonchim[sample.names],
    row.names = NULL,
    check.names = FALSE
)
write.table(track, file = "track_control.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- counts file (sample\tcount, one row per sample) ----
counts <- data.frame(sample = sample.names, count = nonchim[sample.names])
write.table(counts, file = counts_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_merge.R done -", length(asv_seqs), "ASVs across",
    length(sample.names), "samples\n")
