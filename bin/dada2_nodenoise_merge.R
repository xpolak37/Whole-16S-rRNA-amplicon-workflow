#!/usr/bin/env Rscript
# Sequence-table assembly, low-abundance filtering and chimera removal for the
# DADA2 no-denoise path.
#
# Split-mode stage 3 of 3 for dada2_nodenoise (filter -> derep -> merge).
#
# Why this is not just a copy of dada2_merge.R
# -------------------------------------------
# The denoised path merges dada() output: ~700 ASVs per sample, so the merged
# table is small. This path merges RAW uniques, which is a different problem --
# a cohort of 400 PacBio samples carries tens of millions of distinct sequences,
# and makeSequenceTable builds a DENSE samples x sequences matrix. Splitting
# derep per sample alone would therefore not make this path scale; the merge
# stage would simply become the new bottleneck.
#
# The fix is to discard, before building the table, sequences that provably
# cannot survive the minAbundance filter:
#
#   The cohort script keeps a sequence when total(s) >= minAbundance * total_reads.
#   Call that right-hand side T. If a sequence's count is below t in EVERY sample
#   then its total is strictly below t * n_samples. Choosing
#
#       t = floor(T / n_samples)   =>   t * n_samples <= T
#
#   means any such sequence has total < T and would have been dropped anyway.
#
# So a sequence is a candidate iff it reaches t in at least one sample. That test
# needs one pass; a second pass is then required to total each candidate across
# ALL samples, including the samples where it sits below t -- those reads still
# count toward its total. Two passes over gzipped per-sample files, holding only
# candidates in memory, gives exactly the cohort script's answer.
#
# Args (sibling-style long flags):
#   --uniques <dir>         directory of <sample>.uniques.tsv.gz
#   --filter_stats <dir>    directory of <sample>.filter_stats.tsv
#   --derep_stats <dir>     directory of <sample>.derep_stats.tsv
#   --nproc <n>             worker threads
#   --min_fold <num>        removeBimeraDenovo minFoldParentOverAbundance
#   --minAbundance <num>    drop ASVs below this proportion of total reads
#   --counts_name <file>    name of the counts output file
#
# Outputs (cwd):
#   ASV_table.tsv, ASV_sequences.fasta, track_control.tsv, <counts_name>

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
derep_dir    <- get_arg(args, "--derep_stats")
nproc        <- as.integer(get_arg(args, "--nproc", "1"))
min_fold     <- as.numeric(get_arg(args, "--min_fold", "3.5"))
minAbundance <- as.numeric(get_arg(args, "--minAbundance", "0.001"))
counts_name  <- get_arg(args, "--counts_name", "dada2_nodenoise.counts.tsv")

if (is.null(uniques_dir)) stop("dada2_nodenoise_merge.R: --uniques is required")
if (is.null(filt_dir))    stop("dada2_nodenoise_merge.R: --filter_stats is required")
if (is.null(derep_dir))   stop("dada2_nodenoise_merge.R: --derep_stats is required")

read_stats <- function(dir, pattern) {
    files <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) stop("dada2_nodenoise_merge.R: no ", pattern, " in ", dir)
    do.call(rbind, lapply(files, read.delim))
}

filter_stats <- read_stats(filt_dir,  "\\.filter_stats\\.tsv$")
derep_stats  <- read_stats(derep_dir, "\\.derep_stats\\.tsv$")

# Deterministic, locale-independent ordering shared with the other split stages.
sample.names <- sort(filter_stats$SampleID, method = "radix")
filter_stats <- filter_stats[match(sample.names, filter_stats$SampleID), ]

uniq_path <- function(s) file.path(uniques_dir, paste0(s, ".uniques.tsv.gz"))
present   <- vapply(sample.names, function(s) file.exists(uniq_path(s)), logical(1))
if (!any(present)) stop("dada2_nodenoise_merge.R: no per-sample uniques found")

read_uniques <- function(s) {
    read.delim(gzfile(uniq_path(s)), colClasses = c("character", "integer"))
}

# total_reads must match the cohort script's sum(seqtab), which spans every raw
# unique -- not just the candidates kept below.
total_reads <- sum(as.numeric(derep_stats$reads))
n_present   <- sum(present)
threshold   <- minAbundance * total_reads
t_cut       <- floor(threshold / n_present)

cat("Merging", length(sample.names), "samples (", n_present, "with reads )\n")
cat("total reads:", format(total_reads, scientific = FALSE),
    " minAbundance:", minAbundance,
    " global threshold:", format(threshold, scientific = FALSE), "\n")
cat("safe per-sample prefilter cutoff:", t_cut, "\n")
if (t_cut <= 1) {
    cat("NOTE: cutoff <= 1, so no sequences can be prefiltered; the merge will\n",
        "     hold every raw unique. Expect high memory on large cohorts.\n")
}

# ---- pass 1: candidate sequences ----
candidates <- character(0)
for (s in sample.names[present]) {
    tab <- read_uniques(s)
    candidates <- union(candidates, tab$Sequence[tab$Abundance >= t_cut])
}
if (length(candidates) == 0) {
    stop("dada2_nodenoise_merge.R: no sequence reaches the minAbundance threshold")
}
cat("pass 1: ", length(candidates), " candidate sequences\n")

# ---- pass 2: exact totals for candidates across every sample ----
# Restricted to candidates, then handed to makeSequenceTable so that dada2's own
# abundance ordering applies exactly as it would in the cohort script.
drp_list <- lapply(sample.names[present], function(s) {
    tab <- read_uniques(s)
    hit <- tab$Sequence %in% candidates
    setNames(tab$Abundance[hit], tab$Sequence[hit])
})
names(drp_list) <- sample.names[present]

seqtab <- makeSequenceTable(drp_list)
cat("pass 2: sequence table", nrow(seqtab), "x", ncol(seqtab), "\n")

# ---- low-abundance filter (identical arithmetic to the cohort script) ----
asv_total      <- colSums(seqtab)
asv_proportion <- asv_total / total_reads
keep           <- asv_proportion >= minAbundance
cat("Dropping", sum(!keep), "/", length(keep),
    "candidate ASVs below minAbundance =", minAbundance, "\n")
seqtab <- seqtab[, keep, drop = FALSE]

seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    minFoldParentOverAbundance = min_fold,
                                    multithread = nproc)

# ---- ASV table ----
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV_", seq_along(asv_seqs))
final_seqtab <- as.data.frame(t(seqtab.nochim))
rownames(final_seqtab) <- NULL

if (!all(present)) {
    zero <- as.data.frame(matrix(0L, nrow = length(asv_seqs), ncol = sum(!present)))
    colnames(zero) <- make.names(sample.names[!present])
    final_seqtab <- cbind(final_seqtab, zero)
    final_seqtab <- final_seqtab[, make.names(sample.names), drop = FALSE]
}

final_seqtab <- cbind(SeqID = asv_ids, Sequence = asv_seqs, final_seqtab)
write.table(final_seqtab, file = "ASV_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

fasta_lines <- unlist(lapply(seq_along(asv_seqs), function(i) {
    c(paste0(">", asv_ids[i]), asv_seqs[i])
}))
writeLines(fasta_lines, "ASV_sequences.fasta")

# ---- track_control ----
denoised <- setNames(rep(0L, length(sample.names)), sample.names)
denoised[derep_stats$SampleID] <- as.integer(derep_stats$reads)

nonchim  <- setNames(rep(0L, length(sample.names)), sample.names)
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

counts <- data.frame(sample = sample.names, count = nonchim[sample.names])
write.table(counts, file = counts_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("dada2_nodenoise_merge.R done -", length(asv_seqs), "ASVs across",
    length(sample.names), "samples\n")
