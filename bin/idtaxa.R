#!/usr/bin/env Rscript
# IDTAXA classification for ASV fasta.
# Positional args: asv_fasta classifier nproc out_taxa_tsv out_taxa_conf_tsv
#   asv_fasta          ASV sequences (FASTA)
#   classifier         path to .RData with `trainingSet_custom`
#   nproc              parallel workers
#   out_taxa_tsv       SeqID + 7 rank columns
#   out_taxa_conf_tsv  same + minimum-rank confidence column

suppressMessages(suppressWarnings({
    library(Biostrings)
    library(DECIPHER)
}))

args               <- commandArgs(trailingOnly = TRUE)
fasta              <- args[1]
classifier         <- args[2]
nproc              <- as.integer(args[3])
out_taxa_tsv       <- args[4]
out_taxa_conf_tsv  <- args[5]

load(classifier)

dna <- readDNAStringSet(fasta)

tax_info <- IdTaxa(test = dna, trainingSet = trainingSet_custom,
                   strand = "both", processors = nproc)

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

max_len <- max(sapply(tax_info, function(x) length(x$taxon[-1])))

asv_tax_df <- as.data.frame(do.call(rbind, lapply(tax_info, function(x) {
    taxa <- x$taxon[-1]
    taxa[startsWith(taxa, "unclassified_")] <- "unassigned"
    taxa[startsWith(taxa, "uncultured")]    <- "unassigned"
    c(taxa, rep("unassigned", max_len - length(taxa)))
})), stringsAsFactors = FALSE)

colnames(asv_tax_df) <- ranks[seq_len(ncol(asv_tax_df))]
asv_tax_df <- cbind(SeqID = as.character(dna), asv_tax_df)
rownames(asv_tax_df) <- NULL

# Taxonomy column = prefix-joined string for downstream MetaStandard
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
rank_cols <- intersect(ranks, colnames(asv_tax_df))
asv_tax_df$Taxonomy <- apply(asv_tax_df[, rank_cols, drop = FALSE], 1, function(row) {
    keep <- row[row != "unassigned"]
    if (length(keep) == 0) return("")
    paste(paste0(prefixes[seq_along(keep)], keep), collapse = ";")
})

confidence_df <- sapply(tax_info, function(x) min(x$confidence))

asv_tax_conf_df <- asv_tax_df
asv_tax_conf_df$Confidence <- confidence_df

write.table(asv_tax_df,      file = out_taxa_tsv,
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(asv_tax_conf_df, file = out_taxa_conf_tsv,
            sep = "\t", row.names = FALSE, quote = FALSE)
