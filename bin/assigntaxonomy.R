suppressMessages(suppressWarnings({
    library(Biostrings)
    library(dada2)
}))

args           <- commandArgs(trailingOnly = TRUE)
fasta          <- args[1]
classifier     <- args[2]
nproc          <- as.integer(args[3])
denoising_tool <- args[4]

# load fasta
dna <- readDNAStringSet(fasta)

# assign taxonomy
tax_info <- assignTaxonomy(
    seqs        = as.character(dna),
    refFasta    = classifier,
    multithread = nproc,
    outputBootstraps = TRUE,
    verbose     = TRUE
)

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# extract taxonomy
asv_tax_df <- as.data.frame(tax_info$tax, stringsAsFactors = FALSE)
asv_tax_df[is.na(asv_tax_df)] <- "unassigned"

# pad columns to match ranks if species not present
while (ncol(asv_tax_df) < length(ranks)) {
    asv_tax_df[[ranks[ncol(asv_tax_df) + 1]]] <- "unassigned"
}
colnames(asv_tax_df) <- ranks[seq_len(ncol(asv_tax_df))]

# Collapse taxonomy ranks into a single semicolon-separated string
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")

# Collapse taxonomy ranks into a single semicolon-separated string with prefixes
taxonomy_string <- apply(asv_tax_df, 1, function(x) {
    x <- x[x != "unassigned"]
    paste(paste0(prefixes[1:length(x)], x), collapse = ";")
})

# extract bootstrap confidence - take minimum across ranks
confidence_df <- apply(tax_info$boot, 1, function(x) min(x, na.rm = TRUE))

# Build final dataframe
asv_tax_conf_df <- data.frame(
    SeqID      = if (denoising_tool == "unoise") names(dna) else as.character(dna),
    Taxonomy   = taxonomy_string,
    Confidence = confidence_df,
    stringsAsFactors = FALSE
)

write.table(asv_tax_conf_df,
            file = paste0(denoising_tool, "_taxa_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)