suppressMessages(suppressWarnings({
library("ShortRead")
library("dada2")
library("Biostrings")
library(reshape2)
library(DECIPHER)
}))

# Get command-line arguments (excluding default R arguments)
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_path  <- args[1]
output_path <- args[2]
minQ <- args[3]
minLen <- args[4]
maxLen <- args[5]
maxN <- args[6]
maxEE <- args[7]
decipher_classifier <- args[8]

# Following this particular tutorial: https://github.com/benjjneb/LRASManuscript/blob/master/LRASms_fecal.Rmd
# Reads
reads_path <- list.files(input_path, pattern="fastq.gz", full.names=TRUE)

# Filtering
filts <- file.path(reads_path, "filtered", basename(reads_path))
out <- filterAndTrim(reads_path, filts, minQ=minQ, minLen=minLen, maxLen=maxLen, maxN=maxN, rm.phix=FALSE, maxEE=maxEE,multithread=20)

# Denoising
err <- learnErrors(filts, errorEstimationFunction=dada2:::PacBioErrfun, BAND_SIZE=32, multithread=20)
drp <- derepFastq(filts, verbose=FALSE)
dd <- dada(drp, err=err, BAND_SIZE=32, multithread=20) 

# ASV table
seqtab <- makeSequenceTable(dd)

# Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, minFoldParentOverAbundance=3.5, multithread=TRUE) 

# Taxonomy
load(decipher_classifier)
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet_custom, strand="both", processors=20)
 
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Find the maximum taxonomy depth
max_len <- max(sapply(tax_info, function(x) length(x$taxon[-1])))

# Process and pad each taxonomy vector
asv_tax_df <- as.data.frame(do.call(rbind, lapply(tax_info, function(x) {
    taxa <- x$taxon[-1]
    taxa[startsWith(taxa, "unclassified_")] <- "unassigned"
    taxa[startsWith(taxa, "uncultured")] <- "unassigned"
    # pad with NA if shorter
    c(taxa, rep("unassigned", max_len - length(taxa)))
})), stringsAsFactors = FALSE)

# assign ASV IDs as row names
colnames(asv_tax_df) <- ranks
rownames(asv_tax_df) <- as.character(dna)
asv_tax_df <- tibble::rownames_to_column(asv_tax_df, "SeqID")

confidence_df <- sapply(tax_info, function(x) {
    conf <- x$confidence
    conf <- min(conf)
    conf
})

asv_tax_conf_df <- asv_tax_df
asv_tax_conf_df$confidence <- confidence_df


# track version
cat(paste("DECIPHER",packageVersion("DECIPHER")), "\n", 
    file = file.path(path_project_dir,"run_info/tools.txt"), append = TRUE)

cat(paste("Classifier", basename(decipher_classifier)), "\n", 
    file = file.path(path_project_dir,"run_info/tools.txt"), append = TRUE)

# track
track <- cbind(out, sapply(dd, getN), rowSums(seqtab.nochim))

# saving results
write.table(track,file.path(output_path,"track_control.tsv"),row.names = TRUE,sep="\t",quote=FALSE)
write.table(seqtab.nochim,file.path(output_path,"ASV_table.tsv"),row.names = TRUE,sep="\t",quote=FALSE)
write.table(asv_tax_df,file=file.path(output_path,"taxa_table.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(asv_tax_conf_df,file=file.path(output_path,"taxa_table_conf.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
