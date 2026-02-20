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

if (length(args)==7){
    perform_denoising=TRUE
} else {
    perform_denoising=FALSE
    minAbundance <- args[8]
}

# Following this particular tutorial: https://github.com/benjjneb/LRASManuscript/blob/master/LRASms_fecal.Rmd
# Reads
reads_path <- list.files(input_path, pattern="fq", full.names=TRUE)
sample.names <- gsub("_trimmed_cleaned-oriented.fq","",basename(reads_path))

# Filtering
cat("Found ",length(reads_path)," samples. \n")
cat("Starting filtering step\n")

filts <- file.path(output_path, "filtered",paste0(sample.names,"_filt.fastq.gz"))
names(filts) <- sample.names
out <- filterAndTrim(reads_path, filts, minQ=minQ, minLen=minLen, maxLen=maxLen, maxN=maxN, rm.phix=FALSE, maxEE=maxEE,multithread=20,verbose=FALSE)


cat("Starting denoising step\n")
# Denoising
if (perform_denoising){
    err <- learnErrors(filts, errorEstimationFunction=dada2:::PacBioErrfun, BAND_SIZE=32, multithread=20)
    drp <- derepFastq(filts, verbose=FALSE)
    dd <- dada(drp, err=err, BAND_SIZE=32, multithread=20)    
} else {
    # only dereplicate
    dd <- derepFastq(filts, verbose=FALSE)
}

# ASV table
 seqtab <- makeSequenceTable(dd)

cat("Removing chimeras\n")
# Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, minFoldParentOverAbundance=3.5, multithread=TRUE) 
# *minFoldParentOverAbundance = Only sequences greater than this-fold more abundant than a sequence can be its "parents".

# Saving results
write.table(t(seqtab.nochim),file.path(output_path,"ASV_table.tsv"),row.names = TRUE,sep="\t",quote=FALSE)
cat("ASV_table saved\n")

## Extract ASV sequences (stored as column names)
asv_seqs <- colnames(seqtab.nochim)

# Create FASTA content
fasta_lines <- unlist(
  lapply(seq_along(asv_seqs), function(i) {
    c(paste0(">ASV_", i), asv_seqs[i])
  })
)

# Write to file
writeLines(fasta_lines, file.path(output_path, "ASV_sequences.fasta"))

# track
getN <- function(x) sum(getUniques(x))
if (length(sample.names)==1) track <- cbind(out, sum(getUniques(dd)), rowSums(seqtab.nochim))
else track <- cbind(out, sapply(dd, getN), rowSums(seqtab.nochim))
write.table(track,file.path(output_path,"track_control.tsv"),row.names = TRUE,sep="\t",quote=FALSE)

cat("track control saved\n")