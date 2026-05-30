# HAWAII MICROBIOME PROJECT
# Sediment samples from Derek Esibill
# Summer 2024
# DADA2 Analysis
# Last updated: 12/13/24 by ANB

# PLEASE READ: This sequence run was on 11/27/24 conducted on a MiSeq at UMass Boston (P. Kearns).
# The reverse reads are low quality, so this pipeline analyzes the forward reads ONLY.
# The plan is to re-sequence at a future date, or to shuttle the money towards metaG.

# Load libraries-----------------
library(dada2)
library(ShortRead)
library(patchwork)

# Set working directory-----------
path="~/Dropbox/HAWAII-MICROBIOME-PROJ/collab-hawaii-microbiome/RAW-SEQUENCES"
list.files(path)

# Inspect reads----
# Need to make sure this matches what your sequencs look like

# Add the .gz if your files are still gunzipped
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)

# Inspect read quality profiles
P1 <- plotQualityProfile(fnFs[1:5]) # change to your number of samples, max
P2 <- plotQualityProfile(fnRs[1:5])
P1 + P2 # This doesn't look good. Will do a quick analysis but then may need to return here later

# Filter & Trim----
# Moving forward, will only be analyzing the forward reads
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trimming primers
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-30

# Here, you are trimming your reads for quality
# Removing primers
# Truncating the ends of the reads
out <- filterAndTrim(fnFs, filtFs,
                     trimLeft = FWD_PRIMER_LEN,
                     truncLen = 150,  # Example values
                     maxN = 0,  # Discard reads with any Ns
                     maxEE = 2,  # Allow max 2 expected errors
                     truncQ = 2,  # Truncate reads at the first quality score <= 2
                     compress = TRUE,  # Compress the output files
                     verbose = TRUE, 
                     multithread = TRUE) 
head(out)
saveRDS(out,"out.rds")

# Learn the error rates----
# DADA 2's algorithm is based on modeling the error rates in the sequencing process. 
# It leverages a parametric error model (err), which is different for each individual 
# amplicon dataset. The learnErrors method learns this error model from the data, by
# alternating estimation of the error rates and inference of sample composition until 
# they converge on a jointly consistent solution.
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF,"errF.rds")

# Visualize the errors
plotErrors(errF, nominalQ=TRUE)

# Dereplicate----
derepFs <- derepFastq(filtFs, verbose=TRUE)
saveRDS(derepFs, "derepFS.rds")
names(derepFs) <- sample.names

# To read back in
# derepFs <- readRDS("derepFs.rds")
# derepRs <- readRDS("derepRs.rds")
# Name the derep-class objects by the sample names

# Run the DADA2 algorithm----
# If you don't replicate, use "filtFs" 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "dadaFs.rds")
dadaFs[[1]]

# To read back in
# dadaFs <- readRDS("dadaFs.rds")

# Construct a sequence table----
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras----
# Chimeras are common but shouldn't be more than 3-4% 
# If high, then remove primers and redo analysis
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # This means that ~2.5% were chimeric
saveRDS(seqtab.nochim, "seqtab-nochim.rds")

# Do a quick check in----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

# Assigning taxonomy----
# Be sure to put the database in your working directory or another known location
# You can download proper database here: https://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, "taxa-table.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Preparing for Phyloseq----------------
# Constructing a sample table that we can use for Phyloseq import
seqs <- getSequences(seqtab.nochim)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

# Sample names:
write.table(sample.names, "sample_names.txt", sep="\t", quote=F, col.names=NA)



