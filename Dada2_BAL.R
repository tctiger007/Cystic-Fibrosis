################ 16S DADA2 tutorial
################ Please refer to the online tutorial
################ https://benjjneb.github.io/dada2/tutorial.html

library(stringr)
library(dada2)
library(rstudioapi)
library(DECIPHER)
path = dirname(getActiveDocumentContext()$path)
path = paste0(str_remove(path, "BAL_species/Code"), "Early CF Lab Computer/BAL/Run287")
 # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#plotQualityProfile(fnFs[5:6])
#plotQualityProfile(fnRs[5:6])

################# Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) 
# On Windows set multithread=FALSE
head(out)

################ Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)

################ Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# Inspecting the returned dada-class object
dadaFs[[1]]

################ Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


############### Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# saveRDS(seqtab, "seqtab.")


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
############### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

############### Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


############### Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, 
                       "../Data/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE)
taxa <- addSpecies(taxa, "../Data/silva_species_assignment_v138.1.fa.gz")

# inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa, "../Data/tax_BAL_decipher.csv")
seqtab.t <- t(seqtab.nochim)
write.csv(seqtab.t, "../Data/otu_BAL_decipher.csv")





