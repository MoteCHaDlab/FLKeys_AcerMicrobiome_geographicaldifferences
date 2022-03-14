library(dada2); packageVersion("dada2")

setwd("~/acer_microbiome_2020")

# Location of raw reads
path <- "/Users/klingesj/acer_microbiome_2020/demultiplexed"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Plot read quality- you expect read quality to drop off towards the end
#plotQualityProfile(fnFs[1:4])
#plotQualityProfile(fnRs[1:4])
#In gray-scale is a heat map of the frequency of each quality score at each 
#base position. The mean quality score at each position is shown by the green 
#line, and the quartiles of the quality score distribution by the orange lines.
#The red line shows the scaled proportion of reads that extend to at least that
#position (this is more useful for other sequencing technologies, as Illumina 
#reads are typically all the same length, hence the flat red line).

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1)
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

#write.table(sample.names, file = "names.txt", sep="\t")

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(210,210),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)

write.table(out, file = "filter_out.txt", sep="\t")

filtpath <- "/Users/klingesj/acer_microbiome_2020/demultiplexed/filtered"
filtFs <- list.files(filtpath, pattern="F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpath, pattern="R_filt.fastq.gz", full.names = TRUE)
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# If you like, visualize these error rates. You want the estimated error rates (black lines) to be a good 
# fit to the observed error rates (points) and for the error rates drop with increased quality
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

# Note: if this loop doesn't work, you've done something wrong. You've probably made an error earlier in the code, 
# probably with the location of your forward and reverse reads and it's generated duplicate sample names
# Normal to get warning message that there are duplicate sequences
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE) # core sample interference algorithm
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE) # core sample interference algorithm
  merger <- mergePairs(ddF, derepF, ddR, derepR) # Merge paired end reads
  mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)
# Make sequence table from merged reads
st.all <- makeSequenceTable(mergers) # Normal to get warning message saying the sequences being tabled vary in length

# Inspect distribution of read lengths
table(nchar(getSequences(st.all)))
hist(nchar(getSequences(st.all)))

# Remove any ASVs that are considerably off target length
seqtab_trimmed <- st.all[,nchar(colnames(st.all)) %in% seq(250,255)]

# Inspect distribution of read lengths after removal of off-target reads
table(nchar(getSequences(seqtab_trimmed)))
hist(nchar(getSequences(seqtab_trimmed)))

# Remove chimeric sequences
seqtab <- removeBimeraDenovo(seqtab_trimmed, method="consensus", multithread=TRUE, verbose = T) #Identified 776 bimeras out of 2892 input sequences
sum(st.all)-sum(seqtab) # How many chimeras were removed? 696841
sum(seqtab)/sum(st.all) #0.9660389 

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2],
                          nonchim=rowSums(seqtab),
                          final_perc_reads_retained=round(rowSums(seqtab)/out[,1]*100, 1))

summary_tab 
write.table(summary_tab, file = "reads_lost.txt", sep="\t")

# Save chimera-free ASV table as downstream tasks may cause R to crash
saveRDS(seqtab, "seqtab.rds")

# Assign taxonomy based on silva reference database at genus level, you must have the appropriate Silva database downloaded
tax_silva <- assignTaxonomy(seqtab, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Assign taxonomy based on silva reference database at species (100%) level
silva_sp <- addSpecies(tax_silva, "silva_species_assignment_v132.fa.gz")

library(phyloseq)
# Export sequence table with genus and species assignments as phyloseq objects
ps_object <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(tax_silva))
ps_object_sp <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), tax_table(silva_sp))

# Save as RDS objects
saveRDS(ps_object, file = "ps_object.rds")
saveRDS(ps_object_sp, file = "ps_object_sp.rds")
