# Gut_microbiota
A repo containing R scripts for the analysis of gut microbiome data from insect sequence data
```
library(dada2)
```

__loading the fastq files into R__
```
file_path <- "stingless_bee"
list.files(file_path)
```
      ## [1] "10K_DS40_L001_R1_001.fastq.gz"  "10K_DS40_L001_R2_001.fastq.gz" 
      ## [3] "11K_HA41_L001_R1_001.fastq.gz"  "11K_HA41_L001_R2_001.fastq.gz" 
      ## [5] "12K_HA42_L001_R1_001.fastq.gz"  "12K_HA42_L001_R2_001.fastq.gz" 
      ## [7] "13K_HA43_L001_R1_001.fastq.gz"  "13K_HA43_L001_R2_001.fastq.gz" 
      ## [9] "14K_HA44_L001_R1_001.fastq.gz"  "14K_HA44_L001_R2_001.fastq.gz" 
      ## [11] "15K_HA45_L001_R1_001.fastq.gz"  "15K_HA45_L001_R2_001.fastq.gz" 
__sorting the data__
```
dataF <- sort(list.files(file_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
dataR <- sort(list.files(file_path, pattern="_R2_001.fastq.gz", full.names = TRUE))
```
__Extracting the file names__
```
list.sample.names <- sapply(strsplit(basename(dataF), "_"), `[`, 1)
list.sample.names
```
__visualizing the quality of the plots__
```
plotQualityProfile(dataR[1:39])
plotQualityProfile(dataF[1:39])
```
__Assigning where the filtered data should be stored "filtered_directory__
```
filt.dataF <- file.path(file_path, "filtered", paste0(list.sample.names, "_F_filt.fastq.gz"))
filt.dataR <- file.path(file_path, "filtered", paste0(list.sample.names, "_R_filt.fastq.gz"))
names(filt.dataF) <- list.sample.names
names(filt.dataR) <- list.sample.names
```

__Filtering and trimming data__
```
out <- filterAndTrim(dataF, filt.dataF, dataR, filt.dataR, truncLen=c(290,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
```

__Establishing the error rates in the data for both the forwards and reverses__
```
errF <- learnErrors(filt.dataF, multithread=TRUE)
errR <- learnErrors(filt.dataR, multithread=TRUE)
```
__Plotting the error rates__
```
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```
__Denoising__
```
dadaF <- dada(filt.dataF, err=errF, multithread=TRUE)
dadaR <- dada(filt.dataR, err=errR, multithread=TRUE)
```
__inspecting the resulting dada2 class object__
```
dadaF[[1]]
dadaR[[1]]
```
__merging the forward and reverse reads__
```
merge.reads <- mergePairs(dadaF, filt.dataF, dadaR, filt.dataR, verbose=FALSE)
head(merge.reads[[1]])
```
__generating a sequence table__
```
seqtab <- makeSequenceTable(merge.reads)
dim(seqtab)
table(nchar(getSequences(seqtab))) #establishing the number of samples distributed per length
```
__Removing chimeras__
```
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #Establishing the percentage of non_chimeric reads
```
__Tracking the number of reads that passed the various steps in the pipeline__
```
getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(merge.reads, getN), rowSums(seqtab.nochim))
colnames(track.nbr.reads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track.nbr.reads) <- list.sample.names
head(track.nbr.reads)
```
__Taxonomic classification__
```
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
```
