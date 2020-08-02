# 16S Gut microbiota analysis pipeline using Dada2
This repo contains an R script that can be used in the analysis of 16S gut microbiome data. The scripts were developed by analysing the gut microbiome data from stingless bees. The 16S data was generated from an illumina platform by sequencing the V3-V4 region. The sequences used in developing this  script are paired end reads which have already been demultiplexed but still contain primers. The DADA2 workflow was adopted in this analysis because it is highly sensitive and specific as compared to other OTU picking algorithms, it can resolve single-nucleotide differences from amplicon data and classify them into Amplicon Sequence Variants (ASVs) and it contains an error correction model that helps in improving the quality of reads by helping best infer the original true bilogical sequences present in the data. More about Dada2 can be found in the Dada2 paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/, Dada2 tutorial https://benjjneb.github.io/dada2/tutorial.html and the Dada2 manual https://www.bioconductor.org/packages/3.3/bioc/manuals/dada2/man/dada2.pdf.

__Working around with the data__
Load the dada2 package into R. If you have not installed dada2, follow the turorial here http://benjjneb.github.io/dada2/dada-installation.html
```
library(dada2)
```

__Reading the fastq files into R__
```
file_path <- "stingless_bee/"
list.files(file_path)
```
      ## [1] "10K_DS40_L001_R1_001.fastq.gz"  "10K_DS40_L001_R2_001.fastq.gz" 
      ## [3] "11K_HA41_L001_R1_001.fastq.gz"  "11K_HA41_L001_R2_001.fastq.gz" 
      ## [5] "12K_HA42_L001_R1_001.fastq.gz"  "12K_HA42_L001_R2_001.fastq.gz" 
      ## [7] "13K_HA43_L001_R1_001.fastq.gz"  "13K_HA43_L001_R2_001.fastq.gz" 
      ## [9] "14K_HA44_L001_R1_001.fastq.gz"  "14K_HA44_L001_R2_001.fastq.gz" 
      ## [11] "15K_HA45_L001_R1_001.fastq.gz"  "15K_HA45_L001_R2_001.fastq.gz" 
If the list of files that have been read match those you intend to analyze and your dada2 package was loaded you can start your dada2 analysis

__sorting the data into forward and reverse reads__
```
dataF <- sort(list.files(file_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
dataR <- sort(list.files(file_path, pattern="_R2_001.fastq.gz", full.names = TRUE))
list.sample.names <- sapply(strsplit(basename(dataF), "_"), `[`, 2) #specifying the sample names for your dataset
list.sample.names
```
Here, you separate your forward and reverse reads and choose the string on the name of the file to be used as your sample names.

__visualizing the quality of the plots__
Visualizing the quality of your sequence data. Since most of the quality plots are closely similar in the two categories, you dont have to visualize the quality plots from all the samples in each of the categories. 
```
plotQualityProfile(dataR[1:5])
plotQualityProfile(dataF[1:5])
```
on this plot, the bases are on the X axis while the quality score is on the Y. Faint grey represents a heat map of the frequency of each quality score at each base position. Green shows the mean quity score per base and orange represents the quartiles of the quality score distribution. The forward reads generally usually have high quality bases while the reverse reads have more spurious reads as the sequencing process advances.  This plot also shows that our reads have primers as seen with the decreasing quality at the 5' end of each of the quality plots.

__Removing primers__
```
fwd_primer <- "CCTACGGGNGGCWGCAG"
rev_primer <- "GACTACHVGGGTATCTAATCC"
```
```
#Function for creating all the possible orientations of these primers
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
```
```
# creating all the possible orientations using our forward and reverse primers
fwd_primer_orients <- allOrients(fwd_primer)
rev_primer_orients <- allOrients(rev_primer)
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer))) # reverse complement of the primers
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
```
```
# Function for counting the number of reads containing our primer orientations 
count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}
```
```
# counting the sequence strings with primers from our data 
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]))
```
```
#specifying the path to cutadapt on the server
cutadapt <- "/opt/apps/cutadapt/1.18/bin/cutadapt" # change to the cutadapt path on your machine
system2(cutadapt, args = "--version") # test if the version of cutadapt you loaded matches that on the server
```
```
creating an output directory to store the trimmed files and specifying the names of the samples
cut_dir <- file.path(file_path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(dataF))
rev_cut <- file.path(cut_dir, basename(dataR))

names(fwd_cut) <- list.sample.names
names(rev_cut) <- list.sample.names
```
```
# creating a loop over the list of files and running cutadapt on each file.
for (i in seq_along(dataF)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   dataF[i], dataR[i]),
          stdout = cut_logs[i])  
}
```
```
# sanity check: checking for the presence of primers in the first cutadapt-ed sample.
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = rev_cut[[1]]))
```
__Assigning where the filtered data should be stored "filtered" directory__
```
filt.dataF <- file.path(file_path, "filtered", paste0(list.sample.names, "_F_filt.fastq.gz"))
filt.dataR <- file.path(file_path, "filtered", paste0(list.sample.names, "_R_filt.fastq.gz"))
names(filt.dataF) <- list.sample.names
names(filt.dataR) <- list.sample.names
```

__Filtering and trimming data__
```
out <- filterAndTrim(fwd_cut, filt.dataF, rev_cut, filt.dataR, truncLen=c(257,186),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
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
