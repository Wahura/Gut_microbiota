library(dada2)

#loading the fastq files into R
file_path <- "test_data"
head(list.files(file_path))

#sorting the data
dataF <- sort(list.files(file_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
head(dataF)
dataR <- sort(list.files(file_path, pattern="_R2_001.fastq.gz", full.names = TRUE))
head(dataR)
#Extracting the file names
list.sample.names <- sapply(strsplit(basename(dataF), "_"), `[`, 2)
list.sample.names

#visualizing the quality of the plots
plotQualityProfile(dataF[1:5])
plotQualityProfile(dataR[1:5])

library(ShortRead)
library(Biostrings)
library(stringr)

#primer trimming
fwd_primer <- "CCTACGGGNGGCWGCAG"
rev_primer <- "GACTACHVGGGTATCTAATCC"

# Function to check the orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# get the posisible orientations of the forward and reverse primers
fwd_primer_orients <- allOrients(fwd_primer)
fwd_primer_orients
rev_primer_orients <- allOrients(rev_primer)
rev_primer_orients
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer))) # reverse complement of the primers
fwd_primer_rev
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
rev_primer_rev

#function for counting the reads containing primers
count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

# counting the primers on one set of paired end FASTQ files
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]))

rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]))

# path to cutadapt
cutadapt <- "/opt/apps/cutadapt/1.18/bin/cutadapt" # change to the cutadapt path on your machine
system2(cutadapt, args = "--version") # confirm if the version matches that on the server

# Create an output directory called cutadapt to store the clipped files
cut_dir <- file.path(file_path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(dataF))
rev_cut <- file.path(cut_dir, basename(dataR))

names(fwd_cut) <- list.sample.names
head(fwd_cut)
names(rev_cut) <- list.sample.names
head(rev_cut)

# function for creating cutadapt trimming log files
cut_logs <- path.expand(file.path(cut_dir, paste0(list.sample.names, ".log")))

# Function specifying the cutadapt functions to be used in this analysis
cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
                   "-G", rev_primer, "-A", fwd_primer_rev,
                   "-n", 2,"-m",1, "-j",32, "--discard-untrimmed")

# creating a loop over the list of files and running cutadapt on each file.
for (i in seq_along(dataF)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   dataF[i], dataR[i]),
          stdout = NULL)  
}

#checking if the forward and reverse primers have been trimmed
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]))

rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = rev_cut[[1]]))

#plotting the quality plots to determine the filtering lengths
plotQualityProfile(fwd_cut[1:5])
plotQualityProfile(rev_cut[1:5])

#assigning where the filtered data should be stored "filtered_directory_
filt.dataF <- file.path(file_path, "filtered", paste0(list.sample.names, "_F_filt.fastq.gz"))
filt.dataR <- file.path(file_path, "filtered", paste0(list.sample.names, "_R_filt.fastq.gz"))
names(filt.dataF) <- list.sample.names
head(filt.dataF)
names(filt.dataR) <- list.sample.names
head(filt.dataR)

#filtering and trimming data
out <- filterAndTrim(fwd_cut, filt.dataF, rev_cut, filt.dataR, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#checking the quality of the quality trimmed reads
plotQualityProfile(filt.dataF[1:5])
plotQualityProfile(filt.dataR[1:5])

#Establishing the error rates in the data for both the forwards and reverses
errF <- learnErrors(filt.dataF, multithread=TRUE)
errR <- learnErrors(filt.dataR, multithread=TRUE)

#plotting the error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#denoising
dadaF <- dada(filt.dataF, err=errF, multithread=TRUE)
dadaR <- dada(filt.dataR, err=errR, multithread=TRUE)

dadaF[[1]]
dadaR[[1]]

#merging the forward and reverse reads
merge.reads <- mergePairs(dadaF, filt.dataF, dadaR, filt.dataR, verbose=FALSE)
head(merge.reads[[1]])

seqtab <- makeSequenceTable(merge.reads)
dim(seqtab)

#inspecting the number of samples distributed per length
table(nchar(getSequences(seqtab)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(merge.reads, getN), rowSums(seqtab.nochim),
                         final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))

colnames(track.nbr.reads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "per_retained")
rownames(track.nbr.reads) <- list.sample.names
head(track.nbr.reads)

#Taxonomic classification
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz") #classification at species level
taxa.print <- taxa # Reassigning taxa to a new name for downstream analysis

#Defining the rownames for the three table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

head(asv_headers)

#generating a sequence table having the defined row names
library(tidyverse)
seqs <- getSequences(seqtab.nochim)
asv_fasta <- c(rbind(asv_headers, seqs))
head(asv_fasta)
write(asv_fasta, "stingless_ASV.fasta")

#creating a sequence table dataframe
names(seqs) <- sub(">", "", asv_headers)
seqs <- as.data.frame(seqs)
seqs <- seqs %>% rownames_to_column(var = "OTU")

#generating a feature table with newly defined row names
count_asv_tab <- t(seqtab.nochim)
row.names(count_asv_tab) <- sub(">", "", asv_headers)
write.table(count_asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#generating a taxonomy table with the newly defined row names
rownames(taxa.print) <- gsub(pattern=">", replacement="", x=asv_headers)
head(taxa.print)
write.csv(taxa.print, file="ASVs_taxonomy.csv")

#Creating a matrix of the taxonomy and feature tables for creation of a phyloseq object
library(phyloseq)
#taxonomy table
TAX = tax_table(taxa.print)

#feature table
OTU = otu_table(count_asv_tab, taxa_are_rows = TRUE)

library(dplyr)
#Reading the sample meta data into R
sdata <- read.csv("stingless_bee_sample_metadata.csv", sep = ',', header = TRUE)
colnames(sdata) <- c("sample_id", "species")
sdata1 <- sdata %>% remove_rownames %>% column_to_rownames(var="Sample.id")
samdata = sample_data(sdata1)

#creating a phyloseq object
physeq = phyloseq(OTU, TAX, samdata)
physeq

#filtering the unwanted sequences
physeq0 <- subset_taxa(physeq, (Order!="Chloroplast") | is.na(Order))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Phylum!="Chloroflexi") | is.na(Phylum))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Family!="Mitochondria") | is.na(Family))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Kingdom!="Archaea") | is.na(Kingdom))
ntaxa(physeq0)

#removing the negative control
newPhyloObject = subset_samples(physeq0, sample_names(physeq0) != "NC68")

filtered_physeq <- prune_taxa(taxa_sums(newPhyloObject) > 5, newPhyloObject)
filtered_physeq
physeq1 <- subset_taxa(filtered_physeq, (is.na(Genus)))
ntaxa(physeq1)
physeq2 <- subset_taxa(filtered_physeq, (!is.na(Genus)))
ntaxa(physeq2)

library(metagMisc)
#Extracting the filtered taxonomy and feature tables for barplot plotting
tax_table <- phyloseq_to_df(physeq1, addtax = T, addtot = F, addmaxrank = F)

library(janitor)
cumulation <- tax_table %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]

silva_classified <- phyloseq_to_df(physeq2, addtax = T, addtot = F, addmaxrank = F)

library(tidyverse)
to_blast <- merge(seqs, tax_table, by = 'OTU', all = FALSE)
to_blast <- to_blast %>% select(OTU, seqs)
blast_abundance <- tax_table[,c(1, 9:63)]

library(seqRFLP)
to_blast_dada2_BSF_sequences <- dataframe2fas(to_blast, file = "to_blast_dada2_BSF_sequences.fasta")

#Running blast
blastn = "/opt/apps/blast/2.10.0+/bin/blastn"
blast_db = "16SMicrobial_v4/16SMicrobial"
input = "to_blast_dada2_BSF_sequences.fasta"
evalue = 1e-6
format = 6
max_target = 1

colnames <- c("qseqid",
              "sseqid",
              "evalue",
              "bitscore",
              "sgi",
              "sacc")
              
blast_out <- system2("/opt/apps/blast/2.10.0+/bin/blastn", 
                     args = c("-db", blast_db, 
                              "-query", input, 
                              "-outfmt", format, 
                              "-evalue", evalue,
                              "-max_target_seqs", max_target,
                              "-ungapped"),
                     wait = TRUE,
                     stdout = TRUE) %>%
  as_tibble() %>% 
  separate(col = value, 
           into = colnames,
           sep = "\t",
           convert = TRUE)

#Removing .1-9 string to get the rightful accession numbers
blast_out$sacc <- gsub(".[.1-9]$", "", blast_out$sseqid)

library(taxonomizr)
sacc <- as.vector(blast_out$sacc)
taxaId<-accessionToTaxa(sacc,"accessionTaxa.sql",version='base')
print(taxaId)
blast_taxa<-getTaxonomy(taxaId,'accessionTaxa.sql', rownames = FALSE)
print(blast_taxa)
blast_taxa <- as.data.frame(blast_taxa)
write.csv(blast_taxa, "blast_taxonomy.csv")
blast_taxa <- read.csv("blast_taxonomy.csv", sep = ',', header = TRUE)
blast_taxa$OTU <- blast_out$qseqid

#Removing the first staxids and making the OTU collumn the first
blast_taxa <- subset(blast_taxa, select = -X)
new_df <-blast_taxa %>%
  select(OTU, everything())

#checking for the presence of duplicates
anyDuplicated(new_df$OTU)
blast_results <- new_df[!duplicated(new_df$OTU),]

#changing the collumn names
colnames(blast_results) <- c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#merging the blast taxonomic clssification to blast abundance table
merged_data <- merge(blast_results, blast_abundance, by = "OTU", all = FALSE)

#merging the silva classified taxonomy with the blast classified ones
silva_blast <- as.data.frame(bind_rows(silva_classified, merged_data))

Featured_table <- silva_blast[,c(7,9:63)]
#grouping the data
group <- Featured_table %>%
  group_by(Genus)%>%
  summarise_each(funs(sum), "HA41","HA42","HA43","HA44","HA45","HA46","HA47","HA50","HG56","HG57","LG74","LG75","LG76",
                 "LG77","LN71","LN72","LN73","MB10","MB1","MB3","MB4","MB5","MB6", "MB7","MB8","MFB16","MFB17","MFB18","MFB19",
                 "MFB20","MFB26","MFB27","MFB28","MFB29","MFB30","MFR11","MFR12","MFR13","MFR14","MFR15","MFR21","MFR22","MFR23","MFR24",
                 "MFR25","DS35","DS36","DS37","DS38","DS39","DS40","DS31","DS32","DS33","DS34")
View (group)

#creating multiple dataframes for the different species
K_HA <- group[,c(1:9)]
K_HG <- group[,c(1,10,11)]
K_LG <-group[,c(1,12:15)]
K_LN <- group[,c(1,16:18)]
K_MB <-group[,c(1,19:26)]
K_MFB <-group[,c(1,27:36)]
K_MFR <- group[,c(1,37:46)]
K_DS <-group[,c(1,47:56)]

K_HA_total <- K_HA %>% adorn_totals(c("col"))
K_HA_total <- mutate(K_HA_total, K_HA=rowSums(K_HA_total[10])/8)
K_HA_total <- K_HA_total[,c(1,11)]

K_HG_total <- K_HG %>% adorn_totals(c("col"))
K_HG_total <- mutate(K_HG_total, K_HG=rowSums(K_HG_total[4])/2)
K_HG_total <- K_HG_total[,c(1,5)]

K_LG_total <- K_LG %>% adorn_totals(c("col"))
K_LG_total <- mutate(K_LG_total, K_LG=rowSums(K_LG_total[6])/4)
K_LG_total <- K_LG_total[,c(1,7)]

K_LN_total <- K_LN %>% adorn_totals(c("col"))
K_LN_total <- mutate(K_LN_total, K_LN=rowSums(K_LN_total[5])/3)
K_LN_total <- K_LN_total[,c(1,6)]

K_MB_total <- K_MB %>% adorn_totals(c("col"))
K_MB_total <- mutate(K_MB_total, K_MB=rowSums(K_MB_total[10])/8)
K_MB_total <- K_MB_total[,c(1,11)]

K_MFB_total <- K_MFB %>% adorn_totals(c("col"))
K_MFB_total <- mutate(K_MFB_total, K_MFB=rowSums(K_MFB_total[12])/10)
K_MFB_total <- K_MFB_total[,c(1,13)]

K_MFR_total <- K_MFR %>% adorn_totals(c("col"))
K_MFR_total <- mutate(K_MFR_total, K_MFR=rowSums(K_MFR_total[12])/10)
K_MFR_total <- K_MFR_total[,c(1,13)]

K_DS_total <- K_DS %>% adorn_totals(c("col"))
K_DS_total <- mutate(K_DS_total, K_DS=rowSums(K_DS_total[12])/10)
K_DS_total <- K_DS_total[,c(1,13)]

#merging the above dataframes
merged <- merge(K_HA_total, K_HG_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_LG_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_LN_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_MB_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_MFB_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_MFR_total, by = "Genus", all = TRUE)
merged <- merge(merged, K_DS_total, by = "Genus", all = TRUE)

cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]

#specifying the taxa to be tabulated
to_represent <- c("Lactobacillus", "Snodgrassella", "Saccharibacter","Bifidobacterium", "Neokomagataea", "	Saccharibacter","Bombella","	Ameyamaea",
                  "Wolbachia","Nguyenibacter", "Zymobacter", "Acinetobacter", "Gluconacetobacter", "Enterococcus", "Acetobacter", "Alkanindiges", "Chryseobacterium") 

#aggregating the rest of the phyla as others
grouped_data <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% to_represent), "Others")), sum)
View(grouped_data) 

bar <- adorn_percentages(grouped_data, denominator = "col", na.rm = TRUE)

#gathering the data
bar <- bar %>%
  gather(value = "abundance", key = "sample_names", -Genus)

#ordering the data for plotting
bar$Genus <- reorder(bar$Genus, bar$abundance)
bar$Genus <- factor(bar$Genus, levels=rev(levels(bar$Genus)))

# Choosing the colours to use in the barplot
myPalette <- c('#89C5DA', "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#plotting the barplot
ggplot(bar,aes(x = sample_names, y = abundance))+geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = myPalette)

#Beta diversity estimation
#Drawing the venn diagrams
#Removing all genuses with zero hits in the dataframes
K_HA<- K_HA_total[!(K_HA_total$K_HA == 0),]
K_HG<- K_HG_total[!(K_HG_total$K_HG == 0),]
K_LG<- K_LG_total[!(K_LG_total$K_LG == 0),]
K_LN<- K_LN_total[!(K_LN_total$K_LN == 0),]
K_MB<- K_MB_total[!(K_MB_total$K_MB == 0),]
K_MFB <- K_MFB_total[!(K_MFB_total$K_MFB == 0),]
K_MFR<- K_MFR_total[!(K_MFR_total$K_MFR == 0),]
K_DS<- K_DS_total[!(K_DS_total$K_DS == 0),]

HA <- K_HA %>% select(Genus)
HG <- K_HG %>% select(Genus)
LG <- K_LG %>% select(Genus)
LN <- K_LN %>% select(Genus)
MB <- K_MB %>% select(Genus)
MFB<- K_MFB %>% select(Genus)
MFR<- K_MFR %>% select(Genus)
DS<- K_DS %>% select(Genus)


library(gplots)
vd <- list(HA, LG, MFB, MFR, DS)
names(vd) = c("HA", "LG","MFB", "MFR", "DS")
venn(vd)

#Extracting sequences to be included in the study for plotting phylogenetic trees
seq <- merge(seqs, silva_blast, by = "OTU", all = FALSE)
seq <- seq %>% select(OTU,seqs)

#converting the filtered sequences to fasta format ad writing the fasta file to the working directory
phylo_sequences <- dataframe2fas(seq, file = "phylo_sequences.fasta")

#Reading the sequences back to R
phylo_sequences <- readDNAStringSet("phylo_sequences.fasta")
names(phylo_sequences)

library(DECIPHER)
#Running multiple sequence alignment
alignment <- AlignSeqs(phylo_sequences, anchor = NA)

library(phangorn)
#constructng the phylogenetic tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) #unrooted tree
treeUPGMA <- upgma(dm) #rooted tree
fit = pml(treeNJ, data=phang.align)

#fitting the tree using the GTR model
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) 

saveRDS(fitGTR, "stingless_bee_phangorn_tree.RDS")
phangorn <- readRDS("stingless_bee_phangorn_tree.RDS")

#Extracting the tree from the GTR model
phylo_tree <- phangorn$tree
phylogenetic_tree <- phy_tree(phylo_tree)

taxonomy_table <- silva_blast[,c(1:8)] #Extracts all the taxonomic ranks from the dataframe
features <-silva_blast[,c(1,9:63)] #Extracts all the abundance collumn for the different samples from the dataframe

#taxonomy table
my_taxonomy <- taxonomy_table %>% remove_rownames %>% column_to_rownames(var="OTU")
my_taxonomy <- as.matrix(my_taxonomy)
TAX2= tax_table(my_taxonomy)

#feature table
my_feature_table <- features %>% remove_rownames %>% column_to_rownames(var="OTU")
my_feature_table <- as.matrix(my_feature_table)
OTU2 = otu_table(my_feature_table, taxa_are_rows = TRUE)

#creating a phyloseq object
physeq3 = phyloseq(TAX2, OTU2, samdata,phylogenetic_tree)
physeq3

#sample_data(physeq3)[,2] <- sample_data(physeq3)[,1]

total = median(sample_sums(physeq3))#finds median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
standardized_physeq = transform_sample_counts(physeq3, standf)#apply to phyloseq object
ntaxa(standardized_physeq)
sample_sums(standardized_physeq)

#ordinating the phyloseq object
library(vegan)
ordu = ordinate(physeq3, "PCoA", "bray")
plot_ordination(physeq3, ordu, color="species")+ geom_point(size=2) +
  scale_color_manual(values = myPalette)

#alpha diversity estimation
physeq4 <- phyloseq(TAX2, OTU2, samdata, phylogenetic_tree)
physeq4
#checking out the total read counts in the samples
reads <- sample_sums(physeq4)
reads

summary(sample_sums(physeq4))

library(microbiome)
#Extracting the otu table from the phyloseq object and plotting the rarefaction curve
otu_tab <- t(abundances(physeq4))
p <- vegan::rarecurve(otu_tab, 
                      step = 50, label = FALSE, 
                      sample = min(rowSums(otu_tab), 
                                   col = "blue", cex = 0.6))

set.seed(9242)  

#calculatin an even sampling depth for all the samples
rarefied <- rarefy_even_depth(physeq4, sample.size = 106927)
rarefied

#calculating the alpha diversity
diversity <- alpha(rarefied, index = "all")
diversity <- rownames_to_column(diversity, "sample_id")

#Extracting the sample metadata from the phyloseq object
sdata1 <- meta(physeq4)
sdata1 <- rownames_to_column(sdata1, "sample_id")

#Extracting the shannon diversity index
shannon <- diversity %>% select(sample_id, diversity_shannon)

shannon_editted <- merge(shannon, sdata1, by = "sample_id", all = TRUE)


#confirming if the shannon indices are normally distributed
shapiro.test(shannon_editted$diversity_shannon)

library(ggpubr)

#plotting the boxplots for the shannon index data
p <- ggplot(shannon_editted, aes(x=species, y=diversity_shannon)) + geom_boxplot(aes(fill = species)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + stat_compare_means()

p





