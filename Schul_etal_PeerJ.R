##White Band Disease 2017 - 2019 Samples for PeerJ Publication

#Monica D Schul
#Anya Brown

#What version of R
R.version()

#Install Packages
install.packages("rlang")
install.packages("fastqcr")
install.packages('ggplot2')
install.packages('vegan')
install.packages('tidyverse')
install.packages('Rmisc', dependencies = TRUE)
install.packages('dplyr')
install.packages('plyr')
install.packages("permute")
install.packages("lattice")
install.packages("rlang")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")


#### Load libraries ####
library(fastqcr)
library(dada2) 
library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)
library(Rmisc)
library(dplyr)
library(dplyr); packageVersion("dplyr")
library(plyr)

## set working directory
setwd("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021")


##Set a path to the 2019 sequences
########
path1 <-("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/Sequences/Sequences2019") #set path to the folder where your sequences are
fns <- list.files(path1)
fns
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
# 
# # Specify the full path to the fnFs and fnRs
fnFs <- sort(list.files(path1, pattern="_R1_cut.fastq", full.names = TRUE))

fnRs <- sort(list.files(path1, pattern="_R2_cut.fastq", full.names = TRUE))
# # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 
# ##Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1])
# 
# #Perform filtering and trimming
# #Assign the filenames for the filtered fastq.gz files.
# #Make directory and filenames for the filtered fastqs
filt_path <- file.path(path1, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
# 
# 
# # Filter the forward and reverse reads
# 
# #for our sequencing data we get 150 basepairs from the forward and the reverse reads, this is different from the tutorial, here we kept the whole sequence length and did not truncate it 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
# #What is truncQ?
# #What is rm.phix?
# 
# 
# 
# #v .6
# # Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# # Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)
# 
# #visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
# 
# #Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# # Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# 
# #Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# 
# #Inspecting the dada-class object returned by dada:
# dadaFs[[1]]
# 
# #Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# # Inspect the merger data.frame from the first sample
head(mergers[[1]])
tail(mergers[[1]])
# #Construct sequence table
seqtab2019 <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
# dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2019)))

getN <- function(x) sum(getUniques(x))
track2019 <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2019))
colnames(track2019) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track2019) <- sample.names
head(track2019)
write.table(track2019, "dada_read_stats2019.txt",sep="\t",col.names=NA)

saveRDS(seqtab2019, "~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/Sequences/seqtab2019.rds") ##update path as needed


#### 2nd sequencing run
path2 <- ("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/SeptSequences") ##update path as needed
list.files(path2)
fnFs <- sort(list.files(path2, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path2, pattern="_R2_cut.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Perform filtering and trimming
filt_path <- file.path(path2, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
# Learn the Error Rates, it TAKES TIME!
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
# Dereplicate the filtered fastq files
existsF <- file.exists(filtFs)
existsR <- file.exists(filtRs)
derepFs <- derepFastq(filtFs[existsF], verbose=TRUE)
derepRs <- derepFastq(filtRs[existsR], verbose=TRUE)
names(derepFs) <- sample.names[sample.names != "D44-V4-26"]
names(derepRs) <- sample.names[sample.names != "D44-V4-26"]
# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# Inspecting the dada-class object returned by dada:
dadaFs[[1]]
# Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# Construct sequence table
seqtabSept <- makeSequenceTable(mergers)
dim(seqtabSept)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabSept)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtabSept))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats2_Cayman.txt",sep="\t",col.names=NA)
track <- read.table("dada_read_stats2_Cayman.txt")
saveRDS(seqtabSept, "~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/Sequences/seqtabSept.rds") ##update path as needed


####2017 sequences
path3 <-("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/Sequences/Sequences2017") #set path to the folder where your sequences are
fns <- list.files(path3)
fns
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
# 
# # Specify the full path to the fnFs and fnRs
fnFs <- sort(list.files(path3, pattern="_R1_cut.fastq", full.names = TRUE))

fnRs <- sort(list.files(path3, pattern="_R2_cut.fastq", full.names = TRUE))
# # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 
# ##Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1])
# 
# #Perform filtering and trimming
# #Assign the filenames for the filtered fastq.gz files.
# #Make directory and filenames for the filtered fastqs
filt_path <- file.path(path3, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
# 
# 
# # Filter the forward and reverse reads
# 
# #for our sequencing data we get 150 basepairs from the forward and the reverse reads, this is different from the tutorial, here we kept the whole sequence length and did not truncate it 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
# #What is truncQ?
# #What is rm.phix?
# 
# 
# 
# #v .6
# # Forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
# # Reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)
# 
# #visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
# 
# #Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# # Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
# 
# #Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# 
# #Inspecting the dada-class object returned by dada:
# dadaFs[[1]]
# 
# #Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# # Inspect the merger data.frame from the first sample
head(mergers[[1]])
tail(mergers[[1]])
# #Construct sequence table
seqtab2017 <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
# dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2017)))

getN <- function(x) sum(getUniques(x))
track2017 <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2017))
colnames(track2017) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track2017) <- sample.names
head(track2017)
write.table(track2017, "dada_read_stats2017.txt",sep="\t",col.names=NA)

saveRDS(seqtab2017, "~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/Sequences/seqtab2017.rds") ##update path as needed

##Merge all the sequences to make 1 table
st1 <- readRDS("seqtab2017.rds")
st2 <- readRDS("seqtab2019.rds")
st3 <- readRDS("seqtabSept.rds")
st.all <- mergeSequenceTables(st1, st2, st3, repeats="sum") # You will get the message "Duplicated sample names detected in the sequence table row names." to let you know that there are duplicate names across samples - it is not an error, just a message.

#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)## 120 15004
sum(seqtab.nochim)/sum(st.all) #0.9610253


# Track reads through the pipeline
# As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline
rowSums(seqtab.nochim)
# need to write this out to add to dada read stats

# SAVE the non-chimeric sequence variant table SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS
saveRDS(seqtab.nochim, file="~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBDCaymanSamplesCombined2021.rds") ##update path as needed
# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
seqtab.nochim <- readRDS("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/WBDCaymanSamplesCombined2021.rds") ##update path as needed
getwd()
#Assign taxonomy
#Make sure the appropriate database is available in the folder that you direct it to, if you have it in the folder were you have set the directory (setwd() above), then you can just add the database
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)


# Removing sequence rownames for display only
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save the table 
write.table(taxa,"taxa_table_10042021.txt",sep="\t")
taxa<- read.delim("taxa_table_10042021.txt")

### Replacing NAs with the lowest taxonomic level

taxa_df <- as.data.frame(taxa)
dim(taxa_df) #15004     6
#replace NAs with the deepest taxonomy available
taxa_df$Kingdom <- as.character(taxa_df$Kingdom)
king_na <- which(is.na(taxa_df$Kingdom))
taxa_df[king_na, "Kingdom"] <- 'Unknown'

taxa_df$Phylum <- as.character(taxa_df$Phylum)
phy_na <- which(is.na(taxa_df$Phylum))
taxa_df[phy_na, "Phylum"] <- taxa_df$Kingdom[phy_na] 

taxa_df$Class <- as.character(taxa_df$Class)
cl_na <- which(is.na(taxa_df$Class))
taxa_df[cl_na, "Class"] <- taxa_df$Phylum[cl_na]

taxa_df$Order <- as.character(taxa_df$Order)
ord_na <- which(is.na(taxa_df$Order))
taxa_df[ord_na, "Order"] <- taxa_df$Class[ord_na]

taxa_df$Family <- as.character(taxa_df$Family)
fam_na <- which(is.na(taxa_df$Family))
taxa_df[fam_na, "Family"] <- taxa_df$Order[fam_na]

taxa_df$Genus <- as.character(taxa_df$Genus)
gen_na <- which(is.na(taxa_df$Genus))
taxa_df[gen_na, "Genus"] <- taxa_df$Family[gen_na]
#put the data into phyloseq

write.table(taxa_df,"taxa_table_nonas_10042021.txt",sep="\t")
taxa_df <-read.table("taxa_table_nonas_10042021.txt")

########### ########### ########### 
########### Phyloseq ##########
########### ########### ########### 
#put the data into phyloseq

library(phyloseq)

#read in the sample data or medata data file
samp.dat <- read.csv("metadat_comparison_2021.csv")
head(samp.dat)
#check dimensions
dim(samp.dat) ##111  9
samp.dat <- subset(samp.dat, Rep_incl == 1)
#add  labels to rownames to import into phyloseq
rownames(samp.dat) <- samp.dat$SampleID

#Read in cleaned seqtab file
seqtab.nochim <- readRDS("~/coral disease/Microbial Sampling_Away_Healthy_Lesion/Data/with2017/PostReviews_May/WBD October 2021/WBDCaymanSamplesCombined2021.rds") ##update path as needed

#read in taxa data frame
taxa_df <-read.table("taxa_table_nonas_10042021.txt")
#convert taxa dataframe to a matrix for puting into phyloseq
taxa_mat <- as.matrix(taxa_df)


library(phyloseq)
#put in phyloseq
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), #sequences
               sample_data(samp.dat), #sample data
               tax_table(taxa_mat)) #taxa

#remova chloroplasts and mitochondria and Eukaryota
ps_nomito <- subset_taxa(ps, Family !="Mitochondria")
ps_nomito_chloro <- subset_taxa(ps_nomito, Order !="Chloroplast")
ps_nomito_chloro<- subset_taxa(ps_nomito_chloro, Kingdom !="Eukaryota")

#Check sample data
sums <- sample_sums(ps_nomito_chloro)
range(sums) ## 0 232727
mean(sums)  #  76269.54
write.csv(sums, "sample_sums_beforeprune_10042021.csv")

#### Filter Steps #####
#remove samples with sequence sums less than 800 
ps_nomito_chloro = prune_samples(sample_sums(ps_nomito_chloro) > 800, ps_nomito_chloro)

sums_post <- sample_sums(ps_nomito_chloro)

range(sums_post) #1127 232727
mean(sums_post)  #80391.81
write.csv(sums_post, "sample_sums_postprune_10042021.csv")


##Create a color pallet
##############Color Palettes#############

#color palettes
my20colors<-c("#c26162","#d689c2","#d64142","#6db643","#9c58cb","#cdab39","#626fdf","#5cba78","#bd3fa0","#497535","#de76dd","#949345","#6e58a2","#cd6529","#46aed7","#d94681","#49b39c","#9e4e77","#bd814c","#798fd5")
my50colors<-c("#e48ab2","#54c44e","#9e40b6","#a8be30","#9162db","#429a38","#db70db","#85bb4e","#5365d5","#ceb535","#618ae9","#e9a034","#7553a1","#39c685","#b53492","#7ac283","#df3e7d","#58bea1","#c83046","#43ccd7","#dc4e2e","#56a6d9","#e17b31","#446ba9","#b6851a","#c99ae7","#a9a93e","#dc62aa","#487229","#b878be","#778526","#9292d5","#a9471c","#2da0a1","#de6161","#2d7a55","#a64364","#5c9a5b","#985382","#a5b069","#9f4d45","#d0aa5c","#e18387","#70692c","#f08c68","#9b7a2d","#c0745a","#dba577","#8c5c2e","#bc7335") 
my50colors.v2<-rev(my50colors)
my50colors.v3<-c("#5365d5","#ceb535","#618ae9","#e9a034","#7553a1","#39c685","#b53492","#7ac283","#df3e7d","#58bea1","#c83046","#43ccd7","#dc4e2e","#56a6d9","#e17b31","#446ba9","#b6851a","#c99ae7","#a9a93e","#dc62aa","#487229","#b878be","#778526","#9292d5","#a9471c","#2da0a1","#de6161","#2d7a55","#a64364","#5c9a5b","#985382","#a5b069","#9f4d45","#d0aa5c","#e18387","#70692c","#f08c68","#9b7a2d","#c0745a","#dba577","#8c5c2e","#bc7335","#e48ab2","#54c44e","#9e40b6","#a8be30","#9162db","#429a38","#db70db","#85bb4e") 
my50colors.v4<-c("#db70db","#85bb4e","#5365d5","#ceb535","#618ae9","#e9a034","#7553a1","#39c685","#b53492","#7ac283","#df3e7d","#58bea1","#c83046","#43ccd7","#dc4e2e","#56a6d9","#e17b31","#446ba9","#b6851a","#c99ae7","#a9a93e","#dc62aa","#487229","#b878be","#778526","#9292d5","#a9471c","#2da0a1","#de6161","#2d7a55","#a64364","#5c9a5b","#985382","#a5b069","#9f4d45","#d0aa5c","#e18387","#70692c","#f08c68","#9b7a2d","#c0745a","#dba577","#8c5c2e","#bc7335","#e48ab2","#54c44e","#9e40b6","#a8be30","#9162db","#429a38") 
my50colors.v5<-c("#a8be30","#9162db","#429a38","#db70db","#85bb4e","#5365d5","#ceb535","#618ae9","#e9a034","#7553a1","#39c685","#b53492","#7ac283","#df3e7d","#58bea1","#c83046","#43ccd7","#dc4e2e","#56a6d9","#e17b31","#446ba9","#b6851a","#c99ae7","#a9a93e","#dc62aa","#487229","#b878be","#778526","#9292d5","#a9471c","#2da0a1","#de6161","#2d7a55","#a64364","#5c9a5b","#985382","#a5b069","#9f4d45","#d0aa5c","#e18387","#70692c","#f08c68","#9b7a2d","#c0745a","#dba577","#8c5c2e","#bc7335","#e48ab2","#54c44e","#9e40b6") 

my70colors<-c(my20colors, my50colors)
my70colors.v2<-rev(my70colors)
my100colors <- c(my50colors.v3,my50colors.v4)


#Set color palette for samples
colors <- c("AppHealthy" = "#E69F00", "Healthy" = "#D55E00","NearDisease" = "#999999")
colorsgenotype <- c("BLACK"="gray","BLUE"="blue","GREEN"="green","RED"="red","YELLOW"="yellow")


## Figure 2, Making a relative abundance plot
#relative abundance
ps_rel<-transform_sample_counts(ps_nomito_chloro, function(OTU) OTU/sum(OTU))

#filter taxa to mean abundance greater than 0.5%
ps10<-filter_taxa(ps_rel, function(x) mean(x) >.005, TRUE)
ntaxa(ps10) #21

#Abundance > 1%
pdf("RelAbundplot.pdf", height = 5, width = 7)
ps2.ord=plot_bar(ps10, fill="Order", x = "Group_Gen" ) #to make a barplot the colors based on Order

ps2.ord +geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+ theme_bw()+theme(legend.position = "top", axis.text.x = element_text(angle = 90), strip.background = element_blank()) +facet_grid(~Treatment+Year, scales = "free") + scale_fill_manual(values = my100colors) + xlab("Sample") + ylab("Relative Abundance")
dev.off()


############ Alpha Diversity #############
##############################################
#Alpha Diversity - number of unique OTUS (Observed) or taking into account the proportin of OTUs across samples (Shannon)

library(data.table)
ps_nomito_chloro_2 <-  prune_samples(sample_sums(ps_nomito_chloro) > 800, ps_nomito_chloro)

ps_nomito_chloro_r <- rarefy_even_depth(ps_nomito_chloro) #this is the only time I rarefy..to estimate the number of unique groups given the same sampling effort. THis is contentious, but I think for looking at ecological questions it make sense. I do not use rarefied data for any other type of analysis, expecially not when taking a compositional approach

#estimate richness on unrarefied data
rich <- estimate_richness(ps_nomito_chloro_2, measures = c("Shannon", "InvSimpson"))


#For merging, name ID columns the same thing
rich$SampleID <- rownames(rich)

#merged the alpha diversity measures and the sample data frame by their shared column - SampleID

alpha <- merge(rich, samp.dat, by = "SampleID")
str(alpha)

alpharmelt <- reshape2::melt(id.vars = c("SampleID","Treatment","Group","Genotypes","Year","Colony","Location","Rep_incl","Group_Gen"), data= alpha)
colnames(alpharmelt)[10] <- "Alpha"
colnames(alpharmelt)[11] <- "Diversity"


####### Figure 3 data ######

colors2 <- c("AppHealthy 2019" = "#E69F00", "Healthy 2019" = "red4","NearDisease 2019" = "#999999", "Healthy 2017" = "#ff8423")
#tell ggplot the order of factors to plot
alpharmelt$Treatment <- factor(alpharmelt$Treatment, levels = c("Healthy","AppHealthy","NearDisease"))
alpharmelt$Treatment2 <- paste(alpharmelt$Treatment,alpharmelt$Year)

alpharmelt$Treatment2 <- factor(alpharmelt$Treatment2, levels = c("Healthy 2017", "Healthy 2019","AppHealthy 2019","NearDisease 2019"))

pdf("Fig3_ShannonandSimpson.pdf",width = 8, height = 6)
ps= ggplot(alpharmelt, aes(x = Treatment2, y = Diversity)) + geom_boxplot() +  geom_point(size = 4, aes(color = Treatment2)) + theme_bw()
ps+scale_color_manual(values = colors2) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16)) + facet_grid(vars(Alpha), scales = "free") + xlab("Treatment")
dev.off()  


alpha$Group <- as.factor(alpha$Group)
alpha_r2019 <- subset(alpha, Year == 2019)

shan.lm <- lm(Shannon~Treatment+ Group, data = alpha_r2019)
summary(shan.lm)  # Group p >0.05

plot(resid(shan.lm))
anova(shan.lm)


shan.lm
shan.lm2 <- lm(Shannon~Treatment, data = alpha_r2019) #F = 4.6058 P = 0.02195 
plot(resid(shan.lm2))
anova(shan.lm2)

######################
###Tukey Test

anova1 <- aov(shan.lm2)
TukeyHSD(anova1)

########################################
#Kruskal-Wallis chi-squared

sh.krus <- kruskal.test(Shannon ~ Treatment, data= alpha_r2019)
sh.krus
#Kruskal-Wallis chi-squared = 5.4361, df = 2, p-value = 0.066

#analysis not in manuscript
#######################################################
## Dunn test
install.packages("FSA")
library(FSA)
#order by median
alpha_r2019$Treatment = factor(alpha_r2019$Treatment, levels = c("AppHealthy","NearDisease","Healthy"))
PH = dunnTest(Shannon ~ Treatment,
              data=alpha_r2019,
              method="bh")

######################
###Tukey Test

anova1 <- aov(shan.lm2)
TukeyHSD(anova1)

########################################
#Kruskal-Wallis chi-squared

sh.krus <- kruskal.test(Shannon ~ Treatment, data= alpha_r2019)
sh.krus
#Kruskal-Wallis chi-squared = 5.4361, df = 2, p-value = 0.066
#analysis not in manuscript
#######################################################
## Dunn test
install.packages("FSA")
library(FSA)
#order by median
alpha_r2019$Treatment = factor(alpha_r2019$Treatment, levels = c("AppHealthy","NearDisease","Healthy"))
PH = dunnTest(Shannon ~ Treatment,
              data=alpha_r2019,
              method="bh")


#############################################################
###InvSimpson
isim.lm <- lm(InvSimpson~Treatment+ Group, data = alpha_r2019) # Group p <0.05
summary(isim.lm)

plot(resid(isim.lm))
anova(isim.lm)

isim.lm2 <- lm(InvSimpson~Treatment, data = alpha_r2019)
anova2<- aov(isim.lm2)
TukeyHSD(anova2)


####################################################
########Shannon diversity for Healthy data ######

#tell ggplot the order of factors to plot
alpha_rHealthy <- subset(alpha, Treatment == "Healthy")

sh.lm.years <- lm(Shannon~Year+Genotypes, data = alpha_rHealthy)
plot(resid(sh.lm.years))
anova(sh.lm.years)



####################################################################
#############################################################
##Beta Diversity
sample_data(ps_rel)
ps_rel2019 <- subset_samples(ps_rel, Year == 2019)
ps_relHealthy <- subset_samples(ps_rel, Treatment == "Healthy")


#####Betadispersion #########
#how communities vary within a treatment (like, the distance from the sample point in multivariate space to the centroid or center of all the points in that treatment)####

library(vegan) #using a different package, so we need to take out phyloseq objects and make them into martrices or tables that can be used outside of phyloseq

######## 2019 ########
#make the otu table for just 2019 data
otu_rel2019 <- as.matrix(otu_table(ps_rel2019))

#make the taxa table for 2019
taxa_rel2019 <- tax_table(ps_rel2019)

#Make sample data table
metadat_rel2019 <- as.data.frame(as.matrix(sample_data(ps_rel2019)))
metadat_rel2019

#make a distance matrix, and specificy the method - in this case, bray curtis 
m2 <- vegdist(otu_rel2019, method="bray")
mod <-betadisper(m2, metadat_rel2019$Treatment) #model that that we use to test the differences in distance to the centroid within our groups in the column Treatment


anova(mod) #run an anova on the distances based on 2019 treatments 
# F=1.8352, p=0.1804
mod


#to get the distance values
disp<- as.data.frame(mod$distances)

disp$SampleID <- rownames(disp)
colnames(disp)[1] <- "distance"

#merge 
disp2 <- merge(disp,metadat_rel2019, by = "SampleID")
str(disp2)

disp2$Treatment <- factor(disp2$Treatment, levels = c("Healthy","AppHealthy","NearDisease"))


#### Figure 4 part A####

dis.plot <- ggplot(disp2,aes(x = Treatment, y = distance)) + geom_boxplot() + geom_point(size = 8, aes(color= Treatment)) + theme(text = element_text(size = 16)) + theme_bw()+scale_color_manual(values=colors) + ylab("Distance to centroid")
dis.plot

########## Dispersion of healthy data across years ########

#####boxplot across years for healthy data 


library(vegan) #using a different package, so we need to take out phyloseq objects and make them into martrices or tables that can be used outside of phyloseq

#make the otu table for healthy
otu_relHealthy <- as.matrix(otu_table(ps_relHealthy))

#make the taxa table
taxa_relHealthy <- tax_table(ps_relHealthy)

#Make sample data table
metadat_relHealthy <- as.data.frame(as.matrix(sample_data(ps_relHealthy)))

#make a distance matrix, and specificy the method - in this case, bray curtis 
m3 <- vegdist(otu_relHealthy, method="bray")
metadat_relHealthy$Year <- factor(metadat_relHealthy$Year, levels = c("2017","2019"))
mod2 <-betadisper(m3, metadat_relHealthy$Year) #model that that we use to test the differences in distance to the centroid within our groups in the column Treatment


anova(mod2) # based on healthy samples
#F=1.2496, p=0.2801

#to get the distance values
disph<- as.data.frame(mod2$distances)

disph$SampleID <- rownames(disph)
colnames(disph)[1] <- "distance"

#merge 
disp3 <- merge(disph,metadat_relHealthy, by = "SampleID")
str(disp3)

####include in figure 4 ####
h.plot <- ggplot(disp3,aes(x = Year, y = distance)) + geom_boxplot() + geom_point(size = 8, aes(color= Year)) + theme_bw()+scale_color_manual(values=c("#ff8423","#D55E00"))
h.plot 

year.lm <- lm(distance~Year+Genotypes, data = disp3)
summary(year.lm)
anova(year.lm)
plot(resid(year.lm))

##Final figure 4
###Combine part A and B
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("ggpubr")
library(ggpubr)

ggarrange(dis.plot, h.plot + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

######Permanova to test differences among treatments #####
##########################################################

##### 2019 Data####
#differences across treatments in 2019
m22019 <- vegdist(otu_rel2019, method="bray")

#analyse Treatment with 2019 data
perm2019<-adonis(m22019~Treatment,strata = metadat_rel2019$Group, data= metadat_rel2019, perm = 999)
print(perm2019) 

#Call:
#  adonis(formula = m22019 ~ Treatment, data = metadat_rel2019,      permutations = 999, strata = metadat_rel2019$Group) 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999
#Terms added sequentially (first to last)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Treatment  2    2.9208 1.46042  6.6207 0.34626  0.001 ***
#Residuals 25    5.5146 0.22059         0.65374           
#Total     27    8.4355                 1.00000           
---
  ##  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
  ###### Healthy data across years ####
#analyze year, with just Health data
m2h <- vegdist(otu_relHealthy, method="bray")
#analyse Treatment with 2019 data
permh<-adonis(m2h~Year,data= metadat_relHealthy, strata = metadat_relHealthy$Genotypes, perm = 999)
print(permh) 


#Call:
#  adonis(formula = m2h ~ Year, data = metadat_relHealthy, permutations = 999,      strata = metadat_relHealthy$Genotypes) 

#Blocks:  strata 
#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#         =25Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Year       1    0.4889 0.48890  2.2358 0.12261  0.073 .
#Residuals 16    3.4987 0.21867         0.87739         
#Total     17    3.9876                 1.00000         
---
  #  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  

##Figure 5, PCOA 
library(ggforce)
sample_data(ps_rel)$Treatment2 <- paste(sample_data(ps_rel)$Treatment, sample_data(ps_rel)$Year)


colors2 <- c("AppHealthy 2019" = "#E69F00", "Healthy 2019" = "red4","NearDisease 2019" = "#999999", "Healthy 2017" = "#ff8423")


#### PCoA with ellipse around treatment and genotypes as shape

sample_data(ps_rel)$Treatment2 <- factor(sample_data(ps_rel)$Treatment2, levels = c("Healthy 2019", "AppHealthy 2019","NearDisease 2019", "Healthy 2017"))
cay_pcoa <- ordinate(ps_rel, "PCoA", "bray", weighted=TRUE)

p2 <- plot_ordination(ps_rel, cay_pcoa, color = "Treatment2") # our phyloseq object (ps_rel), and the nmds information (cay_nmds) and we are going to color the points by treatment
p2 <- p2+  geom_point(size=6, alpha=1, stroke=2, aes(shape = as.factor(Genotypes)))
p2 <- p2+ xlim(-0.6, 0.6)+ylim(-0.4, 0.8)#we make the points, select the size adn transparency
p2  <- p2 + theme_bw() + theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") #change the background to bw, the text to 16 and get rid of the ugly gridlines and change the legend position
p2 <- p2+scale_color_manual(values=c(colors2))
p2<-p2 + guides(color = guide_legend(nrow = 2, byrow=TRUE), shape =guide_legend(nrow = 3, byrow=TRUE))+
  geom_mark_ellipse()
p2
dev.off()

####DESEQ2 Figure 6
##DeSeq2 Treatment Healthy vs Disease
library("DESeq2")
packageVersion("DESeq2")

#dataset with raw sequences (not relative abundance)
HealthVSDisease <-subset_samples(ps_nomito_chloro, Treatment !="AppHealthy")
HealthyVSDisease2019 <- subset_samples(HealthVSDisease, Year !="2017")
sample_data(HealthVSDisease)$Year <- as.factor(sample_data(HealthVSDisease)$Year)


#set this up so it's phylseq object ~ Treatment
dds2 = phyloseq_to_deseq2(HealthyVSDisease2019, ~ Treatment)
dds2 <- dds2[ rowSums(counts(dds2)) > 10, ]
cts2 <- counts(dds2)
# geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
geoMeans2 <- apply(cts2, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds2 <- estimateSizeFactors(dds2, geoMeans=geoMeans2)

#this runs the test
dds2 <- DESeq(dds2, test="Wald", fitType="parametric")
resultsNames(dds2) #I use this to check how the comparison is made - so if it says Healthy Disease, and my treatment levels were Healthy and Disease, the first name will be the positive number, and the second the negative (so enriched vs not enriched)
# [1] "Intercept"                        "Treatment_NearDisease_vs_Healthy"

#isolate the results at an alpha level of 0.05
res2 = results(dds2, cooksCutoff = FALSE)
results(dds2)
alpha = 0.05
sigtab6 = res2[which(res2$padj < alpha), ] 

#create a datatable with the significant results and the taxa names. Replace the "ps_coral_stx" for whatever your dataset is

sigtab6 = cbind(as(sigtab6, "data.frame"), as(tax_table(HealthVSDisease)[rownames(sigtab6), ], "matrix"))
head(sigtab6)
dim(sigtab6) #check the dimensions to make sure it makes sense

sigtab6_df <- as.data.frame(sigtab6)
sigtab6_df$OTU <- rownames(sigtab6_df)
tax <- as.data.frame(as.matrix(tax_table(HealthVSDisease))) #replace ps_coral_stx with your dataset
tax$OTU <- rownames(tax)
sig.tax <- merge(sigtab6_df, tax, by = "OTU") #this is a dataframe with all of the taxonomic information
dim(sigtab6_df) #check the dimensions so the merge makes sense
dim(sig.tax)
# 
length(which(sigtab6_df$log2FoldChange>0)) #number of ASVs enriched in treatment level A
length(which(sigtab6_df$log2FoldChange<0)) #number of ASVs enriched in treatment level B

#Plot

#I've been making plots that show multiple levels - this case Order and Genus
sigtab6_df$OrderGenus <- paste(sigtab6_df$Order, sigtab6_df$Genus)
write.csv(sigtab6_df,"~/coral disease/WBD study/September 2022/Healthy2019_vs_Disease_log2fold.csv", row.names = FALSE)

#this statement allows for coloring points based on enrichment at the level of the treatment. If you use this, change "Removal" and "Prsent" to match your dataset. "Removal" is the level A of the treatment - so it's the positive values of Log2ratio, and "Present" are the negative values of the log2ratio in this case


sigtab6_df_20 <- subset(sigtab6_df, abs(log2FoldChange) >=20)
sigtab6_df_5 <- subset(sigtab6_df, abs(log2FoldChange) >=5)

colors2 <- ifelse(sigtab6_df$log2FoldChange>0,"Disease","Healthy")
colors_treatHD <- c("Healthy" = "#D55E00","Disease" = "#999999")
#DeSeq plot
#this plot is based on my color schemeing of making the colors show which groups are enriched in my two treatments. you can do it this way if you want, or color by something else
Healthy_vs_Disease_logChange<-ggplot(sigtab6_df, aes(y=log2FoldChange, x = OrderGenus)) +
  xlab("Genera") + geom_point(size=6, pch = 21, aes(fill=colors2)) +
  coord_flip() +geom_abline(color = "red", slope = 0)  + 
  theme_bw() + theme(text =element_text(size = 16))+
  scale_fill_manual(values = colors_treatHD)
Healthy_vs_Disease_logChange
dev.off()


#################################################################################
###################################################################################
#################################################################################
###Supplemental figures

#Figure S2
cay_pcoa <- ordinate(ps_rel, "PCoA", "bray", weighted=TRUE)

p <- plot_ordination(ps_rel, cay_pcoa)# our phyloseq object (ps_rel), and the nmds information (cay_nmds) and we are going to color the points by treatment
p <- p+  geom_point(size=6, alpha=1, stroke=2, aes(shape=Treatment2))#we make the points, select the size adn transparency
p  <- p + theme_bw() + theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") #change the background to bw, the text to 16 and get rid of the ugly gridlines and change the legend position

p + guides(color = guide_legend(nrow = 2, byrow=TRUE), shape =guide_legend(nrow = 3, byrow=TRUE))+
  geom_mark_ellipse(aes(fill=Genotypes))+facet_wrap(~Group)+scale_fill_manual(values = colorsgenotype)
dev.off()

####Figure S3
######### ANCOM ######
#####################
### ANCOM FUNCTION###
## run this first###


###Need to run this first in order to run ANCOM on data
library(exactRankTests)
library(nlme)
library(ggplot2)

ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
}



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}
######

library(data.table)
#Need to run ANCOM function first
#ANCOM
#from phyloseq -use all of the data, not relative abundance, and  not filtered beyond mitochondria and choloroplasts 
#Make the otutable


#### subset data for ANCOM ####
taxa_names(ps_nomito_chloro) <- paste0("Seq", seq(ntaxa(ps_nomito_chloro)))
ps_2019 <- subset_samples(ps_nomito_chloro, Year == "2019")
ps_select_taxa <- subset_samples(ps_nomito_chloro, Year == "2019", Genus = c("Vibrio","Thalassotalea","Thalassolituus","MD3.55","Algicola","Catenococcus","Litoricola","Flavobacteriales","HIMB11"))
no_Rick <- subset_taxa(ps_nomito_chloro, Order !="Rickettsiales")
ps_H <- subset_samples(ps_nomito_chloro, Treatment == "Healthy")



######### 2019 data - comparing treatments ####
otu2019 <- as.matrix(otu_table(ps_2019))
Sample.ID <- rownames(otu2019)
otus2 <- cbind(Sample.ID, otu2019)
otud <- otus2

metadat19 <- sample_data(ps_2019) #get the sample data
metadat19 <- as.data.frame(as.matrix(metadat19)) #make into into a matrix
colnames(metadat19)[1] <- "Sample.ID" #make sure the sample names column is called Sample.ID


names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM
metadat19$Group <- as.factor(metadat19$Group)
library(dplyr)
metadat <- dplyr::select(metadat19, c("Sample.ID","Treatment","Group")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax

comparison_test_group=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=TRUE, #can't remember
                                 repeated=F, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= F, #other factors to include
                                 repeat.var=F, #repeated measure
                                 long = F, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #90% cut=0ff


res <- comparison_test_group$W.taxa #taxa that sifnificantly vary across factor level of interest
res2 <- res[which(res$detected_0.9 == "TRUE"),] #to get the groups that show the highest W and are different at the 90% level 

res3 <- res[which(res$detected_0.8 == "TRUE"),] #to get the groups that show the highest W and are different at the 80% level 

colnames(res2)[1] <- "OTU"
colnames(res3)[1] <- "OTU"
tax2019 <- as.data.frame(as.matrix((tax_table(ps_2019))))
tax2019$OTU <- rownames(tax2019)
sig.otus <- merge(tax2019, res2, by = "OTU")
sig.otus2 <- merge(tax2019, res3, by = "OTU")
write.csv(sig.otus, "ANCOM_det90_2019_treatments_adjgroup_Oct10052021.csv")

write.csv(sig.otus2, "ANCOM_det80_2019_treatments_adjgroup_Oct10052021.csv")

sig.otus <- read.csv("ANCOM_det90_2019_treatments_adjgroup_Oct10052021.csv")

ps_2019_rel <- transform_sample_counts(ps_2019, function(OTU) OTU/sum(OTU))

#dat_rel <- tax_glom(ps_rel2019, taxrank = "Genus")

#melt the relative abundance data and skip the next bit
datm_rel <- psmelt(ps_2019_rel)


sig.otus1 <- subset(sig.otus, W_stat > 122)
colnames(sig.otus1)[2] <- "Genus"


#straight from the melt function
datm_rel$Genus <- gsub("-",".", datm_rel$Genus) 
sig.otus <- dplyr::select(sig.otus, OTU, W_stat)
#merge by OTU
trycrel2b <- merge(sig.otus,datm_rel, by = "OTU")
#check
dim(trycrel2b)

#this combines the Family, genus and OTU
trycrel2b$GenusOTU <- paste(trycrel2b$Family, trycrel2b$Genus, trycrel2b$OTU)

#average abundance
trycrelsum <- Rmisc::summarySE(trycrel2b, measurevar = "Abundance",groupvars = c("Treatment","OTU","Genus", "Family", "W_stat"))


#highest W stats fo rplotting
trycrelsum <- subset(trycrelsum, W_stat > 130)

#themes
themes <-  theme_bw() + theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

summary(trycrelsum)
trycrelsum$Treatment <- factor(trycrelsum$Treatment, levels = c("Healthy","AppHealthy","NearDisease"))
trycrelsum$Genus <- factor(trycrelsum$Genus, levels = c("Vibrio","Thalassotalea","Thalassolituus","MD3.55","Algicola","Catenococcus","Litoricola","Flavobacteriales","HIMB11"))

#Figure 6 ##
pdf("Figure6_TreatmentASV.pdf". width = 10)
ggplot(trycrelsum, aes(x = Treatment, y = Abundance)) +
  geom_point(aes(color = Treatment), size =3, position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymax = Abundance + se, ymin = Abundance-se, width = 0.2)) +
  ggtitle("Significant ASVs Treatment")+
  themes + facet_wrap(~Genus + Family, scales = "free_y")+ 
  theme( axis.text.x = element_blank(),strip.text.x = element_text(size =8), legend.position = "bottom")+
  scale_color_manual(values = colors, labels = c("Healthy", "Apparently Healthy", "Disease")) + 
  ylab("Relative Abundance")
dev.off()

#Supplemental table of data #

trycrelsum2 <- reshape2::dcast(trycrelsum, OTU+Genus+Family+W_stat~Treatment, value.var = "Abundance")
trycrelsum3 <- reshape2::dcast(trycrelsum, OTU+Genus+Family+W_stat~Treatment, value.var = "se")
write.csv(trycrelsum2,"Wstat_means_2019_otu_Jan.csv")
write.csv(trycrelsum3,"Wstat_se_2019_otu_Jan.csv")
getwd()



###Figure S4
##DeSeq2 Treatment Healthy 2019 vs Healthy 2017
library("DESeq2")
packageVersion("DESeq2")

#dataset with raw sequences (not relative abundance)
no_Rick <- subset_taxa(ps_nomito_chloro, Order !="Rickettsiales")
no_RickH <-subset_samples(no_Rick, Treatment == "Healthy")
sample_data(no_RickH)$Year <- as.factor(sample_data(no_RickH)$Year)


#set this up so it's phylseq object ~ Treatment
dds2 = phyloseq_to_deseq2(no_RickH, ~ Year)
dds2 <- dds2[ rowSums(counts(dds2)) > 10, ]
cts2 <- counts(dds2)
# geoMeans = apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
geoMeans2 <- apply(cts2, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds2 <- estimateSizeFactors(dds2, geoMeans=geoMeans2)

#this runs the test
dds2 <- DESeq(dds2, test="Wald", fitType="parametric")
resultsNames(dds2) #I use this to check how the comparison is made - so if it says Healthy Disease, and my treatment levels were Healthy and Disease, the first name will be the positive number, and the second the negative (so enriched vs not enriched)
# [1]  "Intercept"         "Year_2019_vs_2017"

#isolate the results at an alpha level of 0.05
res2 = results(dds2, cooksCutoff = FALSE)
results(dds2)
alpha = 0.05
sigtab6 = res2[which(res2$padj < alpha), ] 

#create a datatable with the significant results and the taxa names. Replace the "ps_coral_stx" for whatever your dataset is

sigtab6 = cbind(as(sigtab6, "data.frame"), as(tax_table(no_RickH)[rownames(sigtab6), ], "matrix"))
head(sigtab6)
dim(sigtab6) #check the dimensions to make sure it makes sense

sigtab6_df <- as.data.frame(sigtab6)
sigtab6_df$OTU <- rownames(sigtab6_df)
tax <- as.data.frame(as.matrix(tax_table(no_RickH))) #replace ps_coral_stx with your dataset
tax$OTU <- rownames(tax)
sig.tax <- merge(sigtab6_df, tax, by = "OTU") #this is a dataframe with all of the taxonomic information
dim(sigtab6_df) #check the dimensions so the merge makes sense
dim(sig.tax)
# 
length(which(sigtab6_df$log2FoldChange>0)) #number of ASVs enriched in treatment level A
length(which(sigtab6_df$log2FoldChange<0)) #number of ASVs enriched in treatment level B

#Plot

#I've been making plots that show multiple levels - this case Order and Genus
sigtab6_df$OrderGenus <- paste(sigtab6_df$Order, sigtab6_df$Genus)
write.csv(sigtab6_df,"C:/Users/schul/Documents/PhD/Papers/WBD/Healthy2019_vs_Healthy2017_log2fold.csv", row.names = FALSE)

#this statement allows for coloring points based on enrichment at the level of the treatment. If you use this, change "Removal" and "Prsent" to match your dataset. "Removal" is the level A of the treatment - so it's the positive values of Log2ratio, and "Present" are the negative values of the log2ratio in this case


sigtab6_df_20 <- subset(sigtab6_df, abs(log2FoldChange) >=20)
sigtab6_df_5 <- subset(sigtab6_df, abs(log2FoldChange) >=5)


colorsYear <- ifelse(sigtab6_df$log2FoldChange>0,"Healthy 2017","Healthy 2019")
colors_Year <- c("Healthy 2017" = "#ff8423","Healthy 2019" = "#D55E00")
#DeSeq plot
#this plot is based on my color schemeing of making the colors show which groups are enriched in my two treatments. you can do it this way if you want, or color by something else
Healthy2019_vs_Healthy2017_logChange<-ggplot(sigtab6_df, aes(y=log2FoldChange, x = OrderGenus)) +
  xlab("Order and Genus") + geom_point(size=6, pch = 21, aes(fill=colorsYear)) +
  coord_flip() +geom_abline(color = "red", slope = 0)  + 
  theme_bw() + theme(text =element_text(size = 16))+
  scale_fill_manual(values = colors_Year)
Healthy2019_vs_Healthy2017_logChange
dev.off()



