###### Using the Dada2 pipeline of MiSeq data raw reads ####
### Schul et al. 2021 #####
# # Forward and reverse fastq filenames have format: SAMPLENAME_R1_cut.fastq and SAMPLENAME_R2_cut.fastq
# # Samplename is everything before the first underscore
# # based on https://benjjneb.github.io/dada2/tutorial.html

#### Load libraries ####
library(dada2) #
packageVersion("dada2") #version 1.16.0
library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)
install.packages('Rmisc', dependencies = TRUE)
library(Rmisc)
library(dplyr)
library(dplyr); packageVersion("dplyr")
library(plyr)
#source(biocLite())
#biocLite("DESeq2")



path <- #set path to the folder where your sequences are
fns <- list.files(path)
fns
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
# 
# # Specify the full path to the fnFs and fnRs
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names = TRUE))
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
 filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
 if(!file_test("-d", filt_path)) dir.create(filt_path)
 filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
 filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
# 
# 
# # Filter the forward and reverse reads
# 
# #for our sequencing data we get 150 basepairs from the forward and the reverse reads, this is different from the tutorial, here we kept the whole sequence length and did not truncate it 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
# head(out)
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
seqtab <- makeSequenceTable(mergers) ## The sequences being tabled vary in length.
# dim(seqtab)
# 
# # Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# 
# 
# 
# #Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# #Track reads through the pipeline
# #As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
 getN <- function(x) sum(getUniques(x))
 track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
 colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
 rownames(track) <- sample.names
 head(track)
 write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)
# 
# #####SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name
# saveRDS(seqtab.nochim, file="disease_seqtab.nochim.rds")
# # RELOAD THE SAVED INFO FROM HERE (if you have closed the project):
# 
# #move to the correct folder with teh path and then include the file name all in quootations
# # each of the / indicate a folder, the ~/ is used to start from the home directory, this a little different between PCs and Macs. For a PC - it would like more like this c:/Documents/my/working/directory")

seqtab.nochim <- readr::read_rds("healthydisease_seqtab.nochim.rds")
#Assign taxonomy
#Make sure the appropriate database is available in the folder that you direct it to, if you have it in the folder were you have set the directory (setwd() above), then you can just add the database
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)


# Removing sequence rownames for display only
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save the table 
write.table(taxa.print,"taxa_table.txt",sep="\t")
taxa2<- read.delim("taxa_table_edited_20172019.txt")

### Replacing NAs with the lowest taxonomic level

taxa_df <- as.data.frame(taxa2)
dim(taxa_df) #12133  6
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

write.table(taxa_df,"taxa_table_nonas.txt",sep="\t")

#put the data into phyloseq

library(phyloseq)

#read in the sample data or medata data file
samp.dat <- read.csv("metadat_comparison_update.csv")
#check dimensions
dim(samp.dat) ##111  8

#add  labels to rownames to import into phyloseq
rownames(samp.dat) <- samp.dat$SampleID

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
range(sums) ##0 to 178647
mean(sums)  # 62896.24
write.csv(sums, "sample_sums_beforeprune.csv")

#### Filter Steps #####
#remove samples with sequence sums less than 100 
ps_nomito_chloro = prune_samples(sample_sums(ps_nomito_chloro) > 100, ps_nomito_chloro)

#only include 1 rep per colony, and all 2019 data
ps_nomito_chloro = subset_samples(ps_nomito_chloro, Rep_incl == "1")
#subset 100 and more

sums_post <- sample_sums(ps_nomito_chloro)
range(sums_post) #194 122419
mean(sums_post)  # 51866.6
write.csv(sums_post, "sample_sums_postprune.csv")


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
colors <- c("AppHealthy" = "magenta", "Healthy" = "salmon","NearDisease" = "magenta4")

####### Figure 2 Abundance Plot #####
### to filter taxa or samples
library(ggplot2)
ps_nomito_chloro = prune_samples(sample_sums(ps_nomito_chloro) > 0, ps_nomito_chloro)
sums <-sample_sums(ps_nomito_chloro)
range(sums)
mean(sums)
#filter taxa to mean abundance greater than 10 sequences for plotting
ps10<-filter_taxa(ps_nomito_chloro, function(x) mean(x) >10, TRUE)
ntaxa(ps10) #195

#relative abundance
ps2<-transform_sample_counts(ps10, function(OTU) OTU/sum(OTU))

A <- subset_samples(ps2, Treatment == "AppHealthy")
H <- subset_samples(ps2, Treatment == "Healthy")
D <- subset_samples(ps2, Treatment == "NearDisease")

ps3 <- merge_phyloseq(D,A,H)

#plot by Order, legend on bottom, facet by treatments

sample_data(ps3)$Treatment <- factor(sample_data(ps3)$Treatment, levels = c("Healthy", "AppHealthy","NearDisease"))
ps2.ord=plot_bar(ps3, fill="Order", x = "SampleID" ) #to make a barplot the colors based on Order
ps2.ord +geom_bar(aes(fill=Order), stat="identity",position="stack")+theme(strip.text=element_text(face="bold"))+ theme( axis.text.x = element_blank())+theme(legend.position = "bottom") + facet_grid(~Treatment+Year, scales = "free") + scale_fill_manual(values= my100colors)
ps2.ord$ps3$Sample <-factor(ps2.ord$ps3$Sample, levels = "H","A","D")

###### NMDS ##### 


#relative abundance
ps_rel<-transform_sample_counts(ps_nomito_chloro, function(OTU) OTU/sum(OTU))

##subsetting data for further analysis
sample_data(ps_rel)
ps_rel2019 <- subset_samples(ps_rel, Year == 2019)
ps_relHealthy <- subset_samples(ps_rel, Treatment == "Healthy")

cay_nmds <- ordinate(ps_rel, "NMDS", "bray", weighted=TRUE)

### Figure 5 NMDS 2019 ####
p <- plot_ordination(ps_rel2019, cay_nmds, color="Treatment") # our phyloseq object (ps_rel), and the nmds information (cay_nmds) and we are going to color the points by treatment
p <- p+  geom_point(size=10, alpha=1, stroke=2, aes(shape = as.factor(Treatment))) #we make the points, select the size adn transparency
p  <- p + theme_bw() + theme(text = element_text(size = 18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") #change the background to bw, the text to 16 and get rid of the ugly gridlines and change the legend position
p <- p+ scale_color_manual(values=colors)
p

#### NMDS just healthy ####
sample_data(ps_relHealthy)$Year <- as.factor(sample_data(ps_relHealthy)$Year)
p <- plot_ordination(ps_relHealthy, cay_nmds, color="Year") # our phyloseq object (ps_rel), and the nmds information (cay_nmds) and we are going to color the points by treatment
p <- p+  geom_point(size=10, alpha=1, stroke=2, aes(shape = as.factor(Year))) #we make the points, select the size adn transparency
p  <- p + theme_bw() + theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") #change the background to bw, the text to 16 and get rid of the ugly gridlines and change the legend position
p <- p+ scale_color_manual(values=c("salmon","salmon"))
p

##### NMDS split years####
sample_data(ps_rel)$Treatment <- factor(sample_data(ps_rel)$Treatment, levels = c("Healthy", "AppHealthy","NearDisease"))
p <- plot_ordination(ps_rel, cay_nmds, color = "Treatment") # our phyloseq object (ps_rel), and the nmds information (cay_nmds) and we are going to color the points by treatment
p <- p+  geom_point(size=6, alpha=1, stroke=2, aes(shape = as.factor(Year))) #we make the points, select the size adn transparency
p  <- p + theme_bw() + theme(text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") #change the background to bw, the text to 16 and get rid of the ugly gridlines and change the legend position
p <- p+facet_grid(~Year)+scale_color_manual(values=colors)
p


#####Betadispersion - how communities vary within a treatment (like, the distance from the sample point in multivariate space to the centroid or center of all the points in that treatment)####

library(vegan) #using a different package, so we need to take out phyloseq objects and make them into martrices or tables that can be used outside of phyloseq

######## 2019 ########3
#make the otu table for just 2019 data
otu_rel2019 <- as.matrix(otu_table(ps_rel2019))

#make the taxa table for 2019
taxa_rel2019 <- tax_table(ps_rel2019)

#Make sample data table
metadat_rel2019 <- as.data.frame(as.matrix(sample_data(ps_rel2019)))
metadat_rel2019

#make a distance matrix, and specificy the method - in this case, bray curtis 
m2 <- vegdist(otu_rel2019, method="bray")
metadat_rel2019$Treatment <- factor(metadat_rel2019$Treatment, levels = c("Healthy", "AppHealthy","NearDisease"))
mod <-betadisper(m2, metadat_rel2019$Treatment) #model that that we use to test the differences in distance to the centroid within our groups in the column Treatment


anova(mod) #run an anova on the distances based on 2019 treatments 
# F=0.301, p=0.7429
mod


#to get the distance values
disp<- as.data.frame(mod$distances)

disp$SampleID <- rownames(disp)
colnames(disp)[1] <- "distance"

#merge 
disp2 <- cbind(disp,metadat_rel2019$Treatment)
str(disp2)
colnames(disp2)[3] <- "Treatment"

#### Figure 4 ####
ggplot(disp2,aes(x = Treatment, y = distance)) + geom_boxplot() + geom_point(size = 8, aes(color= Treatment)) + theme_bw()+scale_color_manual(values=colors)


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
#F=0.491, p=0.4936

#to get the distance values
disp<- as.data.frame(mod2$distances)

disp$SampleID <- rownames(disp)
colnames(disp)[1] <- "distance"

#merge 
disp2 <- cbind(disp,metadat_relHealthy$Year)
str(disp2)
colnames(disp2)[3] <- "Year"

ggplot(disp2,aes(x = Year, y = distance)) + geom_boxplot() + geom_point(size = 8, aes(color= Year)) + theme_bw()+scale_color_manual(values=c("salmon","salmon"))


##########################################################
######Permanova to test differences among treatments #####
##########################################################
 
##### 2019 Data####
#differences across treatments in 2019
m22019 <- vegdist(otu_rel2019, method="bray")

#analyse Treatment with 2019 data
perm2019<-adonis(m22019~Treatment,data= metadat_rel2019, perm = 999)
print(perm2019) ###R2= Treatment: 0.32978, Residuals: 0.67022, p=0.001

###### Healthy data across years ####
#analyze year, with just Health data
m2h <- vegdist(otu_relHealthy, method="bray")
#analyse Treatment with 2019 data
permh<-adonis(m2h~Year,data= metadat_relHealthy, perm = 999)
print(permh) ## R2= Year: 0.09372, Residuals, 0.90628, p= 0.153



##############################################
############ Alpha Diversity #############
##############################################
#Alpha Diversity - number of unique OTUS (Observed) or taking into account the proportin of OTUs across samples (Shannon)

library(data.table)
ps_nomito_chloro_2 <-  prune_samples(sample_sums(ps_nomito_chloro) > 100, ps_nomito_chloro)

ps_nomito_chloro_r <- rarefy_even_depth(ps_nomito_chloro) #this is the only time I rarefy..to estimate the number of unique groups given the same sampling effort. THis is contentious, but I think for looking at ecological questions it make sense. I do not use rarefied data for any other type of analysis, expecially not when taking a compositional approach

#estimate richness on unrarefied data
rich_r <- estimate_richness(ps_nomito_chloro_2, measures = "Shannon")
#make into a datatable
#alphdiv_r <- data.table(rich_r)

#For merging, name ID columns the same thing
rich_r$SampleID <- rownames(rich_r)

#merged the alpha diversity measures and the sample data frame by their shared column - SampleID

alpha_r <- merge(rich_r, samp.dat, by = "SampleID")
str(alpha_r)

####### Figure 3 2019 data ######

###
#tell ggplot the order of factors to plot
alpha_r$Treatment <- factor(alpha_r$Treatment, levels = c("Healthy","AppHealthy","NearDisease"))
p= ggplot(alpha_r[which(alpha_r$Year == "2019"),], aes(x = Treatment, y = Shannon)) + geom_boxplot() +  geom_point(size = 4, aes(color = Treatment)) + theme_bw()
p+scale_color_manual(values = colors) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())



########Shannon diversity for Healthy data ######

#tell ggplot the order of factors to plot
alpha_r$Year <- factor(alpha_r$Year, levels = c("2017","2019"))
p= ggplot(alpha_r[which(alpha_r$Treatment == "Healthy"),], aes(x = Year, y = Shannon)) + geom_boxplot() +  geom_point(size = 8, aes(color = Year)) + theme_bw()
p+ scale_color_manual(values=c("salmon","salmon")) + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())



#evaluate significance of Shannon diversity 

metadat_rel <- as.data.frame(as.matrix(sample_data(ps_rel)))
alpha_r <- merge(rich_r, metadat_rel, by = "SampleID")

#### ANOVA based on 2019 treatments ###
div_summary2 <- lm(Shannon~Treatment, data = alpha_r[which(alpha_r$Year==2019),])
anova(div_summary2) # F=5.6663, p= 0.009983


#### ANOVA based on Just healthy data ###
div_summary3 <- lm(Shannon~Year, data = alpha_r[which(alpha_r$Treatment=="Healthy"),])
anova(div_summary3) # F=1.2379, p=0.2823


###################################################
#################Core Microbiome############
###################################################
library("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
library(microbiome)
#data(ps2) # should be with

A.all <- subset_samples(ps_rel, Treatment == "AppHealthy")
H.all <- subset_samples(ps_rel, Treatment == "Healthy")
D.all <- subset_samples(ps_rel, Treatment == "NearDisease")
H.genus <- tax_glom(H.all, taxrank = "Genus")
A.genus <- tax_glom(A.all, taxrank = "Genus")
D.genus <- tax_glom(D.all, taxrank = "Genus")
H2019 <- subset_samples(H, Year == 2019)
h2017 <- subset_samples(H, Year == 2017)
H2019 <- subset_samples(H.genus, Year == 2019)
h2017 <- subset_samples(H.genus, Year == 2017)
head(prevalence(ps2, detection = 00/100, sort = TRUE))
pseq.coreA <- core(A.genus, detection = 0, prevalence = 0.99)
pseq.coreH <- core(H2019, detection = 0, prevalence = 0.99)
pseq.coreH2017 <- core(h2017, detection = 0, prevalence = 0.99)
pseq.coreD <- core(D.genus, detection = 0, prevalence = 0.99)

hmelt <- psmelt(pseq.coreH)
Hmelt2017 <- psmelt(pseq.coreH2017)
Amelt <- psmelt(pseq.coreA)
Dmelt <- psmelt(pseq.coreD)

### Means and SE for Core microbiomes ####
h.2019_prev99 <- summarySE(hmelt, measurevar = "Abundance", groupvars = "Genus")
h.2017_prev99 <- summarySE(Hmelt2017, measurevar = "Abundance", groupvars = "Genus")
A._prev99 <- summarySE(Amelt, measurevar = "Abundance", groupvars = "Genus")
D._prev99 <- summarySE(Dmelt, measurevar = "Abundance", groupvars = "Genus")


########################
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
#group the otus based on family
#this is probaby a little overkill, but it's how I've done it in the past - all you need is the otutable really, so probably that could have been done in one step, but i like this for future plot-making
dat <- tax_glom(ps_nomito_chloro, taxrank = "Genus") #at the Family level

#melt the data, so it's like a dataframe
datm <- psmelt(dat)


#Cast the new datatable with columns that are of interest
datc <- data.table::dcast(datm, SampleID +Treatment + Year ~ Genus, value.var = 'Abundance', fun.aggregate = sum)

#### subset data for ANCOM ####
taxa_names(ps_nomito_chloro) <- paste0("Seq", seq(ntaxa(ps_nomito_chloro)))
ps_2019 <- subset_samples(ps_nomito_chloro, Year == "2019")
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

library(dplyr)
metadat <- select(metadat19, c("Sample.ID","Treatment")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax

comparison_test_treat=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                                 Vardat=map_test, #calling the metadata
                                 adjusted=TRUE, #can't remember
                                 repeated=F, #repeated measure
                                 main.var="Treatment", #main variable or fator
                                 adj.formula= FALSE, #other factors to include
                                 repeat.var=F, #repeated measure
                                 long = F, #longitudinal study
                                 multcorr=2,
                                 sig=0.05, #significance level
                                 prev.cut=0.90) #90% cut=0ff

res <- comparison_test_treat$W.taxa #taxa that sifnificantly vary across factor level of interest
res2 <- res[which(res$detected_0.9 == "TRUE"),] #to get the groups that show the highest W and are different at the 90% level 
colnames(res2)[1] <- "OTU"
tax2019 <- as.data.frame(as.matrix((tax_table(ps_2019))))
tax2019$OTU <- rownames(tax2019)
sig.otus <- merge(tax2019, res2, by = "OTU")


ps_rel2019<- subset_samples(ps_rel, Year == 2019)
dat_rel <- tax_glom(ps_rel2019, taxrank = "Genus")

#melt the relative abundance data and skip the next bit
datm_rel <- psmelt(dat_rel)


sig.otus1 <- subset(sig.otus, W_stat > 122)
colnames(sig.otus1)[2] <- "Genus"


#straight from the melt function
datm_rel$Genus <- gsub("-",".", datm_rel$Genus) 

#merge by genus
trycrel2b <- merge(sig.otus1,datm_rel, by = "Genus")
#check
dim(trycrel2b)
dim(datm_rel)

#average abundance
trycrelsum <- Rmisc::summarySE(trycrel2b, measurevar = "Abundance",groupvars = c("Treatment","Genus", "W_stat"))


#highest W stats fo rplotting
trycrelsum <- subset(trycrelsum, W_stat > 130)

#themes
themes <-  theme_bw() + theme(text = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


trycrelsum$Treatment <- factor(trycrelsum$Treatment, levels = c("Healthy","AppHealthy","NearDisease"))


#Figure 6 ##
ggplot(trycrelsum, aes(x = Treatment, y = Abundance)) + geom_point(aes(color = Treatment), size =3, position = position_dodge(width = 0.2)) + geom_errorbar(aes(ymax = Abundance+se, ymin = Abundance-se)) + ggtitle("Relative Abundance of Genera differences based on ANCOM results")+themes+ facet_wrap(~Genus, scales = "free_y")+ theme( axis.text.x = element_blank(),strip.text.x = element_text(size =8), legend.position = "bottom")+scale_color_manual(values = colors, labels = c("Healthy", "Apparently Healthy", "Disease")) + ylab("Relative Abundance")


#Supplemental table of data #

trycrelsum2 <- reshape2::dcast(trycrelsum_taxnames, OTU+Genus+Family+W_stat~Treatment, value.var = "RelAbund")
trycrelsum3 <- reshape2::dcast(trycrelsum_taxnames, OTU+Genus+Family+W_stat~Treatment, value.var = "se")
write.csv(trycrelsum2,"Wstat_means_2019_otu.csv")
write.csv(trycrelsum3,"Wstat_se_2019_otu.csv")

################################################### 
################# ANCOM Healthy data #####
################################################### 

#OTU H table
otuH <- as.matrix(otu_table(ps_H))
Sample.ID <- rownames(otuH)
otus2 <- cbind(Sample.ID, otuH)
otud <- otus2

metadatH <- sample_data(ps_H) #get the sample data
metadatH <- as.data.frame(as.matrix(metadatH)) #make into into a matrix
colnames(metadatH)[1] <- "Sample.ID" #make sure the sample names column is called Sample.ID


names(otud) <- make.names(names(otud)) #get the names from the table for 
otu_test <- otud #rename otud to otu_test, for syntax in ANCOM

library(dplyr)
metadat <- select(metadatH, c("Sample.ID","Year", "Genotypes")) # use select to only use treatment columns of interest
map_test <- metadat #rename map_TEst
Vardat <- map_test #specify that this for Vardat - ANCOM syntax

comparison_test_H=ANCOM.main(OTUdat=otu_test, #calling the OTU table
                             Vardat=map_test, #calling the metadata
                             adjusted=TRUE, #can't remember
                             repeated=F, #rpeated measure
                             main.var="Year", #main variable or fator
                             adj.formula= FALSE, #other factors to include
                             repeat.var=F, #repeated measure
                             long = F, #longitudinal study
                             multcorr=2,
                             sig=0.05, #significance level
                             prev.cut=0.90) #90% cut=0ff


resH <- comparison_test_H$W.taxa #taxa that sifnificantly vary across factor level of interest
resH2 <- resH[which(resH$detected_0.9 == "TRUE"),] #to get the groups that show the highest W and are different at the 90% level 
colnames(resH2)[1] <- "OTU"
taxH <- as.data.frame(as.matrix((tax_table(ps_H))))
taxH$OTU <- rownames(taxH)
sig.otusHealthy <- merge(taxH, resH2, by = "OTU")

colnames(sig.otusHealthy)[2] <- "Genus"

#going from ps_rel
dat_rel_all <- tax_glom(ps_rel, taxrank = "Genus")
dat_rel_all_H <- subset_samples(dat_rel_all, Treatment == "Healthy")
datm_relH <- psmelt(dat_rel_all_H)

datm_relH$Genus <- gsub("-",".", datm_relH$Genus) 

#select highest W stat
sig.otusHealthy <- subset(sig.otusHealthy, W_stat > 148)

#using teh relative abundance straight from psmelt
trycrel_2Healthy_b <- merge(datm_relH, sig.otusHealthy, by = "Genus")


## Top Genera Mean relative abundance tables ###
trycrelsumH <- Rmisc::summarySE(trycrel_2Healthy_b, measurevar = "Abundance",groupvars = c("Year","Genus", "W_stat"))
trycrelsumH2 <- dcast(trycrelsumH, Genus+W_stat~Year, value.var = "Abundance")
trycrelsumH3 <- dcast(trycrelsumH, Genus+W_stat~Year, value.var = "se")

write.csv(trycrelsumH2,"Wstat_means_Healthy_psmelt_rel_alb.csv")
write.csv(trycrelsumH3,"Wstat_se_Healthy_psmelt_rel_alb.csv")

