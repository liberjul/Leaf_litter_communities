library(ggplot2)
library(vegan)
library(phyloseq)
library(ggpubr)
library(lme4)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
sl <- c(1:9, 11, 13:14, 16:21, 25:36, 38, 40:44, 48:55, 57:58) # Slice for non-negative or failed samples
otu_dat <- read.table("./Data/all_OTUS_R1_clean.txt", sep="\t", header=TRUE) # Read OTU table with contaminants removed
rownames(otu_dat) <- otu_dat[,1] # Rename rows with OTU number
otu_dat <- otu_dat[,order(colnames(otu_dat))] # Reorder by OTU number
map <- t(read.table("./Data/DEM_map_mod3.csv", fill = TRUE, sep=",", header=TRUE)) # Read the map file
colnames(map) <- map[1,] # Colnames are now by Sample ID
rownames(map) <- c("SampleID", rownames(map)[2:23]) # Rownames are now site/sample variables
map_wo_negs <- map[,sl] # Keep non-negative samples for map table
colnames(otu_dat) <- c("OTU_ID", colnames(map)) # OTU_ID as first column, sample IDs as rest of columns
otu_dat_wo_negs <- as.matrix(otu_dat[,sl + 1]) # Keep non-negative samples for OTU table
rare_otu <- t(rrarefy(t(otu_dat_wo_negs), # Rarefy by the lowest read count in non-negative sample
                      min(colSums(otu_dat_wo_negs))))

syn <- c(2, 4, 6, 8, 11, 13, 15, 17, 19, 21) # Synthetic indexes
cot <- c(1, 3, 5, 7, 9, 14, 16, 18, 20) # Cotton indexes

swab_df <- data.frame(Swab_type=c(rep("Cotton", length(cot)),
                                  rep("Synthetic", length(syn))),
                      Read_count = c(colSums(otu_dat[cot+1]), colSums(otu_dat[syn+1])),
                      Leaf = c(1:5, 7:10, 1:10))
hist(swab_df$Read_count)
var(swab_df$Read_count)
mean(swab_df$Read_count)


MDS_dat_swab <- metaMDS(t(rare_otu[,1:19]), distance = "bray")
MDS_points <- MDS_dat_swab$points
MDS_dat_df <- as.data.frame(MDS_points)
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs[,1:19]))

m1 <- glmer.nb(Read_count ~ Swab_type + (1|Leaf), swab_df)
m2 <- glmer(Read_count ~ Swab_type + (1|Leaf), swab_df, family = "poisson")
pchisq(2 * (logLik(m1) - logLik(m2)), df = 1, lower.tail = FALSE)
summary(m1)
summary(m2)
exp(coef(m1)$Leaf$Swab_typeSynthetic)
est <- confint(m1)
est
exp(est)

dist_swab <- vegdist(t(rare_otu[,1:19]))
adonis(dist_swab ~ Swab_type, MDS_dat_df[1:19,])


bd_out_swab <- betadisper(dist_swab, as.factor(MDS_dat_df[1:19,"Swab_type"]))
model <- TukeyHSD(bd_out_swab)
model
