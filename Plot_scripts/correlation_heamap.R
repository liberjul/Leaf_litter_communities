library(ggplot2)
library(ggpubr)
library(gplots)
library(vegan)
library(scales)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communites")
sl <- c(1:9, 11, 13:21, 25:36, 38, 40:44, 48:55, 57:58) # Slice for non-negative or failed samples
otu_dat <- read.table("./Data/all_OTUS_R1_clean.txt", sep="\t", header=TRUE) # Read OTU table with contaminants removed
rownames(otu_dat) <- otu_dat[,1] # Rename rows with OTU number
otu_dat <- otu_dat[,order(colnames(otu_dat))] # Reorder by OTU number
colnames(otu_dat) <- c("OTU_ID", colnames(map)) # OTU_ID as first column, sample IDs as rest of columns
otu_dat_wo_negs <- as.matrix(otu_dat[,sl + 1]) # Keep non-negative samples for OTU table
map <- t(read.table("./Data/DEM_map_mod3.csv", fill = TRUE, sep=",", header=TRUE)) # Read the map file
colnames(map) <- map[1,] # Colnames are now by Sample ID
rownames(map) <- c("SampleID", rownames(map)[2:23]) # Rownames are now site/sample variables
map_wo_negs <- map[,sl] # Keep non-negative samples for map table
rare_otu <- t(rrarefy(t(otu_dat_wo_negs), # Rarefy by the lowest read count in non-negative sample
              min(colSums(otu_dat_wo_negs))))

cor_mat_otu <- cor(rare_otu) # Correlation matrix
samp_names <- map_wo_negs["Short_name",] # Sample names for axis labels
png("corr_heatmap.png", width = 800, height = 800) # Open png for output
heatmap.2(cor_mat_otu, scale= "none", # Make heatmap
          trace = "none", density.info = "none",
          col=viridis_pal(),
          labRow = samp_names,
          labCol = samp_names)
dev.off()
