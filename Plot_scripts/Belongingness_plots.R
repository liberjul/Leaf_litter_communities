library(ggplot2)
library(ggpubr)
library(patchwork)

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

otu_dat_trim <- rare_otu # Copy OTU table
otu_dat_trim[rare_otu < 5] <- 0 # Remove values less than 5
otu_dat_trim <- otu_dat_trim[rowSums(rare_otu)>9,] # Keep rows with more than 9 total
otu_dat_trim # Check table

for (i in 1:dim(otu_dat_trim)[2]){ # scale OTU data by the total reads in substrate
  otu_dat_trim[,i] <- otu_dat_trim[,i]/sum(map_wo_negs["Soil_Leaf_Litter_Leaf_swab",] == map_wo_negs["Soil_Leaf_Litter_Leaf_swab", i])
}

pie_df <- data.frame(dim(otu_dat_trim)[1]) # New dataframe for proportions
for (subst in unique(map_wo_negs["Soil_Leaf_Litter_Leaf_swab",])){ # For each substrate
  pie_col <- rowSums(otu_dat_trim[,map_wo_negs["Soil_Leaf_Litter_Leaf_swab",] == subst]) # Take total reads per substrate
  pie_df <- cbind(pie_df, pie_col) # Add column to df
}
pie_df <- pie_df[,2:5] # take the data columns
colnames(pie_df) <- unique(map_wo_negs["Soil_Leaf_Litter_Leaf_swab",]) # Rename columns with substrate
row_sums <- rowSums(pie_df) # rowsums to scale
for (i in colnames(pie_df)){ # scale abundance by total in reads in dataset
  pie_df[,i] <- pie_df[,i]/row_sums
}

# Make dataframe with axes corresponding to proportions of reads by substrate
am.coords <- as.data.frame(as.matrix(cbind(pie_df["Leaf swab"] + pie_df["Soil"] - pie_df["Leaf"] - pie_df["Litter"],
                        pie_df["Leaf swab"] - pie_df["Soil"] + pie_df["Leaf"] - pie_df["Litter"]), ncol=2))
colnames(am.coords) <- c("Epi_Litter", "End_Soil")

# Dataframe to make labels on the corners
corners <- data.frame(x=c(-1, -1, 1, 1),
                      y=c(-1.1, 1.1, -1.1, 1.1),
                      labels = c("Litter", "Endophytes", "Soil", "Epiphytes"))

g <- ggplot(am.coords, aes(x=Epi_Litter, y=End_Soil)) +
  geom_text(data=corners, aes(x=x,y=y,label=labels)) +
  coord_fixed(ratio=1, xlim=c(-1.2, 1.2), ylim=c(-1.1, 1.1)) +
  geom_point(alpha=0.1) +
  theme_pubr() +
  labs(x=NULL, y=NULL) +
  theme(axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())
g
ggsave("./Figures/density_shared_OTUs.png", g, width = 6, height = 6, units="in")

# Histograms to show belongingness of 2-substrate OTUs where one substrate is litter
N_L_hist <- ggplot(am.coords[am.coords$Epi_Litter == -1 & abs(am.coords$End_Soil) < 1,],
                   aes(End_Soil)) +
  geom_histogram() + labs(x="Litter - Endophyte") + theme_pubr() + xlim(-1, 1)
# S_P_hist <- ggplot(am.coords[am.coords$Epi_Litter == 1 & abs(am.coords$End_Soil) < 1,],
#                    aes(End_Soil)) +
#   geom_histogram() + labs(x="Soil - Epiphyte") + theme_pubr()  + xlim(-1, 1)
L_S_hist <- ggplot(am.coords[am.coords$End_Soil == -1 & abs(am.coords$Epi_Litter) < 1,],
                   aes(Epi_Litter)) +
  geom_histogram() + scale_y_log10() + labs(x="Litter - Soil") + theme_pubr()
# N_P_hist <- ggplot(am.coords[am.coords$End_Soil == 1,], aes(Epi_Litter)) +
#   geom_histogram() + scale_y_log10() + labs(x="Endophyte - Epiphyte") + theme_pubr()
L_P_hist <- ggplot(am.coords[am.coords$End_Soil == am.coords$Epi_Litter & abs(am.coords$End_Soil) < 1,],
                   aes(End_Soil)) +
  geom_histogram() + labs(x="Litter - Epiphyte")  + theme_pubr()  + xlim(-1, 1)
L_P_hist

density_excl <- g | (N_L_hist / L_S_hist /  L_P_hist )#/ S_P_hist / N_P_hist )
density_excl
ggsave("./Figures/density_shared_and_excl_OTUs.png", density_excl, width = 12, height = 8, units="in")
