library(VennDiagram)
library(RColorBrewer)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
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

otu_dat_trim <- rare_otu # Copy OTU table
otu_dat_trim[rare_otu < 5] <- 0 # Remove values less than 5
otu_dat_trim <- otu_dat_trim[rowSums(rare_otu)>9,] # Keep rows with more than 9 total
otu_dat_trim # Check table

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

# Take non-zeroed reads
Ep_names <- rownames(pie_df)[pie_df$`Leaf swab` > 0]
Ep_names <- Ep_names[!is.na(Ep_names)]
En_names <- rownames(pie_df)[pie_df$Leaf > 0]
En_names <- En_names[!is.na(En_names)]
Li_names <- rownames(pie_df)[pie_df$Litter > 0]
Li_names <- Li_names[!is.na(Li_names)]
So_names <- rownames(pie_df)[pie_df$Soil > 0]
So_names <- So_names[!is.na(So_names)]

# Color palette
myCol <- brewer.pal(4, "Pastel2")
# Make and ave diagram
venn.diagram(
  x = list(Ep_names, En_names, Li_names, So_names),
  category.names = c("Epiphytes", "Endophytes", "Litter", "Soil"),
  filename = "./Figures/OTU_venn_diagram.png",
  output=TRUE,
  lwd = 2,
  fill=myCol,
  cex=0.8,
  fontface="bold",
  fontfamily="sans",
  cat.cex=0.8,
  cat.fontface="bold",
  cat.fontfamily="sans"
)
