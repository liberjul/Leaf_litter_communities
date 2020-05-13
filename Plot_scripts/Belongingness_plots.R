library(ggplot2)
library(ggpubr)
library(patchwork)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))

prop_otu_dat <- rare_otu # Assign the table to new df to maintain shape

for (i in 1:dim(rare_otu)[2]){ # for each column
  prop_otu_dat[,i] <- rare_otu[,i]/(sum(rare_otu[,i])*wts_wo_negs[i]) # Calculate proportions and scale by # of samples per substrate
}
prop_otu_dat_trim <- prop_otu_dat[rowSums(prop_otu_dat) > 4 * 0.0001,]

prop_otu_dat_trim
pie_df <- data.frame(dim(prop_otu_dat_trim)[1]) # New dataframe for proportions
for (subst in unique(map_wo_negs["Substrate",])){ # For each substrate
  pie_col <- rowSums(prop_otu_dat_trim[,map_wo_negs["Substrate",] == subst]) # Take total reads per substrate
  pie_df <- cbind(pie_df, pie_col) # Add column to df
}

pie_df <- pie_df[,2:5] # take the data columns

colnames(pie_df) <- unique(map_wo_negs["Substrate",]) # Rename columns with substrate
row_sums <- rowSums(pie_df) # rowsums to scale
for (i in colnames(pie_df)){ # scale abundance by total in reads in dataset
  pie_df[,i] <- pie_df[,i]/row_sums
}
pie_df

# Make dataframe with axes corresponding to proportions of reads by substrate
am.coords <- as.data.frame(as.matrix(cbind(pie_df["Epi"] + pie_df["Soil"] - pie_df["Endo"] - pie_df["Lit"],
                        pie_df["Epi"] - pie_df["Soil"] + pie_df["Endo"] - pie_df["Lit"]), ncol=2))
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
  labs(x="(Epi + Soil) - (Endo + Lit)", y="(Endo + Epi) - (Lit + Soil)")
  # theme(axis.ticks=element_blank(),
  #       axis.line=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.text.y=element_blank())
g
ggsave("./Figures/density_shared_OTUs.png", g, width = 6, height = 6, units="in")

# Histograms to show belongingness of 2-substrate OTUs where one substrate is litter
N_L_hist <- ggplot(am.coords[am.coords$Epi_Litter == -1 & abs(am.coords$End_Soil) < 1,],
                   aes(End_Soil)) +
  geom_histogram() + labs(x="Litter - Endophyte", y = "OTU Count", tag = "B") + theme_pubr() + xlim(-1, 1)
# S_P_hist <- ggplot(am.coords[am.coords$Epi_Litter == 1 & abs(am.coords$End_Soil) < 1,],
#                    aes(End_Soil)) +
#   geom_histogram() + labs(x="Soil - Epiphyte") + theme_pubr()  + xlim(-1, 1)
L_S_hist <- ggplot(am.coords[am.coords$End_Soil == -1 & abs(am.coords$Epi_Litter) < 1,],
                   aes(Epi_Litter)) +
  geom_histogram() + scale_y_log10() + labs(x="Litter - Soil", y = "OTU Count") + theme_pubr()
# N_P_hist <- ggplot(am.coords[am.coords$End_Soil == 1,], aes(Epi_Litter)) +
#   geom_histogram() + scale_y_log10() + labs(x="Endophyte - Epiphyte") + theme_pubr()
L_P_hist <- ggplot(am.coords[am.coords$End_Soil == am.coords$Epi_Litter & abs(am.coords$End_Soil) < 1,],
                   aes(End_Soil)) +
  geom_histogram() + labs(x="Litter - Epiphyte", y = "OTU Count")  + theme_pubr()  + xlim(-1, 1)
L_P_hist

density_excl <- g + labs(tag="A") | (N_L_hist / L_P_hist /  L_S_hist ) #/ S_P_hist / N_P_hist )
density_excl
ggsave("./Figures/density_shared_and_excl_OTUs.png", density_excl, width = 12, height = 8, units="in")
