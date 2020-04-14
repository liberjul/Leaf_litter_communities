library(VennDiagram)
library(RColorBrewer)
library(dplyr)
library(vegan)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))

otu_dat_trim <- rare_otu # Copy OTU table
# otu_dat_trim[rare_otu < 5] <- 0 # Remove values less than 5
# otu_dat_trim <- otu_dat_trim[rowSums(rare_otu)>9,] # Keep rows with more than 9 total
# otu_dat_trim # Check table

pie_df <- data.frame(dim(otu_dat_trim)[1]) # New dataframe for proportions
for (subst in unique(map_wo_negs["Substrate",])){ # For each substrate
  pie_col <- rowSums(otu_dat_trim[,map_wo_negs["Substrate",] == subst]) # Take total reads per substrate
  pie_df <- cbind(pie_df, pie_col) # Add column to df
}

pie_df <- pie_df[,2:5] # take the data columns
pie_df <- pie_df %>% filter_all(any_vars(abs(.) > 5))

colnames(pie_df) <- unique(map_wo_negs["Substrate",]) # Rename columns with substrate
row_sums <- rowSums(pie_df) # rowsums to scale
for (i in colnames(pie_df)){ # scale abundance by total in reads in dataset
  pie_df[,i] <- pie_df[,i]/row_sums
}
dim(pie_df)

# Take non-zeroed reads
Ep_names <- rownames(pie_df)[pie_df$Epi > 0]
Ep_names <- Ep_names[!is.na(Ep_names)]
En_names <- rownames(pie_df)[pie_df$Endo > 0]
En_names <- En_names[!is.na(En_names)]
Li_names <- rownames(pie_df)[pie_df$Lit > 0]
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

shared_otu_df <- data.frame(total= c(dim(pie_df %>% filter(Epi > 0))[1],
                                     dim(pie_df %>% filter(Endo > 0))[1],
                                     dim(pie_df %>% filter(Lit > 0))[1],
                                     dim(pie_df %>% filter(Soil > 0))[1]),
                            uniq = c(sum(rowSums(pie_df) == pie_df$Epi),
                                     sum(rowSums(pie_df) == pie_df$Endo),
                                     sum(rowSums(pie_df) == pie_df$Lit),
                                     sum(rowSums(pie_df) == pie_df$Soil)),
                            epi_shared = c(NA,
                                            dim(intersect(pie_df %>% filter(Epi > 0), pie_df %>% filter(Endo > 0)))[1],
                                            dim(intersect(pie_df %>% filter(Epi > 0), pie_df %>% filter(Lit > 0)))[1],
                                            dim(intersect(pie_df %>% filter(Epi > 0), pie_df %>% filter(Soil > 0)))[1]),
                            endo_shared = c(dim(intersect(pie_df %>% filter(Endo > 0), pie_df %>% filter(Epi > 0)))[1],
                                            NA,
                                            dim(intersect(pie_df %>% filter(Endo > 0), pie_df %>% filter(Lit > 0)))[1],
                                            dim(intersect(pie_df %>% filter(Endo > 0), pie_df %>% filter(Soil > 0)))[1]),
                            lit_shared = c(dim(intersect(pie_df %>% filter(Lit > 0), pie_df %>% filter(Epi > 0)))[1],
                                           dim(intersect(pie_df %>% filter(Lit > 0), pie_df %>% filter(Endo > 0)))[1],
                                           NA,
                                           dim(intersect(pie_df %>% filter(Lit > 0), pie_df %>% filter(Soil > 0)))[1])
                              )
shared_otu_df$shared <- shared_otu_df$total-shared_otu_df$uniq

shared_percent <- vector(length = 4)
for (i in 1:4){
  shared_percent[i] <- sprintf("%.0f (%.0f",
                               shared_otu_df$shared[i],
                               shared_otu_df$shared[i]/shared_otu_df$total[i]*100)
  shared_percent[i] <- paste(shared_percent[i], "%)", sep="")
}
shared_percent
shared_otu_df$shared <- shared_percent

rownames(shared_otu_df) <- c("Epiphyte", "Endophyte", "Litter", "Soil")

out_df <- t(shared_otu_df)
out_df
rownames(out_df) <- c("Total OTUs", "Unique OTUs",
                      "Shared with Epiphytes", "Shared with Endophytes",
                      "Shared with Litter", "Total Shared")
out_df <- replace_na(out_df, "-")
write.csv(out_df, "./Tables/Shared_otus.csv")
