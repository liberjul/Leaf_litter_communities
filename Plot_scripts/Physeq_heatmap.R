library(ggplot2)
library(vegan)
library(phyloseq)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))

prop_otu_dat <- rare_otu # Assign the table to new df to maintain shape

for (i in 1:dim(rare_otu)[2]){ # for each column
  prop_otu_dat[,i] <- rare_otu[,i]/(sum(rare_otu[,i])*wts_wo_negs[i]) # Calculate proportions and scale by # of samples per substrate
}
prop_otu_dat_trim <- prop_otu_dat[rowSums(prop_otu_dat) > 4 * 0.001,]

phy_otu_table <- otu_table(prop_otu_dat_trim, taxa_are_rows = T) # Make phyloseq object from OTU table
taxa_table <- as.matrix(read.csv("./Data/allrank_otus_R1.fasta_classified_unite.csv")) # Read UNITE classified taxa table, convert to matrix
phy_tax_table <- tax_table(taxa_table) # Convert taxa table matrix to physeq object
physeq <- phyloseq(phy_otu_table, phy_tax_table) # Make physeq object from OTU and taxa tables
phy_sample_table <- sample_data(data.frame(Substrate = map_wo_negs["Substrate",], # Convert some variable to sample data table
                               Sample = map_wo_negs["Short_name",],
                               Plant_species = map_wo_negs["Plant_species",],
                               Site = map_wo_negs["Site",],
                               row.names = sample_names(physeq),
                               stringsAsFactors=FALSE))

physeq1 <- merge_phyloseq(physeq, phy_sample_table) # Merge sample data to physeq object
ps_heatmap <- plot_heatmap(physeq1, "NMDS", "bray", "Sample") # Make heatmap
facet_labs <- c(Epi = "Epiphyte", Endo = "Endophyte", Lit = "Litter", Soil = "Soil") # Change facet labels
ps_heatmap <- ps_heatmap + facet_grid(~Substrate,
                                      switch = "x",
                                      scales = "free_x",
                                      space = "free_x",
                                      labeller = labeller(Substrate=facet_labs)) + # Add facets for Substrate
  theme(axis.text.x.bottom = element_text(angle=90, vjust=0.5, hjust=1)) + # Adjust the column labels
  labs(x = "Sample Name", fill = "Proportion\n of Substrate") # Change x label
ps_heatmap
ggsave("./Figures/physeq_heatmap.png", ps_heatmap, width=12, height = 8, units="in") # Save plot
