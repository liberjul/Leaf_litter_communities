# library(ggplot2)
library(vegan)
library(phyloseq)
library(tidyverse)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))
indicator_dat <- read.csv("./Stats/indic_species_5_per.csv") %>%
  filter(OTU != "OTU_91")

prop_otu_dat <- rare_otu # Assign the table to new df to maintain shape

# for (i in 1:dim(rare_otu)[2]){ # for each column
#   prop_otu_dat[,i] <- rare_otu[,i]/(sum(rare_otu[,i])*wts_wo_negs[i]) # Calculate proportions and scale by # of samples per substrate
# }
# prop_otu_dat_trim <- prop_otu_dat[rowSums(prop_otu_dat) > 4 * 0.001,] # Take OTUs with substrate-weighted abundance > 0.1%

prop_otu_dat_trim <- prop_otu_dat[rownames(prop_otu_dat) %in% indicator_dat$OTU,]
  

indicator_dat %>%
  arrange(index) %>%
  mutate(name = paste(OTU, taxonomic_id, sep=" ")) -> sorted_indic

rownames(prop_otu_dat_trim)


prop_otu_dat_trim <- prop_otu_dat_trim[order(rownames(prop_otu_dat_trim)),]  
rownames(prop_otu_dat_trim)

sorted_indic %>%
  arrange(OTU) %>%
  .[,"name"] -> tax_ids

tax_ids
rownames(prop_otu_dat_trim) <- tax_ids

indicator_dat$OTU

dim(prop_otu_dat_trim) # Check dimensions (for number of OTUs included)

phy_otu_table <- otu_table(prop_otu_dat_trim, taxa_are_rows = T) # Make phyloseq object from OTU table
#taxa_table <- as.matrix(read.csv("./Data/allrank_otus_R1.fasta_classified_unite.csv")) # Read UNITE classified taxa table, convert to matrix
#phy_tax_table <- tax_table(taxa_table) # Convert taxa table matrix to physeq object
#physeq <- phyloseq(phy_otu_table, phy_tax_table) # Make physeq object from OTU and taxa tables

physeq <- phyloseq(phy_otu_table)
otu_table(physeq)
phy_sample_table <- data.frame(Substrate = map_wo_negs["Substrate",], # Convert some variables to sample data table
                                           Sample = map_wo_negs["Short_name",],
                                           Plant_species = map_wo_negs["Plant_species",],
                                           Site = map_wo_negs["Site",],
                                           stringsAsFactors=FALSE)
phy_sample_table <- phy_sample_table %>%
  mutate(Name = paste(Site, # Make name variable from site and plant species
                      unlist(strsplit(Plant_species, " "))[c(T, F)],
                      sep = "-")) %>%
  arrange(Substrate, Plant_species, Site) # Sort by variables
rownames(phy_sample_table) <- sample_names(physeq) # Make sample names row names
phy_sample_table <- sample_data(phy_sample_table) # Convert to physeq sample table
phy_sample_table # Check table
physeq1 <- merge_phyloseq(physeq, phy_sample_table) # Merge sample data to physeq object
physeq1
# 
# RadialTheta <- function(pos){ # Function from https://github.com/joey711/phyloseq/blob/master/R/plot-methods.R to manually sort OTUs
#   pos = as(pos, "matrix")
#   xc  = mean(pos[, 1])
#   yc  = mean(pos[, 2])
#   theta = atan2((pos[, 2] - yc), (pos[, 1] - xc))
#   names(theta) <- rownames(pos)
#   return(theta)
# }

# ps.ord <- ordinate(physeq1, "NMDS", "bray") # Perform ordination
# specDF <- scores(ps.ord, choices = c(1, 2), display="species", physeq=physeq1) # Take scores by OTU
# taxa.order = taxa_names(physeq1)[order(RadialTheta(specDF))] # Sort OTUs

indicator_dat %>%
  arrange(index) %>%
  mutate(name = paste(OTU, taxonomic_id, sep=" "))-> sorted_indic
sorted_indic

ps_heatmap <- plot_heatmap(physeq1, method = NULL, distance = NULL, "Name",
                           taxa.order = rev(sorted_indic$name), # With manual OTU sorting
                           low="#66CCFF", high="#000033", na.value="white") # Make heatmap
facet_labs_sub <- c(Epi = "Epiphyte", Endo = "Endophyte", Lit = "Litter", Soil = "Soil") # Change facet labels

ps_heatmap

ps_heatmap <- ps_heatmap + facet_grid(~Substrate,
                                      switch = "x",
                                      scales = "free_x",
                                      space = "free_x",
                                      labeller = labeller(Substrate=facet_labs_sub)) + # Add facets for Substrate
  theme(axis.text.x.bottom = element_text(angle=90, vjust=0.5, hjust=1, size = 8)) + # Adjust the column labels
  labs(fill = "Proportion\n of Substrate") # Change x label
ps_heatmap$scales$scales[[1]]$name <- "Sample" # Change x axis name
ps_heatmap
ggsave("./Figures/physeq_heatmap.png", ps_heatmap, width=10, height = 8, units="in") # Save plot
ggsave("./Figures/physeq_heatmap.pdf", ps_heatmap, width=210, height = 150, units="mm") # Save plot
ggsave("./Figures_Color/Figure_1.pdf", ps_heatmap, width=210, height = 150, units="mm") # Save plot
ggsave("./Figures_Numbered/Figure_1.pdf", ps_heatmap, width=210, height = 150, units="mm")
ggsave("./Figures_Color/Figure_1.png", ps_heatmap, width=210, height = 150, units="mm") # Save plot
ggsave("./Figures_Color/Figure_1.eps", ps_heatmap, width=210, height = 150, units="mm") # Save plot


ps_heatmap <- plot_heatmap(physeq1, method = NULL, distance = NULL, "Name",
                           taxa.order = rev(sorted_indic$name), # With manual OTU sorting
                           low="white", high="black", na.value="white") # Make heatmap
facet_labs <- c(Epi = "Epiphyte", Endo = "Endophyte", Lit = "Litter", Soil = "Soil") # Change facet labels
ps_heatmap

ps_heatmap <- ps_heatmap + facet_grid(~Substrate,
                                      switch = "x",
                                      scales = "free_x",
                                      space = "free_x",
                                      labeller = labeller(Substrate=facet_labs)) + # Add facets for Substrate
  theme(axis.text.x.bottom = element_text(angle=90, vjust=0.5, hjust=1, size = 8)) + # Adjust the column labels
  labs(fill = "Proportion\n of Substrate") # Change x label
ps_heatmap$scales$scales[[1]]$name <- "Sample" # Change x axis name
ps_heatmap
ggsave("./Figures/physeq_heatmap_bw.png", ps_heatmap, width=12, height = 10, units="in") # Save plot
ggsave("./Figures/physeq_heatmap_bw.pdf", ps_heatmap, width=210, height = 150, units="mm") # Save plot
ggsave("./Figures/Figure_1.pdf", ps_heatmap, width=210, height = 150, units="mm") # Save plot
ggsave("./Figures/Figure_1.png", ps_heatmap, width=210, height = 150, units="mm") # Save plot
