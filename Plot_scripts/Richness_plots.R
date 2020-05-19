library(ggplot2)
library(ggpubr)
library(vegan)
library(phyloseq)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))


phy_otu_table <- otu_table(rare_otu, taxa_are_rows = T) # Make phyloseq object from OTU table
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

ps_rich <- plot_richness(physeq1, x = "Substrate", color = "Plant_species", measures=c("Shannon", "InvSimpson"))
ps_rich <- ps_rich + geom_boxplot(alpha=0) +
  theme_pubr() +
  scale_x_discrete(breaks = c("Endo", "Epi", "Lit", "Soil"),
                   labels=c("Endophytes", "Epiphytes", "Litter", "Soil")) +
  scale_color_manual(values=c("#88CCEE","#332288"), guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic")))+
  labs(color = "Plant Species") +
  theme(legend.position = "right",
        axis.text.x.bottom = element_text(angle=90, vjust=0.5, hjust=1))
ps_rich
ggsave("./Figures/Richness_metrics_substrate.png", ps_rich, width=8, height = 6, units="in") # Save plot
ggsave("./Figures_Color/Supp Figure 3.pdf", ps_rich, width = 190, height = 160, units = "mm")
ggsave("./Figures_Numbered/Supp Figure 3.pdf", ps_rich, width = 190, height = 160, units = "mm")
