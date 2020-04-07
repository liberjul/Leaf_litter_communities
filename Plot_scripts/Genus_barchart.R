library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
sl <- c(1:9, 11, 13:14, 16:21, 25:36, 38, 40:44, 48:55, 57:58) # Slice for non-negative or failed samples
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
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

taxa_table <- read.delim("./Data/consensus_taxonomy_constax.txt") # Load CONSTAX classifications
rownames(taxa_table) <- taxa_table$OTU_ID # Rename rows with the OTU_ID

prop_otu_dat <- rare_otu # Assign the table to new df to maintain shape

for (i in 1:dim(rare_otu)[2]){ # for each column
  prop_otu_dat[,i] <- rare_otu[,i]/(sum(rare_otu[,i])*wts_wo_negs[i]) # Calculate proportions and scale by # of samples per substrate
}
prop_otu_dat <- as.data.frame(prop_otu_dat) # Convert matrix to df
prop_otu_dat$OTU_name <- rownames(rare_otu) # Same rownames as OTU table
taxa_tab_trim <- taxa_table[rownames(rare_otu),] # Reorder taxa table by otu table order
prop_otu_dat <- cbind(prop_otu_dat, taxa_tab_trim) # Combine dataframes

prop_otu_long <- prop_otu_dat %>% # Pivot by sample name to make long
  pivot_longer(cols=starts_with("2"), names_to="Sample")

# prop_otu_long # check df
genus_tab <- prop_otu_long %>% group_by(Genus) %>% summarize(wt_read_count = sum(value)) # Create df grouped by genus
genus_tab <- arrange(genus_tab, desc(wt_read_count)) # Reorder based on abundance
genus_top <- genus_tab %>% top_n(30) # Take top 30 genera, using weighted abundance/read count
if ("" %in% genus_top$Genus){ # If unnamed is in the top 30, drop it
  genus_top <- genus_tab %>% top_n(31)
  genus_top <- genus_top[genus_top$Genus != "",]
}
# genus_top # check table

prop_long_gen_filt <- prop_otu_long %>% # take proportions for samples in top 30 genera
  select(Genus, Sample, value) %>%
  filter(Genus %in% genus_top$Genus)
# prop_long_gen_filt # Check prop table

other_tab <- prop_otu_long %>% # Table for the sum of the other genera
  filter(!(Genus %in% genus_top$Genus)) %>%
  group_by(Sample) %>%
  summarize(value = sum(value))

other_tab$Genus <- rep("Other", dim(other_tab)[1]) # Set name to other

prop_long_gen_filt <- rbind(prop_long_gen_filt, other_tab) # Add "Other" table to proportion table
palette_CB9 <- c('#5e0215', '#1d0f6f', '#03149c', '#761c47', '#08926d', '#557741', # Custom color pallete for 31 samples
                 '#3c16c2', '#98433c', '#c64f11', '#70a228', '#bf8c01', '#1788b3',
                 '#79de00', '#8e15b8', '#4c64ae', '#5ce331', '#49e26e', '#fe159a',
                 '#c6955d', '#fd8355', '#63b3e1', '#ebd34d', '#a96bfe', '#a3f08d',
                 '#c0adb4', '#8dd2c8', '#d8e587', '#d59eed', '#fe71f4', '#b1dcf0',
                 '#d8e9c7')

genus_barplot <- ggplot(prop_long_gen_filt, # Ggplot barplot of proportions
                        aes(x = Sample, y = value, fill = Genus)) +
  geom_bar(position = position_fill(reverse=TRUE), stat="identity") + # Stacked bars
  theme_pubr() +
  scale_fill_manual(values= palette_CB9) + # Use palette
  scale_x_discrete(limits = map_wo_negs["SampleID",], # Replace sample names with shortened sample descriptions
                   labels = map_wo_negs["Short_name",]) +
  ylim(0,1) +
  labs(y = "Proportion of Reads Per Sample", fill = "Genus") + # Change variable labels
  theme(legend.position = "right", # Legend, ticks, and text adjustments
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  guides(fill = guide_legend(reverse = TRUE)) # Reverse legend to match

# genus_barplot # Show plot

ggsave("./Figures/top30_genera_gg_constax.png", genus_barplot, width=14, height = 8, units="in") # Save plot
