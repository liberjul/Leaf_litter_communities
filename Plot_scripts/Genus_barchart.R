library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))
colnames(map_wo_negs) <- map_wo_negs["SampleID",]
colnames(rare_otu) <- map_wo_negs["SampleID",]
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
colnames(prop_otu_dat)
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
colnames(prop_otu_long)
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
subs <- vector(length = dim(prop_long_gen_filt)[1])
sites <- subs
hosts <- subs
short <- subs
for (i in 1:dim(prop_long_gen_filt)[1]){
  subs[i] <- map_wo_negs["Substrate", prop_long_gen_filt$Sample[i]]
  sites[i] <- map_wo_negs["Site", prop_long_gen_filt$Sample[i]]
  hosts[i] <- map_wo_negs["Plant_species", prop_long_gen_filt$Sample[i]]
  short[i] <- map_wo_negs["Short_name", prop_long_gen_filt$Sample[i]]
}
prop_long_gen_filt <- prop_long_gen_filt %>%
  mutate(Substrate = subs,
         Site = sites,
         Host = hosts,
         Name = paste(Site,
                      unlist(strsplit(Host, " "))[c(T, F)],
                      sep = "-")) %>%
  arrange(Substrate, Host, Site)

lookup <- tibble(Sample = unique(prop_long_gen_filt$Sample)) %>%
  mutate(ind = row_number())
lookup
prop_long_gen_filt <- prop_long_gen_filt %>%
  inner_join(., lookup)

palette_CB9 <- c('#5e0215', '#1d0f6f', '#03149c', '#761c47', '#08926d', '#557741', # Custom color pallete for 31 samples
                 '#3c16c2', '#98433c', '#c64f11', '#70a228', '#bf8c01', '#1788b3',
                 '#79de00', '#8e15b8', '#4c64ae', '#5ce331', '#49e26e', '#fe159a',
                 '#c6955d', '#fd8355', '#63b3e1', '#ebd34d', '#a96bfe', '#a3f08d',
                 '#c0adb4', '#8dd2c8', '#d8e587', '#d59eed', '#fe71f4', '#b1dcf0',
                 '#d8e9c7')
colnames(prop_long_gen_filt)
facet_labs <- c(Epi = "Epiphyte", Endo = "Endophyte", Lit = "Litter", Soil = "Soil") # Change facet labels
prop_long_gen_filt
genus_barplot <- ggplot(prop_long_gen_filt, # Ggplot barplot of proportions
                        aes(x = reorder(Sample, ind), y = value, fill = Genus)) +
  geom_bar(position = position_fill(reverse=TRUE), stat="identity") + # Stacked bars
  scale_x_discrete(#limits = prop_long_gen_filt$Sample,
                   breaks = prop_long_gen_filt$Sample, # Replace sample names with shortened sample descriptions
                   labels = prop_long_gen_filt$Name) +
  facet_grid(~Substrate,
            switch = "x",
            scales = "free_x",
            space = "free_x",
            labeller = labeller(Substrate=facet_labs)) +
  theme_pubr() +
  scale_fill_manual(values= palette_CB9) + # Use palette
  
  ylim(0,1) +
  labs(y = "Proportion of Reads Per Sample", fill = "Genus", x = "Sample") + # Change variable labels
  theme(legend.position = "right", # Legend, ticks, and text adjustments
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  guides(fill = guide_legend(reverse = TRUE)) # Reverse legend to match

genus_barplot # Show plot

ggsave("./Figures/top30_genera_gg_constax.png", genus_barplot, width=14, height = 8, units="in") # Save plot


### Specific genera proportions
# prop_long_gen_filt %>%
#   filter(Substrate == "Endo",
#          Host == "Acer rubrum",
#          Genus == "Ramularia") %>%
#   summarize(sum(value))
# 
# prop_otu_long %>%
#   filter(Substrate == "Endo",
#          Host == "Acer rubrum",
#          Genus == "Pestalotiopsis") %>%
#   summarize(sum(value))
# 
# acer_endo_by_gen <- prop_otu_long %>%
#   filter(Sample %in% c(249:252, 254),
#          Genus %in% c("Phyllosticta", "Ramularia",
#                       "Colletotrichum", "Plagiostoma", "Pestalotiopsis")) %>%
#   group_by(Genus) %>%
#   summarize(sum(value))
# acer_endo_by_gen
# 
# acer_endo_by_gen <- prop_otu_long %>%
#   filter(Sample %in% c(244, 245, 247),
#          Genus %in% c("Colletotrichum")) %>%
#   group_by(Genus) %>%
#   summarize(sum(value))
# acer_endo_by_gen
# 
# carya_sum <- prop_long_gen_filt %>%
#   filter(Substrate == "Epi",
#          Host == "Carya ovata")%>%
#   .[,"value"] %>%
#   sum()
# carya_sum
# 
# acer_sum <- prop_long_gen_filt %>%
#   filter(Substrate == "Epi",
#          Host == "Acer rubrum")%>%
#   .[,"value"] %>%
#   sum()
# acer_sum
# 
# prop_long_gen_filt %>%
#   filter(Substrate == "Epi",
#          Host == "Carya ovata") %>%
#   group_by(Genus) %>%
#   summarize(gen_sum = sum(value) / carya_sum *100) %>%
#   arrange(-gen_sum)
# 
# prop_long_gen_filt %>%
#   filter(Substrate == "Epi",
#          Host == "Acer rubrum",
#          Genus == c("Erysiphe", "Golubevia")) %>%
#   group_by(Genus) %>%
#   summarize(gen_sum = sum(value)/ acer_sum *100) %>%
#   arrange(-gen_sum)
# 
#   
# 
# prop_otu_long %>%
#   filter(Sample %in% c(244, 245, 247),
#          Genus %in% c("Erysiphe")) %>%
#   summarize(sum(value)) /
#   
#   prop_otu_long %>%
#   filter(Sample %in% c(244, 245, 247)) %>%
#   summarize(sum(value))
# 
# 
# prop_otu_long %>%
#   filter(Sample %in% c(249:252, 254),
#          Genus %in% c("Exobasidium")) %>%
#   group_by(Genus) %>%
#   summarize(sum(value)) /
#   
#   prop_otu_long %>%
#   filter(Sample %in% c(249:252, 254)) %>%
#   summarize(sum(value))
# 
# 
# 
# prop_otu_long %>%
#   filter(#Substrate == "Soil",
#          # Host == "Carya ovata",
#          Genus %in% c("Phyllosticta", "Ramularia",
#                       "Colletotrichum", "Plagiostoma", "Pestalotiopsis")) %>%
#   group_by(Genus) %>%
#   summarize(sum(value))
# 
# prop_otu_long %>%
#   filter(Sample %in% c(249:252, 254))
