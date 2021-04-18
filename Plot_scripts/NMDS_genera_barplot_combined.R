library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(vegan)
library(patchwork)
library(magick)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))
colnames(map_wo_negs) <- map_wo_negs["SampleID",]
colnames(rare_otu) <- map_wo_negs["SampleID",]
# taxa_table <- read.delim("./Data/consensus_taxonomy_constax.txt") # Load CONSTAX classifications
taxa_table <- read.delim("./Data/constax_taxonomy.txt") # Load CONSTAX classifications
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
  dplyr::select(Genus, Sample, value) %>%
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
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.text = element_text(size = 10),
        text = element_text(size = 10),
        legend.key.size = unit(.5, "line"),
        legend.margin = margin(t = 30, unit = "pt"),
        legend.spacing.y = unit(5.0, 'pt')) +
  guides(fill = guide_legend(reverse = TRUE, ncol =1)) # Reverse legend to match

#genus_barplot # Show plot

# ggsave("./Figures/top30_genera_gg_constax.png", genus_barplot, width=14, height = 8, units="in") # Save plot

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)  # Make ellipses plottable
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

MDS_dat <- metaMDS(t(rare_otu)) # Calculate NMDS axes, using Bray-Curtis as default
MDS_points <- MDS_dat$points # Extract coordinates
MDS_stress <- MDS_dat$stress
MDS_dat_df <- as.data.frame(MDS_points) # Convert to a df
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs)) # Add map data
NMDS_bray = data.frame(MDS1 = MDS_points[,1], # Make dataframe for plotting
                       MDS2 = MDS_points[,2],
                       group=MDS_dat_df$Soil_Leaf_Litter_Leaf_swab,
                       species=MDS_dat_df$Plant_species,
                       site=MDS_dat_df$Site)
plot.new()
ord<-ordiellipse(MDS_dat, MDS_dat_df$Soil_Leaf_Litter_Leaf_swab, display = "sites", kind = "se", conf = 0.97, label = T) # Calculate ellipses
df_ell_bc <- data.frame() # Dataframe for storing ellipses
for(g in levels(MDS_dat_df$Soil_Leaf_Litter_Leaf_swab)){
  df_ell_bc <- rbind(df_ell_bc, cbind(as.data.frame(with(MDS_dat_df[MDS_dat_df$Soil_Leaf_Litter_Leaf_swab==g,], # Add ellipse values
                                                         veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                      ,group=g))
}
NMDS_bray.mean<-aggregate(NMDS_bray[,1:2],  # Calculate mean for groups
                          list(group=NMDS_bray$group),mean)

NMDS_bray.mean

scrs <- cbind(as.data.frame(MDS_points), group=NMDS_bray$group)
segs <- merge(scrs, setNames(NMDS_bray.mean,
                             c("group", "MDS1", "MDS2")),
              by = "group", sort = FALSE)

p_bc <- ggplot(data = NMDS_bray, aes(MDS1, MDS2)) + # Make plot
  geom_point(aes(color = group, shape = site),size = 3) +
  geom_point(aes(color = group, fill = group, alpha = species, shape=site),size = 3) +
  geom_path(data=df_ell_bc, aes(x=NMDS1, y=NMDS2,color=group), size=0.5, linetype=2) +
  geom_segment(data = segs,
               mapping = aes(x = MDS1.x, y = MDS2.x,
                             xend = MDS1.y, yend = MDS2.y,
                             color = group)) +
  geom_text(data = NMDS_bray.mean,
            aes(x = MDS1, y = MDS2, label = c("Endophytes","Epiphytes","Litter", "Soil")),
            size = 10*(5/14)) +
  labs(alpha="Host species", color="Substrate", shape="Site",
       x = "NMDS1", y = "NMDS2") +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1), guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"),
                                    override.aes = list(pch = 21,
                                                        color = 1,
                                                        alpha = 1,
                                                        fill = c(NA, 1)))) +
  annotate(geom = "text", hjust = 0,
           x = min(NMDS_bray$MDS1), y = min(NMDS_bray$MDS2),
           label = paste("Stress =", round(MDS_stress, 4)),
           size = 10*(5/14)) +
  theme_pubr() +
  guides(fill=FALSE) +
  # ggtitle("Bray-Curtis") +
  # theme(plot.title = element_text(hjust=0.5),
  #       legend.position = "right",
  #       legend.justification = "left",
  #       legend.text = element_text(size=7),
  #       legend.spacing = unit(0, "pt"),
  #       legend.key.height = unit(10, "pt"),
  #       text = element_text(size = 10),
  #       legend.margin=margin(t = 0, unit='pt')) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left",
        legend.text = element_text(size = 10),
        legend.spacing = unit(0, "pt"),
        legend.key.height = unit(12, "pt"),
        text = element_text(size = 10),
        legend.margin=margin(t = 0, unit='pt')) +
  scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
p_bc

# ggsave("./Figures/bc_NDMS_all_samples.png", p_bc, width = 8, height = 6, units="in")

MDS_dat_swab <- metaMDS(t(rare_otu[,1:18]), distance = "bray")
MDS_points <- MDS_dat_swab$points
MDS_stress <- MDS_dat_swab$stress
MDS_dat_df <- as.data.frame(MDS_points)
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs[,1:18]))
NMDS_swab = data.frame(MDS1 = MDS_points[,1],
                       MDS2 = MDS_points[,2],
                       group=MDS_dat_df$Swab_type,
                       species=MDS_dat_df$Plant_species,
                       site=MDS_dat_df$Site)
plot.new()
ord<-ordiellipse(MDS_dat_swab, MDS_dat_df$Swab_type, display = "sites", kind = "se", conf = 0.97, label = T)
df_ell <- data.frame()
for(g in levels(MDS_dat_df$Swab_type)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(MDS_dat_df[MDS_dat_df$Swab_type==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
NMDS_swab.mean=aggregate(NMDS_swab[,1:2],list(group=NMDS_swab$group),mean)

scrs <- cbind(as.data.frame(MDS_points), group=NMDS_swab$group)
segs <- merge(scrs, setNames(NMDS_swab.mean,
                             c("group", "MDS1", "MDS2")),
              by = "group", sort = FALSE)

p_bray_swab <- ggplot(data = NMDS_swab, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=group, shape=site), size = 3) +
  geom_point(aes(color = group, fill = group, alpha = species, shape=site),size = 3) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,color=group), size=0.5, linetype=2) +
  geom_segment(data = segs,
               mapping = aes(x = MDS1.x, y = MDS2.x,
                             xend = MDS1.y, yend = MDS2.y,
                             color = group),
               size = 0.1) +
  geom_text(data = NMDS_swab.mean,
            aes(x = MDS1, y = MDS2, label = group),
            size = 10*(5/14)) +
  theme_pubr() +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1), guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"),
                                    override.aes = list(pch = 21,
                                                        color = 1,
                                                        alpha = 1,
                                                        fill = c(NA, 1)),
                                    order=1)) +
  guides(color = guide_legend(order=2),
         shape = guide_legend(order=3)) +
  annotate(geom = "text", hjust = 0,
           x = min(NMDS_swab$MDS1), y = min(NMDS_swab$MDS2),
           label = paste("Stress =", round(MDS_stress, 4)),
           size = 10*(5/14)) +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left",
        legend.text = element_text(size = 10),
        legend.spacing = unit(0, "pt"),
        legend.key.height = unit(12, "pt"),
        text = element_text(size = 10),
        legend.margin=margin(t = 0, unit='pt')) +
  guides(fill=FALSE) +
  labs(shape="Site", color="Swab Material", alpha="Host Species",
       x = "NMDS1", y = "NMDS2")
p_bray_swab

g <- (genus_barplot + labs(tag = "A")) / ((p_bc + labs(tag = "B")) +
                                            (p_bray_swab + labs(tag = "C")))
# g
ggsave("./Figures/NMDS_genera_barplot_combined.png", g, width = 15, height = 13, units = "in")
ggsave("./Figures/NMDS_genera_barplot_combined.eps", g, width = 190, height = 190,
       units = "mm", dpi=1000)
ggsave("./Figures_Color/Figure 5.eps", g, width = 190, height = 190,
       units = "mm", dpi=1000)
ggsave("./Figures_Color/Figure 5.eps", g, width = 15, height = 13, units = "in")
ggsave("./Figures_Color/Figure 5.pdf", g, width = 15, height = 13, units = "in")
