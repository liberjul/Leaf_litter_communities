library(ggplot2)
library(vegan)
library(phyloseq)
library(ggpubr)

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

syn <- c(2, 4, 6, 8, 11, 13, 17, 19, 21) # Synthetic indexes
cot <- c(1, 3, 5, 7, 9, 14, 16, 18, 20) # Cotton indexes

syn_r <- which(map_wo_negs["Swab_type",] == "Synthetic")
cot_r <- which(map_wo_negs["Swab_type",] == "Cotton")
taxa_table <- as.matrix(read.csv("./Data/allrank_otus_R1.fasta_classified_unite.csv"))
phy_tax_table <- tax_table(taxa_table)

phy_swab <- phyloseq(otu_table(rare_otu[,union(cot_r, syn_r) + 1], taxa_are_rows = T),
                     phy_tax_table)

phy_sam_table <- sample_data(
  data.frame(
    row.names = sample_names(phy_swab),
    Swab_type=c(rep("Cotton", length(cot)),
                rep("Synthetic", length(syn))),
    stringsAsFactors=FALSE))
class(phy_sam_table)

phy_swab <- merge_phyloseq(phy_swab, phy_sam_table)
phy_swab

richness <- plot_richness(phy_swab, x="Swab_type", color="Swab_type")
richness <- richness +
  labs(x = "Swab Material", color = "Swab Material") +
  theme(legend.position = "none")
richness
ggsave("./Figures/Richness_metrics_swab_type.png", richness, width = 8, height = 6, units = "in")

swab_df <- data.frame(Swab_type=c(rep("Cotton", length(cot)),
                                  rep("Synthetic", length(syn))),
                      Read_count = c(colSums(otu_dat[cot+1]), colSums(otu_dat[syn+1])),
                      Leaf = c(1:5, 7:10, 1:6, 8:10))

swab_bp <- ggplot(swab_df, aes(x=Swab_type, y=Read_count, fill=Swab_type, pch=Swab_type)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x="Swab Material", y="Read Count") +
  theme_pubr() +
  theme(legend.position = "none")
swab_bp
ggsave("./Figures/swab_read_count_boxplot.png", swab_bp, width = 6, height = 6, units = "in")



veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)  # Make ellipses plottable
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

MDS_dat_swab <- metaMDS(t(rare_otu[,1:18]), distance = "bray")
MDS_points <- MDS_dat_swab$points
MDS_dat_df <- as.data.frame(MDS_points)
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs[,1:18]))
NMDS_swab = data.frame(MDS1 = MDS_points[,1],
                  MDS2 = MDS_points[,2],
                  group=MDS_dat_df$Swab_type,
                  species=MDS_dat_df$Plant_species,
                  site=MDS_dat_df$Site)
NMDS_swab.mean=aggregate(NMDS_swab[,1:2],list(group = MDS_dat_df$Soil_Leaf_Litter_Leaf_swab),mean)
plot.new()
ord<-ordiellipse(MDS_dat_swab, MDS_dat_df$Swab_type, display = "sites", kind = "se", conf = 0.97, label = T)
df_ell <- data.frame()
for(g in levels(MDS_dat_df$Swab_type)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(MDS_dat_df[MDS_dat_df$Swab_type==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
NMDS_swab.mean=aggregate(NMDS_swab[,1:2],list(group=NMDS_swab$group),mean)
p_bray_swab <- ggplot(data = NMDS_swab, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=group, shape=site), size=3) +
  geom_point(aes(color = group, fill = group, alpha = species, shape=site),size=3) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,color=group), size=1, linetype=2) +
  theme_pubr() +
  theme(plot.title = element_text(hjust=0.5), legend.position = "right") +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1), guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
  guides(fill=FALSE) +
  labs(shape="Site", color="Substrate", alpha="Host species")
p_bray_swab

ggsave("./Figures/bc_swab_NMDS_plot.png", p_bray_swab, height=6, width=8, units="in")
