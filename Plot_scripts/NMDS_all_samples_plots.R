library(ggplot2)
library(ggpubr)
library(vegan)
library(patchwork)

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
NMDS_bray.mean=aggregate(NMDS_bray[,1:2],list(group=NMDS_bray$group),mean) # Calculate mean for groups
p_bc <- ggplot(data = NMDS_bray, aes(MDS1, MDS2)) + # Make plot
  geom_point(aes(color = group, shape = site),size=3) +
  geom_point(aes(color = group, fill = group, alpha = species, shape=site),size=3) +
  geom_path(data=df_ell_bc, aes(x=NMDS1, y=NMDS2,color=group), size=1, linetype=2) +
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
           label = paste("Stress =", round(MDS_stress, 4))) +
  theme_pubr() +
  guides(fill=FALSE) +
  # ggtitle("Bray-Curtis") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left") +
  scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
p_bc

ggsave("./Figures/bc_NDMS_all_samples.png", p_bc, width = 8, height = 6, units="in")
MDS_stress
MDS_dat <- metaMDS(t(rare_otu), distance = "jaccard") # Calculate NMDS axes, using Jaccard distance
MDS_points <- MDS_dat$points # Extract coordinates
MDS_dat_df <- as.data.frame(MDS_points) # Convert to a df
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs)) # Add map data
NMDS_jac = data.frame(MDS1 = MDS_points[,1],  # Make dataframe for plotting
                      MDS2 = MDS_points[,2],
                      group=MDS_dat_df$Soil_Leaf_Litter_Leaf_swab,
                      species=MDS_dat_df$Plant_species,
                      site=MDS_dat_df$Site)
ord<-ordiellipse(MDS_dat, MDS_dat_df$Soil_Leaf_Litter_Leaf_swab, display = "sites", kind = "se", conf = 0.97, label = T) # Calculate ellipses
df_ell_jac <- data.frame()  # Dataframe for storing ellipses
for(g in levels(MDS_dat_df$Soil_Leaf_Litter_Leaf_swab)){
  df_ell_jac <- rbind(df_ell_jac, cbind(as.data.frame(with(MDS_dat_df[MDS_dat_df$Soil_Leaf_Litter_Leaf_swab==g,], # Add ellipse values
                                                           veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                        ,group=g))
}
NMDS_jac.mean=aggregate(NMDS_jac[,1:2],list(group=NMDS_jac$group),mean) # Calculate mean for groups
# p_jac <- ggplot(data = NMDS_jac, aes(MDS1, MDS2)) + # Make plot
#   geom_point(aes(color = group, shape = site),size=3) +
#   geom_point(aes(color = group, fill = group, alpha = species, shape=site),size=3) +
#   geom_path(data=df_ell_jac, aes(x=NMDS1, y=NMDS2,color=group), size=1, linetype=2) +
#   labs(alpha="Host species", color="Substrate", shape="Site") +
#   scale_shape_manual(values = 21:25) +
#   scale_alpha_manual(values=c(0,1), guide =
#                        guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
#   guides(fill=FALSE) +
#   ggtitle("Jaccard") +
#   theme_pubr() +
#   theme(plot.title = element_text(hjust=0.5),
#         legend.position = "right",
#         legend.justification = "left") +
#   scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
# p_jac

# combined_NMDS <- p_bc + p_jac + plot_layout(guides="collect")
# ggsave("./Figures/Combined_NMDS_w_site.png", combined_NMDS, height=6, width=8, units="in")

bc.d = vegdist(t(rare_otu), method="bray")
sor.d = vegdist(t(rare_otu), method = "bray", binary = T)
jac.d = vegdist(t(rare_otu), method="jaccard")

bc.pcoa = cmdscale(bc.d, eig=T)
ax1.v.bc = bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc = bc.pcoa$eig[2]/sum(bc.pcoa$eig)

sor.pcoa = cmdscale(sor.d, eig=T)
ax1.v.sor = sor.pcoa$eig[1]/sum(sor.pcoa$eig)
ax2.v.sor = sor.pcoa$eig[2]/sum(sor.pcoa$eig)

jac.pcoa = cmdscale(jac.d, eig=T)
ax1.v.jac = jac.pcoa$eig[1]/sum(jac.pcoa$eig)
ax2.v.jac = jac.pcoa$eig[2]/sum(jac.pcoa$eig)

ax1.v.bc
ax1.v.sor
ax1.v.jac

ax2.v.bc
ax2.v.sor
ax2.v.jac

ax1.v.bc + ax2.v.bc
ax1.v.sor + ax2.v.sor
ax1.v.jac + ax2.v.jac

bc.pcoa_df <- data.frame(ax1 = bc.pcoa$points[,1],
                         ax2 = bc.pcoa$points[,2],
                         substrate=MDS_dat_df$Substrate,
                         species=MDS_dat_df$Plant_species,
                         site=MDS_dat_df$Site)

pcoa_bc <- ggplot(data = bc.pcoa_df, aes(ax1, ax2)) + # Make plot
  geom_point(aes(color = substrate, shape = site),size=3) +
  geom_point(aes(color = substrate, fill = substrate, alpha = species, shape=site),size=3) +
  labs(alpha="Host species", color="Substrate", shape="Site",
       x = paste("PCoA1: ",100*round(ax1.v.bc,3),"% var. explained",sep=""),
       y = paste("PCoA2: ",100*round(ax2.v.bc,3),"%var. explained",sep="")) +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1), guide =
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"))) +
  guides(fill=FALSE) +
  theme_pubr() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left") +
  scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
pcoa_bc
ggsave("./Figures/bc_pcoa_all_samples.png", pcoa_bc, width = 8, height=6, units="in")

ado_bray_pcoa <- adonis2(t(rare_otu) ~ substrate * species * site,
                         bc.pcoa_df, method="bray")
ado_bray_pcoa
