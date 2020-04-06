library(ggplot2)
library(ggpubr)
library(vegan)
library(patchwork)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communites")
sl <- c(1:9, 11, 13:21, 25:36, 38, 40:44, 48:55, 57:58) # Slice for non-negative or failed samples
otu_dat <- read.table("./Data/all_OTUS_R1_clean.txt", sep="\t", header=TRUE) # Read OTU table with contaminants removed
rownames(otu_dat) <- otu_dat[,1] # Rename rows with OTU number
otu_dat <- otu_dat[,order(colnames(otu_dat))] # Reorder by OTU number
colnames(otu_dat) <- c("OTU_ID", colnames(map)) # OTU_ID as first column, sample IDs as rest of columns
otu_dat_wo_negs <- as.matrix(otu_dat[,sl + 1]) # Keep non-negative samples for OTU table
map <- t(read.table("./Data/DEM_map_mod3.csv", fill = TRUE, sep=",", header=TRUE)) # Read the map file
colnames(map) <- map[1,] # Colnames are now by Sample ID
rownames(map) <- c("SampleID", rownames(map)[2:23]) # Rownames are now site/sample variables
map_wo_negs <- map[,sl] # Keep non-negative samples for map table
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
MDS_dat_df <- as.data.frame(MDS_points) # Convert to a df
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs)) # Add map data
NMDS_bray = data.frame(MDS1 = MDS_points[,1], # Make dataframe for plotting
                  MDS2 = MDS_points[,2],
                  group=MDS_dat_df$Soil_Leaf_Litter_Leaf_swab,
                  species=MDS_dat_df$Plant_species,
                  site=MDS_dat_df$Site)
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
  labs(alpha="Host species", color="Substrate", shape="Site") +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1)) +
  theme_pubr() +
  guides(fill=FALSE) +
  ggtitle("Bray-Curtis") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left") +
  scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
p_bc

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
p_jac <- ggplot(data = NMDS_jac, aes(MDS1, MDS2)) + # Make plot
  geom_point(aes(color = group, shape = site),size=3) +
  geom_point(aes(color = group, fill = group, alpha = species, shape=site),size=3) +
  geom_path(data=df_ell_jac, aes(x=NMDS1, y=NMDS2,color=group), size=1, linetype=2) +
  labs(alpha="Host species", color="Substrate", shape="Site") +
  scale_shape_manual(values = 21:25) +
  scale_alpha_manual(values=c(0,1)) +
  guides(fill=FALSE) +
  ggtitle("Jaccard") +
  theme_pubr() +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "right",
        legend.justification = "left") +
  scale_color_discrete(labels=c("Endophytes","Epiphytes","Litter", "Soil"))
p_jac

combined_NMDS <- p_bc + p_jac + plot_layout(guides="collect")
ggsave("Combined_NMDS_w_site.png", combined_NMDS, height=6, width=8, units="in")
