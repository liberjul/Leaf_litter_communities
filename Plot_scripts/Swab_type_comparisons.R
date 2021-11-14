library(ggplot2)
library(vegan)
library(phyloseq)
library(ggpubr)
library(tidyverse)
library(lme4)
library(lmerTest)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))
colnames(map_wo_negs) <- map_wo_negs["SampleID",]
colnames(rare_otu) <- map_wo_negs["SampleID",]
otu_dat <- read.table("./Data/all_OTUS_R1_clean.txt", sep="\t", header=TRUE) # Read OTU table with contaminants removed
rownames(otu_dat) <- otu_dat[,1] # Rename rows with OTU number
otu_dat <- otu_dat[,order(colnames(otu_dat))] # Reorder by OTU number

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
    Leaf = c(1:5, 7:10, 1:6, 8:10),
    stringsAsFactors=FALSE))
class(phy_sam_table)

phy_swab <- merge_phyloseq(phy_swab, phy_sam_table)
phy_swab

diver_df <- cbind(estimate_richness(phy_swab, measures = c("Shannon", "InvSimpson", "Observed")),
                  sample_data(phy_swab))

diver_df                       
m1 <- lmer(Shannon^2 ~ Swab_type + (1|Leaf), diver_df, REML=FALSE)
summary(m1) -> sum_m1
sum_m1
abs(confint(m1))^0.5 * sign(confint(m1))
anova(m1) -> aov_m1
aov_m1$`Pr(>F)`

m2 <- lmer(sqrt(InvSimpson) ~ Swab_type + (1|Leaf), diver_df, REML=FALSE)
summary(m2) -> sum_m2
sum_m2
confint(m2)^2 * sign(confint(m2))
anova(m2) -> aov_m2

m3 <- lmer(Observed^2 ~ Swab_type + (1|Leaf), diver_df, REML=FALSE)
summary(m3) -> sum_m3
sum_m3
abs(confint(m3))^0.5 * sign(confint(m3))
anova(m3) -> aov_m3

swab_richness_stats <- data.frame(variable = c("Shannon", "InvSimpson", "Observed"),
                                  value = c(4.5, 37, 370),
                                  text = c(paste0("F = ", round(aov_m1$`F value`, 2),
                                                  ", p = ", round(aov_m1$`Pr(>F)`, 2)),
                                           paste0("F = ", round(aov_m2$`F value`, 2),
                                                  ", p = ", round(aov_m2$`Pr(>F)`, 2)),
                                           paste0("F = ", round(aov_m3$`F value`, 2),
                                                  ", p = ", round(aov_m3$`Pr(>F)`, 2))),
                                  Swab_type = "Cotton")

richness <- plot_richness(phy_swab, x="Swab_type", color="Swab_type", measures=c("Observed", "Shannon", "InvSimpson"))
richness <- richness + geom_boxplot(alpha=0) + #geom_jitter() +
  labs(x = "Swab Material", color = "Swab Material") +
  # scale_x_discrete(breaks = c("Cotton", "Synthetic"),
  #                  labels = c("C", "S")) +
  theme_pubr() +
  theme(legend.position = "none") +
  geom_text(data=swab_richness_stats,
            aes(y=value, x = Swab_type, label=text), color = "black", hjust=0)
richness
ggsave("./Figures/Richness_metrics_swab_type.png", richness, width = 8, height = 6, units = "in")
ggsave("./Figures_Numbered/Figure_4.pdf",
       richness,
       width = 8, height = 6, units = "in")
ggsave("./Figures_Color/Figure_4.pdf",
       richness + scale_color_manual(values = c("#88CCEE","#332288")),
       width = 90, height = 60, units = "mm")
ggsave("./Figures_Color/Figure_4.eps",
       richness + scale_color_manual(values = c("#88CCEE","#332288")),
       width = 90, height = 60, units = "mm")

swab_df <- data.frame(Swab_type=c(rep("Cotton", length(cot)),
                                  rep("Synthetic", length(syn))),
                      Read_count = c(colSums(otu_dat[cot+1]), colSums(otu_dat[syn+1])),
                      Leaf = c(1:5, 7:10, 1:6, 8:10))

swab_df %>%
  group_by(Swab_type) %>%
  summarise(mean = mean(Read_count),
            sd = sd(Read_count)) -> swab_df_sum
swab_bp <- ggplot(swab_df,
                  aes(x=Swab_type, y=Read_count, color=Swab_type)) +
  geom_boxplot() +
  geom_point(size = 2, alpha = 0.3) +
  geom_point(data = swab_df_sum,
             aes(x = Swab_type,
                 y = mean),
             size = 3) +
  labs(x="Swab Material", y="Read Count") +
  theme_pubr() +
  theme(legend.position = "none")
swab_bp
ggsave("./Figures/swab_read_count_boxplot.png", swab_bp, width = 6, height = 4, units = "in")
ggsave("./Figures_Numbered/Figure_S5.pdf", swab_bp, width = 6, height = 4, units = "in")
ggsave("./Figures_Color/Figure_S5.pdf", swab_bp, width = 6, height = 4, units = "in")
ggsave("./Figures/Figure_S5.pdf", swab_bp, width = 6, height = 4, units = "in")



veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)  # Make ellipses plottable
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

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
                       guide_legend(label.theme = element_text(size = 10, angle = 0, face = "italic"),
                                    override.aes = list(pch = 21,
                                                        color = 1,
                                                        alpha = 1,
                                                        fill = c(NA, 1)))) +
  annotate(geom = "text", hjust = 0,
           x = min(NMDS_swab$MDS1), y = min(NMDS_swab$MDS2),
           label = paste("Stress =", round(MDS_stress, 4))) +
  guides(fill=FALSE) +
  labs(shape="Site", color="Swab Material", alpha="Host Species",
       x = "NMDS1", y = "NMDS2")
p_bray_swab

MDS_stress
ggsave("./Figures/bc_swab_NMDS_plot.png", p_bray_swab, height=6, width=8, units="in")
