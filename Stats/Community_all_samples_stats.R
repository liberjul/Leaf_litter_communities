library(ggplot2)
library(vegan)
library(RVAideMemoire)
library(indicspecies)
library(SpiecEasi)
library(huge)
library(MASS)
library(ggpubr)
library(randomForest)
library(rfUtilities)
library(compositions)
library(caret)
library(purrr)

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
write.table(as.data.frame(rare_otu), # Export rarified otu table for all analyses
            file = "./Data/Rare_otu_table.csv", na = "", sep=",", quote = F)
write.table(as.data.frame(map_wo_negs), # Export map without negatives
            file = "./Data/DEM_map_wo_negs.csv", na = "", sep=",", quote = F)

clr_otu <- clr(otu_dat_wo_negs)
clr_otu_df <- as.data.frame(clr_otu)
map_wo_negs_t <- t(map_wo_negs)
identical(colnames(clr_otu_df), rownames(map_wo_negs_t))
# Add classification variable and sample names
otus_clr_Comp <- data.frame(t(clr_otu_df))
otus_clr_Comp$Sample <- rownames(otus_clr_Comp)
otus_clr_Comp$Substrate <- as.factor(map_wo_negs_t[rownames(otus_clr_Comp), "Substrate"]) 
# otus_clr_Comp$Substrate <- as.numeric(
#   factor(
#     otus_clr_Comp$Substrate, labels = 1:4))
otus_clr_Comp$Substrate
head(otus_clr_Comp)

### TRain the RF model
set.seed(4563)
RF_501 <- randomForest(x=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                       y=otus_clr_Comp$Substrate,
                       ntree=501, 
                       #mtry = 20,
                       importance=TRUE,
                       proximity=TRUE)

str(RF_501)
plot(RF_501)

set.seed(3453)
RF_1001 <- randomForest(x=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                        y=otus_clr_Comp$Substrate,
                        ntree=1001, 
                        #mtry = 20,
                        importance=TRUE,
                        proximity=TRUE)

str(RF_1001)
plot(RF_1001)

importance(RF_501)
##Plot  Mean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_501) 
### Model Fit Analysis
oob_error_data_Comp <- data.frame(
  Trees=rep(1:nrow(RF_501$err.rate), times=5),
  Type=rep(c("OOB", "Endo", "Epi", "Lit", "Soil"), each=nrow(RF_501$err.rate)),
  Error=c(RF_501$err.rate[,"OOB"], 
          RF_501$err.rate[,"Endo"], 
          RF_501$err.rate[,"Epi"],
          RF_501$err.rate[,"Lit"],
          RF_501$err.rate[,"Soil"]))

plot_oob_error_Comp = ggplot(data=oob_error_data_Comp, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_bw() +
  annotate("text", x=250, y=0.75, label=paste("OOB estimate: ",
                                              round(RF_501$err.rate[nrow(RF_501$err.rate),1]*100,
                                                    digits=2), "%"), size=3) +
  geom_line(aes(color=Type)) + 
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 7)) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=1)) +
  theme(legend.position="bottom")

plot_oob_error_Comp
ggsave("oob_error_plot.png", plot_oob_error_Comp, width = 6, height = 5, units = "in")

print(RF_501)

sqrt(ncol(otus_clr_Comp))
oob_values <- c()
mtry_values <- seq(from = 50, to = 70, by = 1)
length(mtry_values)
length(oob_values)
oob_values

for(i in mtry_values) {
  print(i)
  temp.model <- randomForest(x=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                             y=as.factor(otus_clr_Comp$Substrate), 
                             mtry=i, ntree=501)
  oob_values <- c(oob_values, temp.model$err.rate[nrow(temp.model$err.rate),1])
  print(temp.model$err.rate[nrow(temp.model$err.rate),1])
}

oob_values
match(min(oob_values),oob_values)
mtry_values[match(min(oob_values),oob_values)]

set.seed(857)
RF_501_optm <- randomForest(x=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                            y=otus_clr_Comp$Substrate,
                            ntree=501, 
                            mtry = mtry_values[match(min(oob_values),oob_values)],
                            importance=TRUE,
                            proximity=TRUE)

RF_501_optm

RF_bact_sig_Comp <- rf.significance(x=RF_501, xdata=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                                    nperm=99,
                                    ntree=501)  
RF_bact_sig_Comp

RF_sig_Comp_optm <- rf.significance(x=RF_501_optm, xdata=otus_clr_Comp[,1:(ncol(otus_clr_Comp)-2)],
                                    nperm=1000,
                                    ntree=501)  
RF_sig_Comp_optm

### Hypothesis testing for all samples and groupings
sim_out <- simper(t(rare_otu), as.factor(map_wo_negs["Soil_Leaf_Litter_Leaf_swab",]), permutations = 100, parallel = 2)
dist_bray <- vegdist(t(rare_otu), method="bray")
dist_jac <- vegdist(t(rare_otu), method="jaccard")

bd_out <- betadisper(dist_bray, as.factor(map_wo_negs["Soil_Leaf_Litter_Leaf_swab",]))
model <- anova(bd_out)
boxplot(bd_out, xlab="Substrate")
TukeyHSD(bd_out)
summary(model)
model

MDS_dat <- metaMDS(t(rare_otu))
MDS_points <- MDS_dat$points
MDS_dat_df <- as.data.frame(MDS_points)
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs))

ado_bray <- adonis2(t(rare_otu) ~ Soil_Leaf_Litter_Leaf_swab + Plant_species + Site, MDS_dat_df, method="bray")
ado_bray
ado_jac <- adonis2(t(rare_otu) ~ Soil_Leaf_Litter_Leaf_swab + Plant_species + Site, MDS_dat_df, method="jaccard")
ado_jac

ado_bray_int <- adonis2(t(rare_otu) ~ Soil_Leaf_Litter_Leaf_swab * Plant_species * Site, MDS_dat_df, method="bray")
ado_bray_int
ado_jac_int <- adonis2(t(rare_otu) ~ Soil_Leaf_Litter_Leaf_swab * Plant_species * Site, MDS_dat_df, method="jaccard")
ado_jac_int
a <- pairwise.perm.manova(dist,MDS_dat_df$Soil_Leaf_Litter_Leaf_swab,nperm=999)
a$p.value


am_model <- anosim(dist, grouping=MDS_dat_df$Soil_Leaf_Litter_Leaf_swab,permutations=1000)
am_model
#norm_otu <- decostand(otu_dat_wo_negs, method = "log", MARGIN = 2)
#colSums(norm_otu)

library(stats)

otu_mp <- multipatt(as.data.frame(t(rare_otu)), map_wo_negs["Soil_Leaf_Litter_Leaf_swab",], control=how(nperm=9999))
summary(otu_mp, indvalcomp=TRUE)
write.table(as.data.frame(otu_mp$sign), file = "OTU_indicspecies.tsv", na = "", sep="\t", quote = F)
otu_mp -> otu_mp_fdr
otu_mp_fdr$sign$p.value <- p.adjust(otu_mp$sign$p.value, "fdr")
summary(otu_mp_fdr)

write.table(as.data.frame(otu_mp_fdr$sign), file = "OTU_indicspecies_fdr.tsv", na = "", sep="\t", quote = F)
otu_mp_fdr$sign
colnames(t(rare_otu))