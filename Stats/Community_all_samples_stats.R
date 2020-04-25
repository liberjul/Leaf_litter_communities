library(ggplot2)
library(vegan)
library(RVAideMemoire)
library(indicspecies)
library(ggpubr)
library(bbmle)
library(randomForest)
library(rfUtilities)
library(compositions)
library(caret)
library(purrr)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
wts_wo_negs <- c(rep(18,18), rep(10,10), rep(8, 8), rep(10, 10)) # Vector for weighting abundance by number of samples per substrate
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

#################
# ONLY RUN ONCE #
#################
# rare_otu <- t(rrarefy(t(otu_dat_wo_negs), # Rarefy by the lowest read count in non-negative sample
#                       min(colSums(otu_dat_wo_negs))))
# write.table(as.data.frame(rare_otu), # Export rarified otu table for all analyses
#             file = "./Data/Rare_otu_table.csv", na = "", sep=",", quote = F)
# write.table(as.data.frame(map_wo_negs), # Export map without negatives
#             file = "./Data/DEM_map_wo_negs.csv", na = "", sep=",", quote = F)

map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))

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
varImpPlot(RF_501) # From caret
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
sim_out <- simper(t(rare_otu), as.factor(map_wo_negs["Substrate",]), permutations = 100, parallel = 2)
dist_bray <- vegdist(t(rare_otu), method="bray")
dist_jac <- vegdist(t(rare_otu), method="jaccard")

bc.pcoa <- cmdscale(dist_bray, eig=T)
ax1.v.bc=bc.pcoa$eig[1]/sum(bc.pcoa$eig)
ax2.v.bc=bc.pcoa$eig[2]/sum(bc.pcoa$eig)
ax2.v.bc

bd_out <- betadisper(dist_bray, as.factor(map_wo_negs["Substrate",]))
model <- anova(bd_out)
model

bd_out <- betadisper(dist_bray, as.factor(map_wo_negs["Plant_species",]))
model <- anova(bd_out)
model

boxplot(bd_out, xlab="Substrate")
TukeyHSD(bd_out)


MDS_dat <- metaMDS(t(rare_otu))
MDS_points <- MDS_dat$points
MDS_dat_df <- as.data.frame(MDS_points)
MDS_dat_df <- cbind(MDS_dat_df, t(map_wo_negs))
round(MDS_dat$eig*100/sum(MDS_dat$eig),1)
### Code from: https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r
n = 10
stress <- vector(length = n)
for (i in 1:n) {
  stress[i] <- metaMDS(t(rare_otu), distance = "jaccard", k = i)$stress
}
names(stress) <- paste0(1:n, "Dim")
# x11(width = 10/2.54, height = 7/2.54)
par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")
###

ado_bray <- adonis2(t(rare_otu) ~ Substrate + Plant_species + Site,
                    MDS_dat_df, method="bray")
ado_bray
write.csv(ado_bray, "./Stats/ado_bray_stats_table.csv")
ado_jac <- adonis2(t(rare_otu) ~ Substrate + Plant_species + Site,
                   MDS_dat_df, method="jaccard")
ado_jac
write.csv(ado_jac, "./Stats/ado_jac_stats_table.csv")
ado_bray_int <- adonis2(t(rare_otu) ~ Substrate * Plant_species * Site,
                        MDS_dat_df, method="bray")
ado_bray_int
write.csv(ado_bray_int, "./Stats/ado_bray_int_stats_table.csv")
ado_jac_int <- adonis2(t(rare_otu) ~ Substrate * Plant_species * Site,
                       MDS_dat_df, method="jaccard")
ado_jac_int
write.csv(ado_jac_int, "./Stats/ado_jac_int_stats_table.csv")

a <- pairwise.perm.manova(dist_bray,MDS_dat_df$Substrate,nperm=999)
a


am_model <- anosim(dist_bray, grouping=MDS_dat_df$Substrate,permutations=1000)
am_model

### Indicator Species Analysis
otu_mp <- multipatt(as.data.frame(t(rare_otu)), map_wo_negs["Substrate",], control=how(nperm=9999))
summary(otu_mp, indvalcomp=TRUE)
write.table(as.data.frame(otu_mp$sign), file = "./Stats/OTU_indicspecies.csv", na = "", sep=",", quote = F)
otu_mp -> otu_mp_fdr
otu_mp_fdr$sign$p.value <- p.adjust(otu_mp$sign$p.value, "fdr")
summary(otu_mp_fdr)

write.table(as.data.frame(otu_mp_fdr$sign), file = "./Stats/OTU_indicspecies_fdr.csv", na = "", sep=",", quote = F)
otu_mp_fdr$sign

### Alpha Diversity Analysis
diver_df <- data.frame(shannon = diversity(t(rare_otu), index = "shannon"),
                       inv_simp = diversity(t(rare_otu), index = "inv"),
                       Substrate = as.factor(map_wo_negs["Substrate",]),
                       Plant_species = as.factor(map_wo_negs["Plant_species",]),
                       Site = as.factor(map_wo_negs["Site",]))
tibble(diver_df)     

m1 <- lm(shannon ~ Substrate + Plant_species + Site, diver_df)
summary(m1)

m2 <- aov(shannon ~ Substrate + Plant_species + Site, diver_df)
summary(m2)

m3 <- aov(inv_simp ~ Substrate + Plant_species + Site, diver_df)
summary(m3)

