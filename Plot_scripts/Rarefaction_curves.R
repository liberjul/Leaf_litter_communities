library(ggplot2)
library(ggpubr)
library(vegan)
library(dplyr)
library(purrr)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
sl <- c(1:9, 11, 13:21, 25:36, 38, 40:44, 48:55, 57:58) # Slice for non-negative or failed samples
otu_dat <- read.table("./Data/all_OTUS_R1_clean.txt", sep="\t", header=TRUE) # Read OTU table with contaminants removed
rownames(otu_dat) <- otu_dat[,1] # Rename rows with OTU number
otu_dat <- otu_dat[,order(colnames(otu_dat))] # Reorder by OTU number
map <- t(read.table("./Data/DEM_map_mod3.csv", fill = TRUE, sep=",", header=TRUE)) # Read the map file
colnames(map) <- map[1,] # Colnames are now by Sample ID
rownames(map) <- c("SampleID", rownames(map)[2:23]) # Rownames are now site/sample variables
map_wo_negs <- map[,sl] # Keep non-negative samples for map table
colnames(otu_dat) <- c("OTU_ID", colnames(map)) # OTU_ID as first column, sample IDs as rest of columns
otu_dat_wo_negs <- as.matrix(otu_dat[,sl + 1]) # Keep non-negative samples for OTU table

out_rc <- rarecurve(t(otu_dat_wo_negs)) # Rarefaction curve vegan
names(out_rc) <- map_wo_negs["Short_name",] # rename output objects
### Some code from http://r-sig-ecology.471788.n2.nabble.com/Rarefaction-curves-in-ggplot-td7580273.html
as_tibble_rc <- function(x){
  # convert rarecurve() output to data.frame
  # bind_rows doesn't work because items are different lengths
  # also need to extract sample sizes from attribute
  # Allocate result dataframe
  nsamples <- map_int(x, length)
  total_samples <- sum(nsamples)
  if(!is.null(names(x))){
    sites <- names(x)
  } else {
    sites <- as.character(1:length(nsamples))
  }
  result <- data_frame(Site = rep("", total_samples),
                       Sample_size = rep(0, total_samples),
                       Species = rep(0, total_samples))
  start <- 1
  for (i in 1:length(nsamples)){
    result[start:(start + nsamples[i]-1), "Site"] <- sites[i]
    result[start:(start + nsamples[i]-1), "Sample_size"] <- attr(x[[i]], "Subsample")
    result[start:(start + nsamples[i]-1), "Species"] <- x[[i]]
    start <- start + nsamples[i]
  }
  result
}

out <- as_tibble_rc(out_rc)
# add grouping variable
sitedata <- data_frame(Site = map_wo_negs["Short_name",],
                       Substrate = map_wo_negs["Substrate",])
alldata <- left_join(out, sitedata, by = "Site")


rc_gplot <- ggplot(data = alldata,
       aes(x = Sample_size, y = Species, color = Substrate, group = Site)) +
  geom_line() +
  theme_pubr() +
  theme(legend.position = "right") +
  labs(x = "Read Count", y = "Observed OTUs") +
  scale_color_discrete(breaks = c("Endo", "Epi", "Lit", "Soil"),
                     labels = c("Endophyte", "Epiphyte", "Litter", "Soil"))

rc_gplot
ggsave("Rarefaction_curve.png", rc_gplot, width=6, height=6, units="in")
