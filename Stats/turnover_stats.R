library(tidyr)
library(dplyr)
library(codyn)

# **************************************************------------
setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Undergrad/Fall 2018/PLP 847/miseq_dat/Leaf_litter_communities")
map_wo_negs <- as.matrix(read.csv("./Data/DEM_map_wo_negs.csv", stringsAsFactors = F))
rare_otu <- as.matrix(read.csv("./Data/Rare_otu_table.csv"))
taxa_table <- read.delim("./Data/consensus_taxonomy_constax.txt") # Load CONSTAX classifications

tax <- vector(length = dim(taxa_table)[1])
for (i in 1:length(tax)){
  a <- taxa_table[i,][taxa_table[i,] != ""]
  tax[i] <- a[length(a)]
}

taxa_table_new <- data.frame(OTU = taxa_table$OTU_ID,
                             Taxonomy=tax)
meta_dat <- as.data.frame(t(map_wo_negs))
meta_dat$SampleID <- rownames(meta_dat)

# CALCULATE TURNOVER ------------------------------------

# reformat the data --------------------------------------
GenLongTurn <- function(otu_ab, taxa, meta){
  data.frame(OTU = as.factor(row.names(otu_ab)), otu_ab) %>%
    tidyr::gather(SampleID, Abundance, -OTU) %>%
    dplyr::left_join(meta[, c("SampleID", "Site","Plant_species","Substrate")], by = "SampleID") %>%
    dplyr::left_join(taxa, by="OTU") %>%
    group_by(Site, Substrate, Taxonomy) %>%
    # summarise_(Abundance=sum(Abundance)) %>%
    summarise_if(is.numeric, ~ sum(.)) %>%
    as.data.frame() -> long_form
  return(long_form)
}

GenLongTurn(rare_otu, taxa_table_new, meta_dat) -> turndf_ITS_exp1_long
head(turndf_ITS_exp1_long)


# Calculating turnover ----------------------------------
CalcTurn <- function(df){
  turnover(df = df,
           time.var = "Substrate",
           species.var = "Taxonomy",
           abundance.var = "Abundance",
           replicate.var = "Site") -> df_tot
  # calcualting appear
  turnover(df = df,
           time.var = "Substrate",
           species.var = "Taxonomy",
           abundance.var = "Abundance",
           replicate.var = "Site",
           metric = "appearance") -> df_app
  # calcualting drops
  turnover(df = df,
           time.var = "Substrate",
           species.var = "Taxonomy",
           abundance.var = "Abundance",
           replicate.var = "Site",
           metric = "disappearance") -> df_dis
  df_tot$metric<-"total"
  names(df_tot)[1]="turnover"
  df_app$metric<-"gain"
  names(df_app)[1]="turnover"
  df_dis$metric<-"drop"
  names(df_dis)[1]="turnover"
  df_all <- rbind(df_tot, df_app, df_dis)
  return(df_all)
}
## H1 endo -> epi -> litter -> soil
turndf_ITS_exp1_long %>%
  mutate(Substrate = as.numeric(recode_factor(Substrate,
                                   "Endo" = "1",
                                   "Epi" = "2",
                                   "Lit" = "3",
                                   "Soil" = "4"))) %>%
  CalcTurn %>%
  mutate(Substrate = recode_factor(Substrate,
                                   "1" = "Endo",
                                   "2" = "Epi",
                                   "3" = "Lit",
                                   "4" = "Soil")) %>%
  ggplot(aes(x = Substrate,
             y = turnover,
             color = metric)) +
  labs(x = "Substrate", y = "Turnover", color = "Metric") +
  geom_point() +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  theme(legend.position = "right") -> h1_plot
h1_plot
## H2 epi -> endo -> soil -> litter
turndf_ITS_exp1_long %>%
  mutate(Substrate = as.numeric(recode_factor(Substrate,
                                              "Epi" = "1",
                                              "Endo" = "2",
                                              "Soil" = "3",
                                              "Lit" = "4"))) %>%
  CalcTurn %>%
  mutate(Substrate = recode_factor(Substrate,
                                   "1" = "Epi",
                                   "2" = "Endo",
                                   "3" = "Soil",
                                   "4" = "Lit")) %>%
  ggplot(aes(x = Substrate,
             y = turnover,
             color = metric)) +
  labs(x = "Substrate", y = NULL, color = "Metric") +
  geom_point() +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  theme(legend.position = "right") -> h2_plot
h2_plot
## H3 litter -> endophyte -> soil -> epi
turndf_ITS_exp1_long %>%
  mutate(Substrate = as.numeric(recode_factor(Substrate,
                                              "Lit" = "1",
                                              "Endo" = "2",
                                              "Soil" = "3",
                                              "Epi" = "4"))) %>%
  CalcTurn %>%
  mutate(Substrate = recode_factor(Substrate,
                                   "1" = "Lit",
                                   "2" = "Endo",
                                   "3" = "Soil",
                                   "4" = "Epi")) %>%
  ggplot(aes(x = Substrate,
             y = turnover,
             color = metric)) +
  labs(x = "Substrate", y = "Turnover", color = "Metric") +
  geom_point() +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  theme(legend.position = "right") -> h3_plot
h3_plot
## H4 epi -> endo -> lit -> soil
turndf_ITS_exp1_long %>%
  mutate(Substrate = as.numeric(recode_factor(Substrate,
                                              "Epi" = "1",
                                              "Endo" = "2",
                                              "Lit" = "3",
                                              "Soil" = "4"))) %>%
  CalcTurn %>%
  mutate(Substrate = recode_factor(Substrate,
                                   "1" = "Epi",
                                   "2" = "Endo",
                                   "3" = "Lit",
                                   "4" = "Soil")) %>%
  ggplot(aes(x = Substrate,
             y = turnover,
             color = metric)) +
  labs(x = "Substrate", y = NULL, color = "Metric") +
  geom_point() +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  # scale_x_discrete(limits = c("Epi", "Endo", "Lit", "Soil"),
  #                  breaks = c("Epi", "Endo", "Lit", "Soil")) +
  theme(legend.position = "right") -> h4_plot
h4_plot


g <- (h1_plot + h2_plot) /
  (h3_plot + h4_plot) + plot_layout(guides = "collect")
ggsave("./Figures/Turnover_4hypos.png", g, width = 12, height = 10, units = "in")