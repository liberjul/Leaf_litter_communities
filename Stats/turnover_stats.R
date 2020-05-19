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

combs <- c()
subs <- unique(turndf_ITS_exp1_long$Substrate)
turn_df_paired <- data.frame()
for (i in subs){
  for (j in subs){
    if (i != j & ! paste(j, i, sep="-") %in% combs){
      turndf_ITS_exp1_long %>%
        mutate(Substrate = as.character(Substrate)) %>%
        filter(Substrate %in% c(i, j)) %>%
        mutate(Substrate = as.numeric(recode_factor(Substrate,
                                                    i = "1", j = "2"))) %>%
        CalcTurn %>%
        mutate(Substrate = recode_factor(Substrate,
                                         "1" = i, "2" = j),
               Paired_sub = paste(i, Substrate, sep="-")) %>%
        bind_rows(turn_df_paired) -> turn_df_paired
      
      word <- paste(i, j, sep="-")
      combs <- c(combs, word)
      print(word) 
    }
  }
}
turn_df_paired

paired_plot <- turn_df_paired %>%
  mutate(metric = recode_factor(metric,
                "drop" = "Drop",
                "gain" = "Gain",
                "total" = "Total")) %>%
  ggplot(aes(x = Paired_sub, y = turnover, color = metric)) +
  geom_point(size = 0.5) +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  scale_color_manual(values = c("#332288","#88CCEE","black")) +
  ylim(0, 1) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1, size = 10)) +
  labs(x = "Paired substrates", y = "Turnover", color = "Metric")
paired_plot

ggsave("./Figures/Turnover_paired.png", paired_plot, width = 8, height = 6, units="in")
ggsave("./Figures_Color/Figure 4.pdf", paired_plot, width = 90, height = 80, units = "mm")

paired_plot_gray <- turn_df_paired %>%
  mutate(metric = recode_factor(metric,
                                "drop" = "Drop",
                                "gain" = "Gain",
                                "total" = "Total")) %>%
  ggplot(aes(x = Paired_sub, y = turnover, color = metric)) +
  geom_point(size = 0.5) +
  geom_boxplot(alpha = 0) +
  theme_pubr() +
  scale_color_manual(values = c("#a6a6a6","#595959","black")) +
  ylim(0, 1) +
  theme(legend.position = "right",
        text = element_text(size = 10),
        axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1, size = 10)) +
  labs(x = "Paired substrates", y = "Turnover", color = "Metric")
paired_plot_gray
ggsave("./Figures_Numbered/Figure 4.pdf", paired_plot_gray, width = 90, height = 80, units = "mm")
