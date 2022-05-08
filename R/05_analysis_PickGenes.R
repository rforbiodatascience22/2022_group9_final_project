# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("ggvenn")
library("ggplot2")
library("patchwork")
rm(list=ls())

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data --------------------------------------------------------------------
BC_data_clean_aug <- read_csv(file = "data/03_BC_data_clean_aug.csv")


# Model data -------------------------------------------------------------------

# Make/Load 4 different logistic models for each subtype
# NOTE: If there're certain glm result files, the program will directly read
# the existing files instead of running glm again. Be aware of this!!!

# Luminal A
if (!file.exists("results/05_LumA_glm.csv")){
  LumA_glm <- subtype_glm("Luminal_A", 
                          BC_data_clean_aug)
  write_csv(LumA_glm,
            file = "results/05_LumA_glm.csv")
}else{
  LumA_glm <- read_csv("results/05_LumA_glm.csv")
}


# Luminal B
if (!file.exists("results/05_LumB_glm.csv")){
  LumB_glm <- subtype_glm("Luminal_B", 
                          BC_data_clean_aug)
  write_csv(LumB_glm,
            file = "results/05_LumB_glm.csv")
}else{
  LumB_glm <- read_csv("results/05_LumB_glm.csv")
}

# HER2 enriched
if (!file.exists("results/05_Her2_glm.csv")){
  Her2_glm <- subtype_glm("HER2_enriched", 
                          BC_data_clean_aug)
  write_csv(Her2_glm,
            file = "results/05_Her2_glm.csv")
}else{
  Her2_glm <- read_csv("results/05_Her2_glm.csv")
}

# Basal-like
if (!file.exists("results/05_Basal_glm.csv")){
  Basal_glm <- subtype_glm("Basal_like", 
                           BC_data_clean_aug)
  write_csv(Basal_glm,
            file = "results/05_Basal_glm.csv") 
}else{
  Basal_glm <- read_csv("results/05_Basal_glm.csv")
}


# Wrangle data -----------------------------------------------------------------

# Find significant proteins
significant_proteins <- list(
  Basal = Basal_glm %>% 
    filter(identified_as == "significant") %>% 
    pull(proteome),
  Her2 = Her2_glm %>% 
    filter(identified_as == "significant") %>% 
    pull(proteome),
  LumA = LumA_glm %>% 
    filter(identified_as == "significant") %>% 
    pull(proteome),
  LumB = LumB_glm %>% 
    filter(identified_as == "significant") %>% 
    pull(proteome))


# Find overlapping proteins
overlap_pro <- intersect(intersect(intersect(significant_proteins[["Basal"]], 
                                             significant_proteins[["Her2"]]), 
                                   significant_proteins[["LumA"]]),
                         significant_proteins[["LumB"]])


# Build new data frame for overlap proteins
BC_overlap_genes <- BC_data_clean_aug %>% 
  select(c(overlap_pro),
         `PAM50 mRNA`)

# Visualize data ------------------------------------------------------------

# Venn diagram 
significant_proteins %>% 
  ggvenn()

# Save plot -----------------------------------
ggsave(file = "results/05_venndiagram.png",
       width = 10, 
       height = 7, 
       dpi = 150)


# Overlap gene expression heatmap for Basal like subtype
Basal_pl1 <- ggplot(data = BC_data_clean_aug %>% 
                      filter(Basal_like == 1) %>% 
                      select(`Complete TCGA ID`,
                             c(overlap_pro)) %>% 
                      pivot_longer(cols = -`Complete TCGA ID`,
                                   names_to = "proteome",
                                   values_to = "expr_level"),
                    mapping = aes(x = `Complete TCGA ID`,
                                  y = proteome,
                                  fill = expr_level)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red",
                       limits = c(-7.5, 6)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Basal-like",
       fill = "Expression level")

# Overlap gene expression heatmap for HER2 enriched subtype
Her2_pl2 <- ggplot(data = BC_data_clean_aug %>% 
                     filter(HER2_enriched == 1) %>% 
                     select(`Complete TCGA ID`,
                            c(overlap_pro)) %>% 
                     pivot_longer(cols = -`Complete TCGA ID`,
                                  names_to = "proteome",
                                  values_to = "expr_level"),
                   mapping = aes(x = `Complete TCGA ID`,
                                 y = proteome,
                                 fill = expr_level)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red",
                       limits = c(-7.5, 6)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Her2-enriched",
       fill = "Expression level")

# Overlap gene expression heatmap for Luminal A subtype
LumA_pl3 <- ggplot(data = BC_data_clean_aug %>% 
                     filter(Luminal_A == 1) %>% 
                     select(`Complete TCGA ID`,
                            c(overlap_pro)) %>% 
                     pivot_longer(cols = -`Complete TCGA ID`,
                                  names_to = "proteome",
                                  values_to = "expr_level"),
                   mapping = aes(x = `Complete TCGA ID`,
                                 y = proteome,
                                 fill = expr_level)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red",
                       limits = c(-7.5, 6)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Luminal A",
       fill = "Expression level")

# Overlap gene expression heatmap for Luminal B subtype
LumB_pl4 <- ggplot(data = BC_data_clean_aug %>% 
                     filter(Luminal_B == 1) %>% 
                     select(`Complete TCGA ID`,
                            c(overlap_pro)) %>% 
                     pivot_longer(cols = -`Complete TCGA ID`,
                                  names_to = "proteome",
                                  values_to = "expr_level"),
                   mapping = aes(x = `Complete TCGA ID`,
                                 y = proteome,
                                 fill = expr_level)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red",
                       limits = c(-7.5, 6)) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Luminal B",
       fill = "Expression level")


# Overlap gene expression heatmap for all subtypes
(Basal_pl1+Her2_pl2)/(LumA_pl3+LumB_pl4) + 
  plot_layout(guides = "collect") & 
  plot_annotation(title = "Protein expression of the 24 common significant genes") 


# Save plot --------------------------------------------------------------------
ggsave(file = "results/05_subtype_heatmap.png",
       width = 10, 
       height = 7, 
       dpi = 150)


# Save files -------------------------------------------------------------------
save(significant_proteins,
     file = "results/06_significant_proteins.rda")
write_csv(x = BC_overlap_genes,
          file = "results/06_BC_overlap_genes.csv")