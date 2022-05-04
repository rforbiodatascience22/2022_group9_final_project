# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("ggvenn")
library("ggplot2")
library("patchwork")
rm(list=ls())

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load files -------------------------------------------------------------------
Basal_glm <- read_csv(file = "results/05_Basal_glm.csv")
Her2_glm <- read_csv(file = "results/05_Her2_glm.csv")
LumA_glm <- read_csv(file = "results/05_LumA_glm.csv")
LumB_glm <- read_csv(file = "results/05_LumB_glm.csv")
proteomes_clean <- read_csv(file = "data/02_proteomes_clean.csv")
BC_data_clean_aug <- read_csv(file = "data/03_BC_data_clean_aug.csv")

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
    pull(proteome)
)
# Plot Venn diagram
ggvenn(significant_proteins)

# Find overlapping proteins
overlap_pro <- intersect(intersect(intersect(significant_proteins[["Basal"]], 
                                             significant_proteins[["Her2"]]), 
                                   significant_proteins[["LumA"]]),
                         significant_proteins[["LumB"]])

# Find unique proteins
Basal_unique <- setdiff(setdiff(setdiff(significant_proteins[["Basal"]],
                                significant_proteins[["Her2"]]),
                        significant_proteins[["LumA"]]),
                        significant_proteins[["LumB"]]) #959

Her2_unique <- setdiff(setdiff(setdiff(significant_proteins[["Her2"]],
                                       significant_proteins[["Basal"]]),
                                significant_proteins[["LumA"]]),
                        significant_proteins[["LumB"]]) #639

LumA_unique <- setdiff(setdiff(setdiff(significant_proteins[["LumA"]],
                                       significant_proteins[["Basal"]]),
                               significant_proteins[["Her2"]]),
                       significant_proteins[["LumB"]]) #618

LumB_unique <- setdiff(setdiff(setdiff(significant_proteins[["LumB"]],
                                       significant_proteins[["Basal"]]),
                               significant_proteins[["LumA"]]),
                       significant_proteins[["Her2"]]) #460


# Visualization ----------------------------------------------------------------
# Overlap gene expression heatmap
ggplot(data = proteomes_clean %>% 
         select(`TCGA ID`,
                c(overlap_pro)) %>% 
         pivot_longer(cols = -`TCGA ID`,
                      names_to = "proteome",
                      values_to = "expr_level"),
       mapping = aes(x = `TCGA ID`,
                     y = proteome,
                     fill = expr_level)) +
  geom_tile(alpha = 0.5) +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

pl1 <- ggplot(data = BC_data_clean_aug %>% 
         filter(Basal_like == 1) %>% 
         select(`Complete TCGA ID`,
                c(overlap_pro)) %>% 
         pivot_longer(cols = -`Complete TCGA ID`,
                      names_to = "proteome",
                      values_to = "expr_level"),
       mapping = aes(x = `Complete TCGA ID`,
                     y = proteome,
                     fill = expr_level)) +
  geom_tile(alpha = 0.5) +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Basal-like")

pl2 <- ggplot(data = BC_data_clean_aug %>% 
                filter(HER2_enriched == 1) %>% 
                select(`Complete TCGA ID`,
                       c(overlap_pro)) %>% 
                pivot_longer(cols = -`Complete TCGA ID`,
                             names_to = "proteome",
                             values_to = "expr_level"),
              mapping = aes(x = `Complete TCGA ID`,
                            y = proteome,
                            fill = expr_level)) +
  geom_tile(alpha = 0.5) +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Her2-enriched")

pl3 <- ggplot(data = BC_data_clean_aug %>% 
                filter(Luminal_A == 1) %>% 
                select(`Complete TCGA ID`,
                       c(overlap_pro)) %>% 
                pivot_longer(cols = -`Complete TCGA ID`,
                             names_to = "proteome",
                             values_to = "expr_level"),
              mapping = aes(x = `Complete TCGA ID`,
                            y = proteome,
                            fill = expr_level)) +
  geom_tile(alpha = 0.5) +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Luminal A")

pl4 <- ggplot(data = BC_data_clean_aug %>% 
                filter(Luminal_B == 1) %>% 
                select(`Complete TCGA ID`,
                       c(overlap_pro)) %>% 
                pivot_longer(cols = -`Complete TCGA ID`,
                             names_to = "proteome",
                             values_to = "expr_level"),
              mapping = aes(x = `Complete TCGA ID`,
                            y = proteome,
                            fill = expr_level)) +
  geom_tile(alpha = 0.5) +
  scale_fill_gradient2(midpoint = 0,
                       low = "blue",
                       mid = "white",
                       high = "red") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(title = "Luminal B")

(pl1+pl2)/(pl3+pl4)
