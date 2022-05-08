# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("ggvenn")
library("ggplot2")
library("patchwork")
rm(list=ls())

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data --------------------------------------------------------------------
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")


# Model data -------------------------------------------------------------------
# Making four different logistic models for each of the subtypes
# Luminal A
LumA_glm <- subtype_glm("Luminal_A", BC_data_clean_aug)

# Luminal B
LumB_glm <- subtype_glm("Luminal_B", BC_data_clean_aug)

# HER2_enriched
Her2_glm <- subtype_glm("HER2_enriched", BC_data_clean_aug)

# Basal-like
Basal_glm <- subtype_glm("Basal_like", BC_data_clean_aug)


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

# Plot Venn diagram ------------------------------------------------------------
significant_proteins %>% 
  ggvenn()

# Save plot --------------------------------------------------------------------
ggsave(file = "results/06_venndiagram.png",
       width = 10, 
       height = 7, 
       dpi = 150)


# Find overlapping proteins ----------------------------------------------------
overlap_pro <- intersect(intersect(intersect(significant_proteins[["Basal"]], 
                                             significant_proteins[["Her2"]]), 
                                   significant_proteins[["LumA"]]),
                         significant_proteins[["LumB"]])

# Is there a way to do this not using base R? ^^^

# Find unique proteins (MIGHT DELETE)
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


# Build new data frame for overlap proteins
BC_overlap_genes <- BC_data_clean_aug %>% 
  select(c(overlap_pro),
         `PAM50 mRNA`)

# Visualization ----------------------------------------------------------------
# Overlap gene expression heatmap
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

(pl1+pl2)/(pl3+pl4) +
  plot_annotation(title = "Protein Expression of 24 common signigicant genes for each patient grouped by tumor type")

# Save plot --------------------------------------------------------------------
ggsave(file = "results/06_subtype_heatmap.png",
       width = 10, 
       height = 7, 
       dpi = 150)


# Save files -------------------------------------------------------------------
save(significant_proteins,
     file = "results/06_significant_proteins.rda")
write_csv(x = BC_overlap_genes,
          file = "results/06_BC_overlap_genes.csv")
