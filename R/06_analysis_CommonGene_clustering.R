# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
library("ggthemes")
rm(list = ls())

# Define functions -------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data ----------------------------------------------------------------
BC_overlap_genes <- read.csv('results/06_BC_overlap_genes.csv')

# PCA analysis -----------------------------------------------------------------
pca_BC_overlap <- BC_overlap_genes %>% 
  select(-PAM50.mRNA) %>% 
  prcomp(center = TRUE, 
         scale. = TRUE)

# Adding original data back
pca_BC_overlap_aug <- pca_BC_overlap %>% 
  augment(BC_overlap_genes)

# K-means analysis -------------------------------------------------------------

# Extract the maximum PC that explain 95% of the cumulative variance
max_PC <- pca_BC_overlap %>%
  tidy("pcs") %>% 
  filter(cumulative <= 0.96) %>% 
  summarise(max = max(PC)) %>% 
  pull()

# K-means clustering
set.seed(7)
k_pca_BC_overlap <- pca_BC_overlap_aug %>% 
  select(str_c(".fittedPC", 
               1:max_PC)) %>%
  kmeans(centers = 4)

# Adding back the augmented PCA data
pca_aug_k_pca_BC_overlap <- k_pca_BC_overlap %>% 
  augment(pca_BC_overlap_aug) %>% 
  rename(cluster_pca_CommonGenes = .cluster)

# Plots of the cumulative variance  --------------------------------------------
plot_bar_PC_CumVar <- pca_BC_overlap %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(subtitle = "a) Variance explained the PCs",
       y = "Cumulative variance") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 2,
                y = 0.95,
                label = "95%",
                vjust = -1)) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  new_theme

# Scatter plot of K-means clustering -------------------------------------------

# Extract the variance explained by PC1
pc1 <- pca_BC_overlap %>%
  tidy("pcs") %>% 
  filter(PC == 1) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent)

# Extract the variance explained by PC2
pc2 <- pca_BC_overlap %>%
  tidy("pcs") %>% 
  filter(PC == 2) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent) 

# Plot coloured according to subtypes
plot_k_pca_BC_overlap_subtypes <- pca_aug_k_pca_BC_overlap %>% 
  ggplot(mapping = aes(x= .fittedPC1, 
                       y = .fittedPC2,
                       colour = PAM50.mRNA)) + 
  geom_point() + 
  labs(subtitle = "b) K-means clustering after PCA",
       x = str_c("PC1 (", pc1, "%)"),
       y = str_c("PC2 (", pc2, "%)"),
       color = "Subtype") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278")) + 
  new_theme +
  theme(legend.position = "right")

# Plot coloured according to cluster
plot_k_pca_BC_overlap_cluster <- pca_aug_k_pca_BC_overlap %>% 
  ggplot(mapping = aes(x= .fittedPC1, 
                       y = .fittedPC2,
                       colour = cluster_pca_CommonGenes)) + 
  geom_point() + 
  labs(x = str_c("PC1 (", pc1, "%)"),
       y = str_c("PC2 (", pc2, "%)"),
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#00688B",
                                "2" = "#CD3278",
                                "3" ="#FFA500",
                                "4" = "#00CD66")) +
  new_theme +
  theme(legend.position = "right")

# Final plot (cumulative variance and K-means clustering) ----------------------
plot_bar_PC_CumVar + 
  (plot_k_pca_BC_overlap_subtypes / plot_k_pca_BC_overlap_cluster) + 
  plot_annotation(title = "BC overlap data")


# Save plot --------------------------------------------------------------------
ggsave(file = "results/06_BC_overlap_PCA.png",
       width = 10, 
       height = 6.25, 
       dpi = 150)

# The accuracy of the predictions for this exact K-means analysis --------------
select(PAM50.mRNA, cluster_pca_CommonGenes) %>% 
  mutate(cluster_pca_CommonGenes = case_when(
    cluster_pca_CommonGenes == 1 ~ 'Luminal A',
    cluster_pca_CommonGenes == 2 ~ 'HER2-enriched',
    cluster_pca_CommonGenes == 3 ~ 'Basal-like',
    cluster_pca_CommonGenes == 4 ~ 'Luminal B'),
    cluster_pca_CommonGenes_correct = case_when(
      cluster_pca_CommonGenes == PAM50.mRNA ~ 1,
      cluster_pca_CommonGenes != PAM50.mRNA ~ 0)) %>% 
  summarise(score_pca_CommonGenes = mean(cluster_pca_CommonGenes_correct))
