# Load libraries ----------------------------------------------------------
library(tidyverse)
library(broom)  
library(cowplot)
library(patchwork) # required to arrange plots side-by-side
library(ggthemes) # for colorblind color scale
library(dplyr)
library(vroom)
library(purrr)
rm(list = ls())

# Load data ----------------------------------------------------------------
BC_overlap_genes <- read.csv('results/06_BC_overlap_genes.csv')


# PCA analysis -----------------------------------------------------------------

pca_BC_overlap <- BC_overlap_genes %>% 
  select(-PAM50.mRNA) %>% 
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_BC_overlap_aug <- pca_BC_overlap %>% 
  augment(BC_overlap_genes)


# POSSIBLE TO FIND THE 1:17 by setting some sort of condition instead??
set.seed(7)

# K-means analysis -------------------------------------------------------------

# Extracting the 17 PCA that explains 95% of the variance
k_pca_BC_overlap <- pca_BC_overlap_aug %>% 
  select(str_c(".fittedPC", 1:17)) %>%
  kmeans(centers = 4)

# Adding back the augmented PCA data
pca_aug_k_pca_BC_overlap <- k_pca_BC_overlap %>% 
  augment(pca_BC_overlap_aug) %>% 
  rename(cluster_pca_CommonGenes = .cluster)


# Visualizing the cumulative variance  -----------------------------------------

pca_BC_overlap %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "BC_overlap_genes") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 3,
                y = 0.90,
                label = "95%",
                vjust = -1)) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  theme_half_open(12) +
  theme(text = element_text(family = "Avenir",
                            size = 12))

# K-means - Scatter plot -------------------------------------------------------

# Plot coloured according to subtypes
plot_k_pca_BC_overlap_subtypes <- pca_aug_k_pca_BC_overlap %>% 
  ggplot(mapping = aes(x= .fittedPC1, 
                       y = .fittedPC2,
                       colour = PAM50.mRNA)) + 
  geom_point() + 
  labs(x = "PC1 (XX%)",
       y = "PC2 (XXX%)",
       color = "Subtype") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278"))


# Plot coloured according to cluster
plot_k_pca_BC_overlap_cluster <- pca_aug_k_pca_BC_overlap %>% 
  ggplot(mapping = aes(x= .fittedPC1, 
                       y = .fittedPC2,
                       colour = cluster_pca_CommonGenes)) + 
  geom_point() + 
  labs(x = "PC1 (XX%)",
       y = "PC2 (XXX%)",
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#00688B",
                                "2" = "#00CD66",
                                "3" ="#FFA500",
                                "4" = "#CD3278"))


# Both plots together
(plot_k_pca_BC_overlap_subtypes/plot_pca_aug_k_pca_BC_overlap) &
  theme_half_open(12) &
  plot_annotation(title = "Common genes",
                  theme = theme(plot.title = element_text(size = 14))) &
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        text = element_text(family = "Avenir",
                            size = 12))



#4) Calculate the accuracy of the predictions ?????
pca_aug_K_pca %>% 
  select(PAM50.mRNA, cluster_org_CommonGenes, cluster_pca_CommonGenes) %>% 
  mutate(cluster_org_CommonGenes = case_when(cluster_org_CommonGenes == 1 ~ 'Basal-like',
                                             cluster_org_CommonGenes == 2 ~ 'Luminal A',
                                             cluster_org_CommonGenes == 3 ~ 'Luminal B',
                                             cluster_org_CommonGenes == 4 ~ 'HER2-enriched'),
         cluster_pca_CommonGenes = case_when(cluster_pca_CommonGenes == 1 ~ 'Luminal A',
                                             cluster_pca_CommonGenes == 2 ~ 'HER2-enriched',
                                             cluster_pca_CommonGenes == 3 ~ 'Basal-like',
                                             cluster_pca_CommonGenes == 4 ~ 'Luminal B'),
         cluster_org_CommonGenes_correct = case_when(cluster_org_CommonGenes == PAM50.mRNA ~ 1,
                                                     cluster_org_CommonGenes != PAM50.mRNA ~ 0),
         cluster_pca_CommonGenes_correct = case_when(cluster_pca_CommonGenes == PAM50.mRNA ~ 1,
                                                     cluster_pca_CommonGenes != PAM50.mRNA ~ 0)) %>% 
  summarise(score_org_CommonGenes = mean(cluster_org_CommonGenes_correct),
            score_pca_CommonGenes = mean(cluster_pca_CommonGenes_correct))


# Calculate BSS/TSS ratio
k_pca$betweenss/k_pca$totss

# Trash?? ----------------------------------------------------------------------


#pca_aug <- pca_fit %>% 
#  augment(BC_overlap_genes)

#Clustering by kmeans 
## Cluster PCA data



#Individual variance
#pca_fit %>%
#  tidy("pcs") %>%
#  filter(PC <= 10) %>%  # Only look at top10th PC
#  ggplot(aes(PC, percent)) +
#  geom_col(fill = "#56B4E9", alpha = 0.8) +
#  scale_x_continuous(breaks = 1:24)+
#  scale_y_continuous(
#    labels = scales::percent_format(),
#    expand = expansion(mult = c(0, 0.01))
#  ) +
#  theme_minimal_hgrid(12)

#2) PC1 & PC2 scatter plot on original data


#3) Compare the clustering results to the ground truth 
#pl1 <- pca_aug_K_pca %>%
#  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = PAM50.mRNA)) +
#  geom_point() +
#  theme(legend.position = "bottom")

#pl2 <- pca_aug_K_pca %>%
#  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org_CommonGenes)) +
#  geom_point() +
#  theme(legend.position = "bottom")

#pl3 <- pca_aug_K_pca %>%
#  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_pca_CommonGenes)) +
#  geom_point() +
#  theme(legend.position = "bottom")

#pl1 + pl3
