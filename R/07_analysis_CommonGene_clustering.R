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


# Model data ---------------------------------------------------------------
#1) PCA
# Run pca
pca_fit <- BC_overlap_genes %>% 
  select(-PAM50.mRNA) %>% 
  prcomp(center = TRUE, scale. = T)

# Augment data
pca_aug <- pca_fit %>% 
  augment(BC_overlap_genes)


#Clustering by kmeans 
## Cluster original data
pca_aug_K_org <- pca_aug %>% 
  select(contains('NP')) %>%
  kmeans(centers = 4) %>% 
  augment(pca_aug) %>% 
  rename(cluster_org_CommonGenes = .cluster)

## Cluster PCA data
pca_aug_K_pca <- pca_aug_K_org %>% 
  select(str_c(".fittedPC", 1:17)) %>%
  kmeans(centers = 4) %>% 
  augment(pca_aug_K_org) %>% 
  rename(cluster_pca_CommonGenes = .cluster)




# Visualize data ------------------------------------------------------------
#1) Variance percent
## Cumulative variance
pca_fit %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:24)+
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))) +
  theme_minimal_hgrid(12) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_text(aes(0, 0.95,
                label = "95%", 
                vjust = -1))

## Individual variance
pca_fit %>%
  tidy("pcs") %>%
  filter(PC <= 10) %>%  # Only look at top10th PC
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:24)+
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)

#2) PC1 & PC2 scatter plot on original data
ggplot(data = pca_aug,
       mapping = aes(x= .fittedPC1, 
                      y = .fittedPC2,
                      col = PAM50.mRNA)) + 
  geom_point() +
  scale_color_colorblind()

#3) Compare the clustering results to the ground truth 
pl1 <- pca_aug_K_pca %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = PAM50.mRNA)) +
  geom_point() +
  theme(legend.position = "bottom")

pl2 <- pca_aug_K_pca %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org_CommonGenes)) +
  geom_point() +
  theme(legend.position = "bottom")

pl3 <- pca_aug_K_pca %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_pca_CommonGenes)) +
  geom_point() +
  theme(legend.position = "bottom")

pl1 + pl2 + pl3

#4) Calculate the accuracy of the predictions
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

### from the result we can see that we have only 0.418 accuracy. 
### seems not so good.


