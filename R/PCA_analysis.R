### From the BC_data I prepare to do some PCA visualization----
### Now we use the BC_clean data to ruduce the dimensions 
### 

library(tidyverse)
library(broom)  
library(cowplot)
library(patchwork) # required to arrange plots side-by-side
library(ggthemes) # for colorblind color scale
library(dplyr)
library(vroom)
library(purrr)
rm(list = ls())

### load the data 
BC_overlap_genes <- read.csv('results/06_BC_overlap_genes.csv')
BC_overlap_genes

### Principal Conponent Analysis
BC_overlap_genes_pca <- BC_overlap_genes %>% select(-c(PAM50.mRNA)) %>% 
prcomp(center = TRUE, scale. = T)

### 
BC_overlap_genes_pca %>% tidy('pcs') %>% 
  ggplot(aes(x = PC,
             y = percent)) + geom_col() +
  theme_bw()

BC_pca_plt <- data.frame(BC_overlap_genes_pca$x,subtype = BC_overlap_genes$PAM50.mRNA)

### augment the data  
BC_aug <- BC_overlap_genes_pca %>% augment(BC_overlap_genes)
BC_aug

### the bar plot
BC_aug %>% ggplot(aes(x= .fittedPC1, 
                      y = .fittedPC2,
                      col = PAM50.mRNA)) + geom_point() +
scale_color_colorblind()

### cluster 
BC_K_org <- BC_aug %>% 
  select(contains('NP')) %>%
  kmeans(centers = 4)
BC_K_org

BC_aug_K_org <- BC_K_org %>% augment(BC_aug) %>% 
  rename(cluster_org = .cluster)
BC_aug_K_org

### compare the Kmean cluster to the grandtruth 
pl1 = ggplot(BC_aug_K_org,aes(x=.fittedPC1,
                              y=.fittedPC2,
                              col = PAM50.mRNA)) + geom_point() +
  theme(legend.position = "bottom")


pl2 = ggplot(BC_aug_K_org,aes(x=.fittedPC1,
                              y=.fittedPC2,
                              col =cluster_org)) + geom_point() +
  theme(legend.position = "bottom")

pl1 + pl2

### calculate the accuracy of the prediction
BC_aug_K_org %>% 
  select(PAM50.mRNA,cluster_org) %>% 
  mutate(predict = case_when(cluster_org == 1 ~ 'HER2-enriched',
                             cluster_org == 2 ~ 'Luminal B',
                             cluster_org == 3 ~ 'Basal-like',
                             cluster_org == 4 ~ 'Luminal A'),
         cluster_org_correct = case_when(predict == PAM50.mRNA ~ 1,
                                         predict != PAM50.mRNA ~ 0)) %>% 
  summarise(score_org = mean(cluster_org_correct))

### from the result we can see that we have only 0.418 accuracy. 
### seems not so good.


