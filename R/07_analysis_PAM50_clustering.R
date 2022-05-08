library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
library("ggthemes")
rm(list = ls())

# Define functions -------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data --------------------------------------------------------------------
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
BC_data_PAM50_clean <- read_csv(file = "data/02_BC_data_PAM50_clean.csv")

# PCA analysis  ----------------------------------------------------------------

# For the BC data set
pca_org <- BC_data_clean_aug %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>%
  prcomp(center = TRUE, 
         scale. = TRUE)

# Adding original data back
pca_org_aug <- pca_org %>% 
  augment(BC_data_clean_aug)

# For the common protein IDs between PAM50 and BC data
pca_red <- BC_data_PAM50_clean %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>%
  prcomp(center = TRUE, 
         scale. = TRUE)

# Adding original data back
pca_red_aug <- pca_red %>% 
  augment(BC_data_PAM50_clean)

# K-means analysis -------------------------------------------------------------

# ***** For full BC_data set *****

# Extract the maximum PC that explain 95% of the cumulative variance
max_PC_org <- pca_org %>%
  tidy("pcs") %>% 
  filter(cumulative <= 0.96) %>% 
  summarise(max = max(PC)) %>% 
  pull()

# Clustering
set.seed(4)
k_org <- pca_org_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  select(c(1:max_PC_org)) %>% 
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_org <- k_org %>%
  augment(pca_org_aug) %>% 
  rename(cluster_org = .cluster)


# ***** For PAM50 version *****

# Extract the maximum PC that explain 95% of the cumulative variance
max_PC_red <- pca_red %>%
  tidy("pcs") %>% 
  filter(cumulative <= 0.96) %>% 
  summarise(max = max(PC)) %>% 
  pull()

# Clustering
set.seed(4)
k_red <- pca_red_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  select(c(1:max_PC_red)) %>% 
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_red <- k_red %>%
  augment(pca_red_aug) %>% 
  rename(cluster_red = .cluster)

# Plots of the cumulative variance  --------------------------------------------

# For the BC data set
plot_pca_org_cum <- pca_org %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "a) Variance explained the PCs") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 8,
                y = 0.94,
                label = "95%",
                vjust = -1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  new_theme

# For the common protein IDs between PAM50 and BC data
plot_pca_red_cum <- pca_red %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "a) Variance explained the PCs") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 3,
                y = 0.94,
                label = "95%",
                vjust = -1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) +
  new_theme

# K-means - Scatter plot -------------------------------------------------------

# ***** For full BC_data set *****

# Extract the variance explained by PC1 and PC2
org_pc1 <- pca_org %>%
  tidy("pcs") %>% 
  filter(PC == 1) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent)

org_pc2 <- pca_org %>%
  tidy("pcs") %>% 
  filter(PC == 2) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent) 

# Scatter plot colored according to subtype
plot_pca_aug_k_org_subtypes <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  labs(x = str_c("PC1 (", org_pc1, "%)"),
       y = str_c("PC2 (", org_pc2, "%)"),
       color = "Subtype",
       subtitle = "b) K-means clustering after PCA") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278")) +
  new_theme +
  theme(legend.position = "right")

# Scatter plot colored according to cluster
plot_pca_aug_k_org_clusters <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_org)) +
  geom_point() +
  labs(x = str_c("PC1 (", org_pc1, "%)"),
       y = str_c("PC2 (", org_pc2, "%)"),
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#CD3278",
                                "2" = "#FFA500",
                                "3"="#00688B",
                                "4" = "#00CD66")) +
  new_theme +
  theme(legend.position = "right")


# ***** For PAM50 version *****

# Extract the variance explained by PC1 and PC2
red_pc1 <- pca_red %>%
  tidy("pcs") %>% 
  filter(PC == 1) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent)

red_pc2 <- pca_red %>%
  tidy("pcs") %>% 
  filter(PC == 2) %>% 
  mutate(percent = round(percent * 100,
                         digits = 1)) %>% 
  pull(percent) 

# Scatter plot colored according to subtype
plot_pca_aug_k_red_subtypes <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  labs(x = str_c("PC1 (", red_pc1, "%)"),
       y = str_c("PC2 (", red_pc2, "%)"),
       color = "Subtype",
       subtitle = "b) K-means clustering after PCA") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278")) +
  new_theme +
  theme(legend.position = "right")
  

# Scatter plot colored according to cluster
plot_pca_aug_k_red_cluster <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_red)) +
  geom_point() +
  labs(x = str_c("PC1 (", red_pc1, "%)"),
       y = str_c("PC2 (", red_pc2, "%)"),
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#FFA500",
                                "2" = "#00CD66",
                                "3" = "#00688B",
                                "4" = "#CD3278")) +
  new_theme +
  theme(legend.position = "right")

# Cumulative variance and k-means clustering plots -----------------------------

# ***** For full BC_data set *****

plot_pca_org_cum + 
  (plot_pca_aug_k_org_subtypes/plot_pca_aug_k_org_clusters) + 
  plot_annotation(title = "BC data")

ggsave(file = "results/07_BC_data_cumulative_kmeans.png",
       width = 10, 
       height = 6.25, 
       dpi = 150)

# ***** For PAM50 version *****

plot_pca_red_cum + 
  (plot_pca_aug_k_red_subtypes/plot_pca_aug_k_red_cluster) + 
  plot_annotation(title = "Protein IDs common between PAM50 data and BC data")


ggsave(file = "results/07_BC_data_PAM50_cumulative_kmeans.png",
       width = 10, 
       height = 6.25, 
       dpi = 150)

# Comparison of match ----------------------------------------------------------

# ***** For full BC_data set *****

pca_aug_k_org %>%
  select(`PAM50 mRNA`, 
         cluster_org) %>%
  mutate(cluster_org = case_when(cluster_org == 1 ~ "Luminal B",
                                 cluster_org == 2 ~ "Luminal A",
                                 cluster_org == 3 ~ "Basal-like",
                                 cluster_org == 4 ~ "HER2-enriched"),
         cluster_pca_correct = case_when(`PAM50 mRNA` == cluster_org ~ 1,
                                         `PAM50 mRNA` != cluster_org ~ 0)) %>% 
  summarise(score_pca_org = mean(cluster_pca_correct))

# ***** For PAM50 version *****

pca_aug_k_red %>%
  select(`PAM50 mRNA`, 
         cluster_red) %>%
  mutate(cluster_red = case_when(cluster_red == 1 ~ "Luminal A",
                                 cluster_red == 2 ~ "HER2-enriched",
                                 cluster_red == 3 ~ "Basal-like",
                                 cluster_red == 4 ~ "Luminal B"),
         cluster_pca_correct = case_when(`PAM50 mRNA` == cluster_red ~ 1,
                                         `PAM50 mRNA` != cluster_red ~ 0)) %>% 
  summarise(score_pca_red = mean(cluster_pca_correct))
