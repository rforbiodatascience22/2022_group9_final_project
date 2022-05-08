library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
library("cowplot")

# Define functions -------------------------------------------------------------
source(file = "./R/99_project_functions.R")

# Load data ----------------------------------------------------------------
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
BC_data_PAM50_clean <- read_csv(file = "data/02_BC_data_PAM50_clean.csv")

# PCA analysis  ----------------------------------------------------------------

# For full BC_data set
pca_org <- BC_data_clean_aug %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>%
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_org_aug <- pca_org %>% 
  augment(BC_data_clean_aug)

# Reduced version (protein IDs common in PAM50 and BC_data_clean_aug)  ### New reference in stead on "reduced"
pca_red <- BC_data_PAM50_clean %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>%
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_red_aug <- pca_red %>% 
  augment(BC_data_PAM50_clean)


# K-means analysis -------------------------------------------------------------

# Yi Huang will fix something here *****

set.seed(4)

# For full BC_data set
k_org <- pca_org_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  select(c(1:64)) %>% 
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_org <- k_org %>%
  augment(pca_org_aug) %>% 
  rename(cluster_org = .cluster)


set.seed(4)

# For reduced version (protein IDs common in PAM50 and BC_data_clean_aug)
k_red <- pca_red_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  select(c(1:15)) %>% 
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_red <- k_red %>%
  augment(pca_red_aug) %>% 
  rename(cluster_red = .cluster)

# Visualizing the cumulative variance  -----------------------------------------

# For full BC_data set (PC 1 to 64 selected for K-means)
plot_pca_org_cum <- pca_org %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "a) BC_data") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 5,
                y = 0.90,
                label = "95%",
                vjust = -1))

# Reduced version (PC 1 to 15 chosen for K-means)
plot_pca_red_cum <- pca_red %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "b) BC_data_PAM50") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 2,
                y = 0.90,
                label = "95%",
                vjust = -1))

# Plot both
plot_pca_org_cum + 
  plot_pca_red_cum +
  plot_layout(guides = 'collect') &
  scale_x_continuous(expand = c(0, 0)) &
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) &
  theme_half_open(12) &
  theme(text = element_text(family = "Avenir",
                            size = 12))

# change title and 95% placement ^^^^

# Save plot ---------------------------------------------------------------
ggsave(file = "results/07_CumVar_Comparison.png",
       width = 8.96, 
       height = 5.42, 
       dpi = 150)



# K-means - Scatter plot -------------------------------------------------------

# Original data (colored according to PAM50 mRNA)
plot_pca_aug_k_org_subtypes <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  labs(x = "PC1 (10.7%)",
       y = "PC2 (7.82%)",
       color = "Subtype") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278"))

# Colored according to clusters (NB: WAY OF EXTRACTING PERCENTAGE OF PCA DIRECTLY?)
plot_pca_aug_k_org_clusters <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_org)) +
  geom_point() +
  labs(x = "PC1 (10.7%)",
       y = "PC2 (7.82%)",
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#CD3278",
                                "2" = "#FFA500",
                                "3"="#00688B",
                                "4" = "#00CD66"))

# Plot of both
(plot_pca_aug_k_org_subtypes/plot_pca_aug_k_org_clusters) &
  theme_half_open(12) &
  plot_annotation(title = "Full breast cancer data",
                  theme = theme(plot.title = element_text(size = 14))) &
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        text = element_text(family = "Avenir",
                            size = 12))


# Save plot ---------------------------------------------------------------
ggsave(file = "results/07_BC_data_K_means.png",
       width = 8.56, 
       height = 6.42, 
       dpi = 150)


# Scatter plot reduced version (subtype)
plot_pca_aug_k_red_subtypes <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  labs(x = "PC1 (34.9%)",
       y = "PC2 (12.7%)",
       color = "Subtype") +
  scale_color_manual(values = c("Basal-like" = "#00688B",
                                "HER2-enriched" = "#00CD66",
                                "Luminal A"="#FFA500",
                                "Luminal B" = "#CD3278")) +
  theme_half_open(12) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        text = element_text(family = "Avenir",
                            size = 12))
  

# Scatter plot reduced version (cluster)
plot_pca_aug_k_red_cluster <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_red)) +
  geom_point() +
  labs(x = "PC1 (34.9%)",
       y = "PC2 (12.7%)",
       color = "Cluster") +
  scale_color_manual(values = c("1" = "#FFA500",
                                "2" = "#00CD66",
                                "3" = "#00688B",
                                "4" = "#CD3278")) +
  theme_half_open(12) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        text = element_text(family = "Avenir",
                            size = 12))

# Comparison of subtype and cluster (reduced version)
(plot_pca_aug_k_red_subtypes + plot_pca_aug_k_red_cluster) &
  plot_annotation(title = "K-means clustering after PCA")



# Save plot ---------------------------------------------------------------
ggsave(file = "results/07_BC_data_PAM50_reduced.png",
       width = 8.56, 
       height = 6.42, 
       dpi = 150)


# Comparison of match ----------------------------------------------------------

# Original data
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

# Reduced data
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

# Calculate BSS/TSS ratio ------------------------------------------------------
