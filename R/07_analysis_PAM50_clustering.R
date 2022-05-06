library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
library("cowplot")

# Define functions -------------------------------------------------------------
source(file = "./R/99_project_functions.R")

BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
BC_data_PAM50_clean <- read_csv(file = "data/02_BC_data_PAM50_clean.csv")

# PCA analysis  ----------------------------------------------------------------

# For full BC_data set
pca_ori <- BC_data_clean_aug %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_ori_aug <- pca_ori %>% 
  augment(BC_data_clean_aug)

# Reduced version (protein IDs common in PAM50 and BC_data_clean_aug)
pca_red <- BC_data_PAM50_clean %>% 
  select(starts_with(c("NP","XP","YP"))) %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_red_aug <- pca_red %>% 
  augment(BC_data_PAM50_clean)

# K-means analysis -------------------------------------------------------------
set.seed(4)

# For full BC_data set
k_org <- pca_ori_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_org <- k_org %>%
  augment(pca_ori_aug) %>% 
  rename(cluster_org = .cluster)

# For reduced version (protein IDs common in PAM50 and BC_data_clean_aug)
k_red <- pca_red_aug %>%
  select(starts_with(c("NP","XP","YP"))) %>%
  select(c(1:14)) %>% 
  kmeans(centers = 4)

# Adding original data back and renaming .cluster
pca_aug_k_red <- k_red %>%
  augment(pca_red_aug) %>% 
  rename(cluster_red = .cluster)

# Visualizing the results  -----------------------------------------------------

# For full BC_data set
plot_pca_ori_cumulative <- pca_ori %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "a) Full data set") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 5,
                y = 0.90,
                label = "95%",
                vjust = -1))


# Reduced version (protein IDs common in PAM50 and BC_data_clean_aug)
plot_pca_red_cumulative <- pca_red %>%
  tidy("pcs") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "turquoise3",
           alpha = 0.7) +
  labs(y = "Cumulative variance",
       subtitle = "b) Reduced version") +
  geom_hline(yintercept = 0.95,
             linetype = "dashed") +
  geom_text(aes(x = 2,
                y = 0.90,
                label = "95%",
                vjust = -1))

# Plot of the full and reduced version
(plot_pca_ori_cumulative + plot_pca_red_cumulative) +
  plot_layout(guides = 'collect') &
  scale_x_continuous(expand = c(0, 0)) &
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) &
  theme_half_open(12) &
  theme(text = element_text(family = "Avenir",
                            size = 12))

# Scatter plot (K-means) -------------------------------------------------------

# Original data

# Coloured according to PAM50 mRNA
pl1 <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

# Coloured according to clusters
pl2 <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")

# Plot of both
pl1 + pl2 +  plot_layout(guides = 'collect') &
  scale_x_continuous(expand = c(0, 0)) &
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) &
  theme_half_open(12) &
  theme(text = element_text(family = "Avenir",
                            size = 12))

# Reduced version
pl3 <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = `PAM50 mRNA`)) +
  geom_point() +
  labs(x = "PC1 (XX%)",
       y = "PC2 (YY%)",
       subtitle = "hhh") +
  scale_color_manual(values = c("Basal-like" = "#00CDCD",
                                "HER2-enriched" = "#9A32CD",
                                "Luminal A"="#66CD00",
                                "Luminal B" = "#CD6839")) +
  theme(legend.position = "bottom")


pl4 <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = cluster_red)) +
  geom_point() +
  labs(x = "PC1 (XX%)",
       y = "PC2 (YY%)",
       subtitle = "hhh") +
  theme(legend.position = "bottom") + 
  scale_color_manual(values = c("1" = "steelblue",
                                "2" = "orange",
                                "3"="purple",
                                "4" = "green"))


pl3 + pl4 + plot_layout(guides = 'collect') &
  theme_half_open(12) &
  theme(text = element_text(family = "Avenir",
                            size = 12))



# 
pca_aug_k_red %>%
  select(`PAM50 mRNA`, cluster_red, cluster_red) %>%
  mutate(cluster_red = case_when(cluster_red == 1 ~ "Luminal A",
                                 cluster_red == 2 ~ "HER2-enriched",
                                 cluster_red == 3 ~ "Basal-like",
                                 cluster_red == 4 ~ "Luminal B"),
         cluster_pca_correct = case_when(`PAM50 mRNA` == cluster_red ~ 1,
                                         `PAM50 mRNA` != cluster_red ~ 0)) %>% 
  summarise(score_pca = mean(cluster_pca_correct))

# Calculate BSS/TSS ratio ------------------------------------------------------
