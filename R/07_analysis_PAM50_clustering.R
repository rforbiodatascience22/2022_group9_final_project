library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
# library("tidymodels")
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


(plot_pca_ori_cumulative + plot_pca_red_cumulative) +
  plot_layout(guides = 'collect') &
  scale_x_continuous(expand = c(0, 0)) &
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.01))) &
  theme_half_open(12) &
  theme(text = element_text(family = "Avenir",
                            size = 12))


# K-means analysis (original data) ---------------------------------------------
set.seed(4)

k_org <- pca_ori_aug %>%
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = 4)

pca_aug_k_org <- k_org %>%
  augment(pca_ori_aug) %>% 
  rename(cluster_org = .cluster)
pca_aug_k_org


pl1 <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

pl2 <- pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")

pl1 + pl2


# K-means analysis (Reduced data) ----------------------------------------------
set.seed(4)
k_red <- pca_red_aug %>%
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = 4)

pca_aug_k_red <- k_red %>%
  augment(pca_red_aug) %>% 
  rename(cluster_org = .cluster)


pl3 <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

pl4 <- pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")

pl3 + pl4


# Calculate BSS/TSS ratio ------------------------------------------------------
k_red$betweenss/k_red$totss
