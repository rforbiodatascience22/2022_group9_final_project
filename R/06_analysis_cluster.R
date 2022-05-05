library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
library("tidymodels")
library(cowplot)

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
PAM50_data <- read_csv(file = "data/01_PAM50.csv")

# Modifying data files ---------------------------------------------------------

# Transforming PAM50 data (NB: SHOULD BOTH BE MOVE TO CLEANING PART?)
#PAM50_selected <- PAM50_data %>%
#  pivot_longer(cols = -RefSeqProteinID) %>%
#  pivot_wider(names_from = RefSeqProteinID) %>%
#  rename("RefSeqProteinID" = 1) %>%
#  as_tibble()


#BC_data_clean_aug_PAM50 <- BC_data_clean_aug %>% 
#  select(c(1:28),
#         names(BC_data_clean_aug)[names(BC_data_clean_aug) %in% names(PAM50_selected)])

# Selecing the ones that are the same for both


# PCA analysis  ----------------------------------------------------------------

# Doing a PCA analysis - only numerical data for protein IDs
pca_fit_ori <- BC_data_clean_aug %>% 
  select(c(29:8022)) %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_fit_ori_aug <- pca_fit_ori %>% 
  augment(BC_data_clean_aug)

# PCA on reduced version
pca_fit_reduced <- BC_data_clean_aug_PAM50 %>% 
  select(c(29:54)) %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(center = TRUE, scale. = TRUE)

# Adding original data back
pca_fit_red_aug <- pca_fit_reduced %>% 
  augment(BC_data_clean_aug_PAM50)

# Visualizing the results  -----------------------------------------------------

# Scatter plot of 1st and 2nd PC - Original data
pca_fit_original_plot <- pca_fit_ori_aug %>% # add original dataset back in
  ggplot(aes(.fittedPC1, 
             .fittedPC2, 
             color = `PAM50 mRNA`)) + 
  geom_point(size = 1.5) +
  theme_half_open(12) + 
  background_grid()

# Scatter plot of 1st and 2nd PC - Reduced data
pca_fit_reduced_plot <- pca_fit_red_aug %>%
  ggplot(aes(.fittedPC1, 
             .fittedPC2, 
             color = `PAM50 mRNA`)) + 
  geom_point(size = 1.5) +
  theme_half_open(12) + 
  background_grid()

# Variance explained by each PC (plots) - Original data
pca_fit_original_PC_plot <- pca_fit_original %>%
  tidy("pcs") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))) +
  theme_half_open(12)


# Variance explained by each PC (plots) - Reduced data
pca_fit_reduced_PC_plot <- pca_fit_reduced %>%
  tidy("pcs") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))) +
  theme_half_open(12)

# Plots all together
(pca_fit_original_plot + pca_fit_reduced_plot +  plot_layout(guides = 'auto'))/
  (pca_fit_original_PC_plot + pca_fit_reduced_PC_plot) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "top")


# K-means analysis (original data) ---------------------------------------------
set.seed(4)

my_k_org <- pca_fit_ori_aug %>%
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = 4)

my_pca_aug_k_org <- my_k_org %>%
  augment(pca_fit_ori_aug) %>% 
  rename(cluster_org = .cluster)
my_pca_aug_k_org


pl1 <- my_pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

pl2 <- my_pca_aug_k_org %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")

pl1 + pl2


# K-means analysis (Reduced data) ----------------------------------------------
set.seed(4)
my_k_red <- pca_fit_red_aug %>%
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = 4)

my_pca_aug_k_red <- my_k_red %>%
  augment(pca_fit_red_aug) %>% 
  rename(cluster_org = .cluster)
my_pca_aug_k_red


pl3 <- my_pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

pl4 <- my_pca_aug_k_red %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")

pl3 + pl4


# TRASH ?? ---------------------------------------------------------------------

my_k_red <- pca_fit_reduced %>%
  augment(BC_data_clean_aug_PAM50) %>% 
  select(.fittedPC1, .fittedPC2) %>%
  kmeans(centers = 4)


pl1 <- my_k_org %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = `PAM50 mRNA`)) +
  geom_point() +
  theme(legend.position = "bottom")

pl2 <- my_pca_aug_k_org_pca %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = cluster_org)) +
  geom_point() +
  theme(legend.position = "bottom")



points <- BC_data_clean_aug %>% 
  select(c(29:8022))

k_means <- BC_data_clean_aug %>% 
  select(c(29:8022)) %>% 
  na.omit() %>% 
  kmeans(centers = 4)




#augment(k_means, BC_data_clean_aug)

#tidy(k_means)

#glance(k_means)
  
kclusts <- 
  tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(points, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, BC_data_clean_aug)
  )

kclusts

clusters <- 
  kclusts %>%
  unnest(cols = c(tidied))

assignments <- 
  kclusts %>% 
  unnest(cols = c(augmented))

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))



p1 <- 
  ggplot(assignments, aes(x = NP_958782, y = NP_000436)) +
  geom_point(aes(color = .cluster), alpha = 0.8) + 
  facet_wrap(~ k)
p1