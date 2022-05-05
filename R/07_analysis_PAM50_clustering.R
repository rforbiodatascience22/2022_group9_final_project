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


pca_fit <- BC_data_clean_aug %>% 
  select(c(29:8022)) %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = TRUE)



pca_fit %>%
  augment(BC_data_clean_aug) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = `PAM50 mRNA`)) + 
  geom_point(size = 1.5) +
  theme_half_open(12) + background_grid()


pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:9) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal_hgrid(12)






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