# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("broom")
library("purrr")
library("vroom")
rm(list=ls())


# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data --------------------------------------------------------------------
proteomes_clean <- read_csv(file = "data/02_proteomes_clean.csv")
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")



# Model data -------------------------------------------------------------------
# Making four different logistic models for each of the subtypes
# Luminal A
LumA_glm <- subtype_glm("Luminal_A", BC_data_clean_aug)

# Luminal B
LumB_glm <- subtype_glm("Luminal_B", BC_data_clean_aug)

# HER2_enriched
Her2_glm <- subtype_glm("HER2_enriched", BC_data_clean_aug)

# Basal-like
Basal_glm <- subtype_glm("Basal_like", BC_data_clean_aug)


# Visualise data ---------------------------------------------------------------
# Density plot of Luminal A
LumA_glm %>%
  ggplot(aes(estimate,
             fill = identified_as)) +
  geom_density(alpha = 0.6) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom") # Missing labels etc.

# Finding only significant proteomes
#BC_data_clean_aug_LuminalA %>% filter(identified_as_LuminalA == "significant")
# Reduction: 9274 to 1,861


# Density plot of Luminal B
LumB_glm %>%
  ggplot(aes(estimate,
             fill = identified_as)) +
  geom_density(alpha = 0.6) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom")

# BC_data_clean_aug_LuminalB %>% filter(identified_as_LuminalB == "significant")
# 9,274 to 1,435


# Density plot of HER2-enriched
Her2_glm %>%
  ggplot(aes(estimate,
             fill = identified_as)) +
  geom_density(alpha = 0.6) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom")


# Density plot of Basal-like
Basal_glm %>%
  ggplot(aes(estimate,
             fill = identified_as)) +
  geom_density(alpha = 0.6) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom")
# 9,274 to 2,503



# Visualize the model
Expresion_Subtype_long_nested %>% 
  filter(identified_as == "significant") %>% 
  ggplot(aes(x = estimate,
             y = fct_reorder(protein, desc(estimate)),
             colour = identified_as)) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low,
                     xmax = conf.high,
                     height = 0.2)) +
  theme_classic(base_family = "Avenir",
                base_size = 8) +
  theme(axis.text.y = element_text(),
        legend.position = "bottom") +
  labs(y = "") 
### We need to think of another way of visualization


# Write data ---------------------------------------------------------------
write_csv(LumA_glm, 
          file = "results/05_LumA_glm.csv")

write_csv(LumB_glm, 
          file = "results/05_LumB_glm.csv")

write_csv(Her2_glm, 
          file = "results/05_Her2_glm.csv")

write_csv(Basal_glm, 
          file = "results/05_Basal_glm.csv")

ggsave(...)