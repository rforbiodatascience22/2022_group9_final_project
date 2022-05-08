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