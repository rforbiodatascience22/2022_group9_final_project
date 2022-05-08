# Load libraries ----------------------------------------------------------
library("tidyverse")
rm(list=ls())

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
BC_data_clean   <- read_csv(file = "data/02_BC_data_clean.csv")

# Wrangle data ------------------------------------------------------------

# Modifying BC_data -Ccolumn: PAM50 mRNA  ##DOESNT WORK --> CREATES NA's###
BC_data_clean_aug <- BC_data_clean %>%
  mutate(Luminal_A = if_else(`PAM50 mRNA` == "Luminal A", 1, 0),
         Luminal_B = if_else(`PAM50 mRNA` == "Luminal B", 1, 0),
         HER2_enriched = if_else(`PAM50 mRNA` == "HER2-enriched", 1, 0),
         Basal_like = if_else(`PAM50 mRNA` == "Basal-like", 1, 0))

# Write data --------------------------------------------------------------
write_csv(x = BC_data_clean_aug,
          file = "data/03_BC_data_clean_aug.csv")
