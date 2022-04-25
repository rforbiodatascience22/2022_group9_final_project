# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
BC_data_clean <- read_tsv(file = "data/02_BC_Data.csv")


# Wrangle data ------------------------------------------------------------
BC_data_clean_aug <- BC_data_clean %>%
  mutate(Subtype = case_when('PAM50 mRNA' == "Luminal A" ~ 0,
                             'PAM50 mRNA' == "Luminal B" ~ 1,
                             'PAM50 mRNA' == "HER2-enriched" ~ 2,
                             'PAM50 mRNA' == "Basal-like" ~ 3))


# Write data --------------------------------------------------------------
write_csv(x = BC_data_clean_aug,
          file = "data/03_BC_data_clean_aug.csv")