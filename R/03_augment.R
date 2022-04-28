# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients_clean  <- read_csv(file = "data/02_patients_clean.csv")
PAM50_clean     <- read_csv(file = "data/02_PAM50_clean.csv")
proteomes_clean <- read_csv(file = "data/02_proteomes_clean.csv")
BC_data_clean   <- read_csv(file = "data/02_BC_data_clean.csv")


# Wrangle data ------------------------------------------------------------
patients_clean_aug  <- patients_clean
PAM50_clean_aug     <- PAM50_clean
proteomes_clean_aug <- proteomes_clean


# Delete 3 normal data and transform 'PAM50 mRNA' column
BC_data_clean_aug <- BC_data_clean %>%
  #filter(str_length(`Complete TCGA ID`) == 12) %>% 
  mutate(Subtype = case_when(`PAM50 mRNA` == "Luminal A" ~ 0,
                             `PAM50 mRNA` == "Luminal B" ~ 1,
                             `PAM50 mRNA` == "HER2-enriched" ~ 2,
                             `PAM50 mRNA` == "Basal-like" ~ 3))


view(BC_data_clean_aug)
BC_data_clean %>%
  str_length(`Complete TCGA ID`)

# Write data --------------------------------------------------------------
write_csv(x = patients_clean_aug,
          file = "data/03_patients_clean_aug.csv")
write_csv(x = PAM50_clean_aug,
          file = "data/03_PAM50_clean_aug.csv")
write_csv(x = proteomes_clean_aug,
          file = "data/03_proteomes_clean_aug.csv")
write_csv(x = BC_data_clean_aug,
          file = "data/03_BC_data_clean_aug.csv")
