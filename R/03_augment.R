# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients_clean  <- read_csv(file = "data/02_patients_clean.csv")
PAM50_clean     <- read_csv(file = "data/02_PAM50_clean.csv")
proteomes_clean <- read_csv(file = "data/02_proteomes_clean.csv")


# Wrangle data ------------------------------------------------------------
patients_clean_aug  <- patients_clean
PAM50_clean_aug     <- PAM50_clean
proteomes_clean_aug <- proteomes_clean



# Write data --------------------------------------------------------------
write_csv(x = patients_clean_aug,
          file = "data/03_patients_clean_aug.csv")
write_csv(x = PAM50_clean_aug,
          file = "data/03_PAM50_clean_aug.csv")
write_csv(x = proteomes_clean_aug,
          file = "data/03_proteomes_clean_aug.csv")
write_csv(x = Gene_Expresion_clean_aug,
          file = "data/03_Gene_Expresion_clean_aug.csv")
write_csv(x = BC_data_clean_aug,
          file = "data/03_BC_data_clean_aug.csv")