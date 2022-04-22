# Load libraries ----------------------------------------------------------
library("tidyverse")
library('dplyr')
library("data.table")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
proteomes_raw <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")
patients_raw <- read_csv(file = "data/_raw/clinical_data_breast_cancer.csv")
PAM50_raw <- read_csv(file = "data/_raw/PAM50_proteins.csv")


# Wrangle data ------------------------------------------------------------
# Wrangle with proteomes_raw (modify format of column names)
proteomes <- proteomes_raw 

adjusted_names <- proteomes %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)

colnames(proteomes)[4:86] <- adjusted_names

# proteomes_long <- proteomes %>%
#   pivot_longer(-c(1:3), 
#                names_to = "TCGA_ID", 
#                values_to = "Expression_level")# %>% ...

# Wrangle with patients_raw (delete some of the columns? / modify colnames format)
patients <- patients_raw # %>% ...

# Wrangle with patients_raw
PAM50 <- PAM50_raw # %>% ...






# Write data --------------------------------------------------------------
write_csv(x = proteomes,
          file = "data/01_proteomes.csv")
write_csv(x = patients,
          file = "data/01_patients.csv")
write_csv(x = PAM50,
          file = "data/01_PAM50.csv")