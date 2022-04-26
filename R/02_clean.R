# Load libraries ----------------------------------------------------------
library("tidyverse")
library("dplyr")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
patients <- read_csv(file = "data/01_patients.csv")
PAM50 <- read_csv(file = "data/01_PAM50.csv")
proteomes <- read_csv(file = "data/01_proteomes.csv") 


# Wrangle data ------------------------------------------------------------
patients_clean <- patients


PAM50_clean <- PAM50


#count_na_func <- function(x) sum(is.na(x)) 

proteomes_clean %>%
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10)

# Only selecting the proteomes that have XX% of data

# Write data --------------------------------------------------------------
write_csv(x = patients_clean,
          file = "data/02_patients_clean.csv")
write_csv(x = PAM50_clean,
          file = "data/02_PAM50_clean.csv")
write_csv(x = proteomes_clean,
          file = "data/02_proteomes_clean.csv")

# "data/01_proteomes.csv"
# "data/01_patients.csv"
# "data/01_PAM50.csv"
# "data/01_BC_Data.csv"
