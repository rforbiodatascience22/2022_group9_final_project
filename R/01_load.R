# Load libraries ----------------------------------------------------------
library("tidyverse")
library("data.table") # check if we use this
library("readr")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
proteomes_raw <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")
patients_raw  <- read_csv(file = "data/_raw/clinical_data_breast_cancer.csv")
PAM50_raw     <- read_csv(file = "data/_raw/PAM50_proteins.csv")


# Wrangle data ------------------------------------------------------------
patients  <- patients_raw
PAM50     <- PAM50_raw
proteomes <- proteomes_raw 


# Write data --------------------------------------------------------------
write_csv(x = patients,
          file = "data/01_patients.csv")
write_csv(x = PAM50,
          file = "data/01_PAM50.csv")
write_csv(x = proteomes,
          file = "data/01_proteomes.csv")

