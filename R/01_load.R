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
# Wrangle with patients_raw  (modify first column)
patients <- patients_raw


# Wrangle with PAM50_raw
PAM50 <- PAM50_raw


# Wrangle with proteomes_raw (modify format of column names)
proteomes <- proteomes_raw 

adjusted_names <- proteomes %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes)[4:86] <- adjusted_names

# Transpose tibble
Gene_Expresion <- proteomes %>% #CANT BECAUSE NOT UNIQUE --> MODIFY BEFORE THIS
  select(c(1,4:86)) %>% 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


# Merge data --------------------------------------------------------------
BC_Data <- left_join(patients,    #WHAT JOIN
                  Gene_Expresion,
                  by = c("Complete TCGA ID", "TCGA ID",))


# Write data --------------------------------------------------------------
write_csv(x = proteomes,
          file = "data/01_proteomes.csv")
write_csv(x = patients,
          file = "data/01_patients.csv")
write_csv(x = PAM50,
          file = "data/01_PAM50.csv")
write_csv(x = BC_Data,
          file = "data/01_BC_Data.csv")
