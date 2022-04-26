# Load libraries ----------------------------------------------------------
library("tidyverse")
library('dplyr')
library("data.table")
library('readr')

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
proteomes_raw <- read_csv(file = "data/_raw/77_cancer_proteomes_CPTAC_itraq.csv")
patients_raw  <- read_csv(file = "data/_raw/clinical_data_breast_cancer.csv")
PAM50_raw     <- read_csv(file = "data/_raw/PAM50_proteins.csv")


# Wrangle data ------------------------------------------------------------
patients <- patients_raw
PAM50 <- PAM50_raw
proteomes <- proteomes_raw 

<<<<<<< HEAD
=======
#Modification of column names
adjusted_names <- proteomes %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes)[4:86] <- adjusted_names

# Transpose tibble and removal of duplicates
Gene_Expresion <- proteomes %>%
  select(-c(2,3,13,71,77)) %>%  #Deletes gene_symbol, gene_name and the 3 duplicates 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)

# How to find those with less than 20% 

#proteomes_raw %>% 
#  mutate(frac_NA = rowSums(is.na(.), na.rm = TRUE)/ncol(proteomes_raw)) %>% 
#  filter(frac_NA < 0.20)

# Merge data --------------------------------------------------------------
BC_Data <- left_join(patients,                #WHAT JOIN
                  Gene_Expresion,
                  by = c("Complete TCGA ID" = "TCGA ID"))


>>>>>>> f73974b21a1ad37754ee3f1ea1527bf705ba01a3
# Write data --------------------------------------------------------------
write_csv(x = patients,
          file = "data/01_patients.csv")
write_csv(x = PAM50,
          file = "data/01_PAM50.csv")
write_csv(x = proteomes,
          file = "data/01_proteomes.csv")
<<<<<<< HEAD
=======
write_csv(x = Gene_Expresion,
          file = "data/01_Gene_Expression.csv")
write_csv(x = BC_Data,
          file = "data/01_BC_Data.csv")

BC_Data
>>>>>>> f73974b21a1ad37754ee3f1ea1527bf705ba01a3
