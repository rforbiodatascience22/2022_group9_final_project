# Load libraries ----------------------------------------------------------
library("tidyverse")
library("dplyr")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
patients  <- read_csv(file = "data/01_patients.csv")
PAM50     <- read_csv(file = "data/01_PAM50.csv")
proteomes <- read_csv(file = "data/01_proteomes.csv") 


# Wrangle data ------------------------------------------------------------
patients_clean  <- patients
PAM50_clean     <- PAM50
proteomes_clean <- proteomes


#count_na_func <- function(x) sum(is.na(x)) 

# Modification of column names in proteomes so they are comparable with patients_clean_aug
adjusted_names <- proteomes_clean %>% 
  select(-c("RefSeq_accession_number","gene_symbol","gene_name")) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes_clean)[4:86] <- adjusted_names


# Creating new data sets with columns consisting of dublicates and too little data removed:
proteomes_clean <- proteomes_clean %>% 
  select(unique(colnames(.))) %>%                          # Removing duplicates
  select(-c("gene_symbol","gene_name")) %>%                # Removing unnecessary describtions of protein
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10) %>%                               # Removing columns consisting of more than 10% NAs
  select(-c("frac_na")) %>% 
  pivot_longer(cols= -1,                                   # Transposing
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


# Merge data --------------------------------------------------------------
# SHOULD
BC_data_clean <- right_join(patients_clean,                #WHAT JOIN
                            proteomes_clean,
                            by = c("Complete TCGA ID" = "TCGA ID"))


# Write data --------------------------------------------------------------
write_csv(x = patients_clean,
          file = "data/02_patients_clean.csv")
write_csv(x = PAM50_clean,
          file = "data/02_PAM50_clean.csv")
write_csv(x = proteomes_clean,
          file = "data/02_proteomes_clean.csv")
write_csv(x = BC_data_clean,
          file = "data/02_BC_data_clean.csv")

# "data/01_proteomes.csv"
# "data/01_patients.csv"
# "data/01_PAM50.csv"
# "data/01_BC_Data.csv"
