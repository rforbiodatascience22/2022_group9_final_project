# Load libraries ----------------------------------------------------------
library("tidyverse")
library("readr")
rm(list=ls())

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
PAM50     <- read_csv(file = "data/01_PAM50.csv")
patients  <- read_csv(file = "data/01_patients.csv")
proteomes <- read_csv(file = "data/01_proteomes.csv")


# Wrangle data ------------------------------------------------------------

# PAM50 - Transposing data
PAM50_clean <- PAM50 %>%
  pivot_longer(cols = -RefSeqProteinID) %>%
  pivot_wider(names_from = RefSeqProteinID) %>%
  rename("RefSeqProteinID" = 1) %>%
  as_tibble()


# For proteomes - Re-naming column names
adjusted_names <- proteomes %>% 
  select(-c("RefSeq_accession_number",
            "gene_symbol",
            "gene_name")) %>% 
  colnames() %>% 
  map(change_format) %>% 
  unlist() %>% 
  unique()

# Remove duplicates and rename proteins
proteomes_clean <- proteomes %>% 
  select(-c("AO-A12D.05TCGA",
            "C8-A131.32TCGA",
            "AO-A12B.34TCGA")) %>%         
  rename_with(~ adjusted_names,
              .cols = -c(1:3)) %>%        
  select(-c("gene_symbol",
            "gene_name")) %>%              
  mutate(frac_na = apply(.,
                         1,
                         count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.01) %>%               
  select(-c("frac_na")) %>% 
  pivot_longer(cols= -1,                            
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>%  
  rename("TCGA ID" = name)


# Merge data -------------------------------------------------------------------
BC_data_clean <- inner_join(patients,           
                            proteomes_clean,
                            by = c("Complete TCGA ID" = "TCGA ID"))


# Subset data (Common genes in PAM50 and BC_data) -------------------------
BC_data_PAM50_clean <- BC_data_clean %>% 
  select(-starts_with(c("NP","XP","YP")),
         colnames(BC_data_clean)[colnames(BC_data_clean) %in% colnames(PAM50_clean)])


# Write data --------------------------------------------------------------
write_csv(x = PAM50_clean,
          file = "data/02_PAM50_clean.csv")
write_csv(x = proteomes_clean,
          file = "data/02_proteomes_clean.csv")
write_csv(x = BC_data_clean,
          file = "data/02_BC_data_clean.csv")
write_csv(x = BC_data_PAM50_clean,
          file = "data/02_BC_data_PAM50_clean.csv")
