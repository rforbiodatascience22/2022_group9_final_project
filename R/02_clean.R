# Load libraries ----------------------------------------------------------
library("tidyverse")
library("dplyr")
rm(list=ls())

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients  <- read_csv(file = "data/01_patients.csv")
proteomes <- read_csv(file = "data/01_proteomes.csv")
PAM50 <- read_csv(file = "data/01_PAM50.csv")

# Wrangle data ------------------------------------------------------------

# For patients: 
# Reduce meaningless variables
patients_clean <- patients %>% 
  select(-starts_with("Days"))

# For proteomes: 
# Create new column names
adjusted_names <- proteomes %>% 
  select(-c("RefSeq_accession_number",
            "gene_symbol",
            "gene_name")) %>% 
  colnames() %>% 
  map(change_format) %>% 
  unlist() %>% 
  unique()


# For PAM50
PAM50_clean <- PAM50 %>%
  pivot_longer(cols = -RefSeqProteinID) %>%
  pivot_wider(names_from = RefSeqProteinID) %>%
  rename("RefSeqProteinID" = 1) %>%
  as_tibble()


# Wrangle proteomes (removing duplicates, renaming proteins)
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


# Merge data --------------------------------------------------------------
BC_data_clean <- inner_join(patients_clean,           
                            proteomes_clean,
                            by = c("Complete TCGA ID" = "TCGA ID"))


# Subset data ------------------------------------------------------------------
BC_data_PAM50_clean <- BC_data_clean %>% 
  select(c(1:28),
         names(BC_data_clean)[names(BC_data_clean) %in% names(PAM50_clean)])


# Write data --------------------------------------------------------------
write_csv(x = patients_clean,
          file = "data/02_patients_clean.csv")
write_csv(x = proteomes_clean,
          file = "data/02_proteomes_clean.csv")
write_csv(x = BC_data_clean,
          file = "data/02_BC_data_clean.csv")
write_csv(x = BC_data_PAM50_clean,
          file = "data/02_BC_data_PAM50_clean.csv")
