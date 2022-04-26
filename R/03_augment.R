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


# Modification of column names in proteomes
adjusted_names <- proteomes_clean_aug %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes_clean_aug)[4:86] <- adjusted_names


# Creating new tibble that is a transposed and reduced version of proteomes_clean_aug
Gene_Expresion <- proteomes_clean_aug %>%
  select(-c(2,3,13,71,77)) %>%  #Deletes gene_symbol, gene_name and the 3 duplicates 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


# Merge data --------------------------------------------------------------
BC_Data <- left_join(patients_clean_aug,                #WHAT JOIN
                     Gene_Expresion,
                     by = c("Complete TCGA ID" = "TCGA ID"))


# Wrangle data ------------------------------------------------------------
BC_data_clean_aug <- BC_data_clean %>%
  mutate(Subtype = case_when('PAM50 mRNA' == "Luminal A" ~ 0,
                             'PAM50 mRNA' == "Luminal B" ~ 1,
                             'PAM50 mRNA' == "HER2-enriched" ~ 2,
                             'PAM50 mRNA' == "Basal-like" ~ 3))


# Write data --------------------------------------------------------------
write_csv(x = BC_Data,
          file = "data/01_BC_Data.csv")
#BC_data_clean <- read_csv(file = "data/02_BC_Data.csv")
write_csv(x = BC_data_clean_aug,
          file = "data/03_BC_data_clean_aug.csv")