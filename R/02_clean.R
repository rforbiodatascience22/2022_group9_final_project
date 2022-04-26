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

# Only selecting the proteomes that have XX% of data
proteomes_clean %>%
  mutate(frac_na = apply(., 1, count_na_func)/ncol(.)) %>% 
  filter(frac_na < 0.10)

# Modification of column names in proteomes
adjusted_names <- proteomes_clean_aug %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes_clean_aug)[4:86] <- adjusted_names


# Creating new tibble that is a transposed and reduced version of proteomes_clean_aug
Gene_Expresion_clean_aug <- proteomes_clean_aug %>%
  select(-c(2,3,13,71,77)) %>%  #Deletes gene_symbol, gene_name and the 3 duplicates 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


# Merge data --------------------------------------------------------------
BC_data_clean_aug <- left_join(patients_clean_aug,                #WHAT JOIN
                               Gene_Expresion_clean_aug,
                               by = c("Complete TCGA ID" = "TCGA ID"))


# Wrangle data ------------------------------------------------------------
BC_data_clean_aug <- BC_data_clean_aug %>%
  mutate(Subtype = case_when('PAM50 mRNA' == "Luminal A" ~ 0,
                             'PAM50 mRNA' == "Luminal B" ~ 1,
                             'PAM50 mRNA' == "HER2-enriched" ~ 2,
                             'PAM50 mRNA' == "Basal-like" ~ 3))



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
