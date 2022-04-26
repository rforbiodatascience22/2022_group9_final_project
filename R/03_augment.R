# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients <- read_csv(file = "data/01_patients.csv")
PAM50 <- read_csv(file = "data/01_PAM50.csv")
proteomes <- read_csv(file = "data/01_proteomes.csv") 

# Wrangle data ------------------------------------------------------------
# Modification of column names in proteomes
adjusted_names <- proteomes %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)
colnames(proteomes)[4:86] <- adjusted_names

# Creating new tibble that is a transposed reduced version of proteomes
Gene_Expresion <- proteomes %>%
  select(-c(2,3,13,71,77)) %>%  #Deletes gene_symbol, gene_name and the 3 duplicates 
  pivot_longer(cols= -1,
               names_repair = "check_unique") %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)


# Merge data --------------------------------------------------------------
BC_Data <- left_join(patients,                #WHAT JOIN
                     Gene_Expresion,
                     by = c("Complete TCGA ID" = "TCGA ID"))


# Write data --------------------------------------------------------------
write_csv(x = BC_Data,
          file = "data/01_BC_Data.csv")
