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


# Wrangle data (1 METHOD) ------------------------------------------------------------
# Wrangle with proteomes_raw (modify format of column names)
proteomes <- proteomes_raw 

adjusted_names <- proteomes %>% 
  select(4:86) %>% 
  colnames() %>% 
  map(change_format)

colnames(proteomes)[4:86] <- adjusted_names

# proteomes_long <- proteomes %>%
#   pivot_longer(-c(1:3), 
#                names_to = "TCGA_ID", 
#                values_to = "Expression_level")# %>% ...

# Wrangle with patients_raw (delete some of the columns? / modify colnames format)
patients <- patients_raw # %>% ...

# Wrangle with patients_raw
PAM50 <- PAM50_raw # %>% ...



# Wrangle data (2 METHOD) ------------------------------------------------------------
# Transpose tibble
Gene_Expresion <- proteomes_raw[,c(1,4:86)] %>% 
  pivot_longer(cols= -1) %>% 
  pivot_wider(names_from = "RefSeq_accession_number",
              values_from = "value") %>% 
  rename("TCGA ID" = name)

# Adds "ID_short" to make merge easier
patients <- patients_raw %>% 
  mutate("ID_short" = substr(patients_raw$`Complete TCGA ID`,6,nchar(patients_raw$`Complete TCGA ID`))) #ALL unique (105)
Gene_Expresion <- Gene_Expresion %>% 
  mutate("ID_short" = substr(Gene_Expresion$`TCGA ID`,0,7))                                             #80/83 are unique (TAKE MEAN)

#Which ones are there several of in proteomes?
table(substr(Gene_Expresion$`TCGA ID`,0,7))

# Merge data
my_data <- full_join(patients,
                  Gene_Expresion,
                  by = "ID_short")



Gene_Expresion[,1:200] %>% 
  select(-c(`TCGA ID`)) %>%
  cor(use="complete.obs")

dim(Gene_Expresion)
  
str(Gene_Expresion)

cor(proteomes_raw$`C8-A131.32TCGA`,proteomes_raw$`C8-A131.01TCGA`,use="complete.obs")
#transpose t() if you want to check correlation between patients --> to check the 3 not unique in proteomes



# Write data --------------------------------------------------------------
write_csv(x = proteomes,
          file = "data/01_proteomes.csv")
write_csv(x = patients,
          file = "data/01_patients.csv")
write_csv(x = PAM50,
          file = "data/01_PAM50.csv")