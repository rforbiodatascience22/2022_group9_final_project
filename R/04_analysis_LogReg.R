# Load libraries ----------------------------------------------------------
library("tidyverse")
library('readr')
library('ggplot2')


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
BC_Data <- read_csv(file = "data/01_BC_Data.csv")


# Wrangle data ------------------------------------------------------------
my_data_clean_aug %>% ...


# Model data
my_data_clean_aug %>% ...


# Visualise data ----------------------------------------------------------

# Correlation between differnent variables
ggplot(data = BC_Data,
       mapping = aes(y = "Age at Initial Pathologic Diagnosis")) +
  geom_boxplot()


# Correlation between differnent variables
BC_Data %>% 
  select("Age at Initial Pathologic Diagnosis",`Integrated Clusters (with PAM50)`,`miRNA Clusters`) %>% 
  pairs()


# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)