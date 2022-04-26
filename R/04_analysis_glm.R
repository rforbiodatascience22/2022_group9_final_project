# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
Gene_Expresion_clean_aug <- read_csv(file = "data/03_Gene_Expresion_clean_aug.csv")
BC_data_clean_aug        <- read_csv(file = "data/03_BC_data_clean_aug.csv")


# Wrangle data ------------------------------------------------------------
# Add subtype to Gene Expression data
Expresion_Subtype <- left_join(Gene_Expresion_clean_aug,
                                    BC_data_clean_aug,
                                    by = c("TCGA ID" = "Complete TCGA ID")) %>%
  column_to_rownames(var = "TCGA ID") %>% 
  select(c(1:12553), "Subtype")


# Create long table
Expresion_Subtype_long <- Expresion_Subtype %>% 
  pivot_longer(cols = -Subtypes,
               names_to = "gene",
               values_to = "log2_expr_level")



# Model data
my_data_clean_aug %>% ...


# Visualise data ----------------------------------------------------------
my_data_clean_aug %>% ...


# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)