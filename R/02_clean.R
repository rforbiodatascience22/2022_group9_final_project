# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
BC_Data <- read_csv(file = "data/01_BC_Data.csv")


# Wrangle data ------------------------------------------------------------
BC_Data <- my_data # %>% ...


#library(dplyr)

# Which threshold should we set?
# if I just write is.na(), it gives me na errors

# Only using those with a fraction of NAs less than XX%
BC_Data %>%
  mutate(frac_NA = rowSums(is.na(BC_Data), na.rm = TRUE)/ncol(BC_Data[1,])) %>% 
  filter(frac_NA < 0.20)

# Columns:
BC_Data %>%
  mutate(frac_NA = colSums(is.na(BC_Data), na.rm = TRUE)/nrow(BC_Data)) %>% 
  filter(frac_NA < 0.20)


colSums(is.na(BC_Data[, -c(1:32)]), na.rm=TRUE)
colSums(is.na(BC_Data[, -c(1:32)]), na.rm=TRUE)/nrow(BC_Data)

BC_Data %>%
  filter(map(~sum(is.na(.))/nrow(BC_Data) < 0.1))

# sum(is.na(BC_Data[1,]))/ncol(BC_Data[1,])

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")

# "data/01_proteomes.csv"
# "data/01_patients.csv"
# "data/01_PAM50.csv"
# "data/01_BC_Data.csv"
