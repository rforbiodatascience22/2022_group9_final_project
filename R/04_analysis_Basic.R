# Load libraries ----------------------------------------------------------
library("tidyverse")
library('readr')
library('ggplot2')

library(patchwork)


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
proteomes_clean_aug <- read_csv(file = "data/03_proteomes_clean_aug.csv")
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")


# Varibles Tables ------------------------------------------------------------

BC_variable <- BC_data_clean_aug %>%
  select(-starts_with(c("NP","XP","YP"))) %>% 
  select(-Subtype) %>% 
  colnames() %>% 
  c(.,"NP_XXXXXX")

BC_variable_explinaion <- c("Patient ID (80 patients)",
                 "(All Female)",
                 "Age When Diagnosed (30-88 years)",
                 "ER Status",
                 "PR Status",                          
                 "HER2 Final Status",                  
                 "Tumor",                 
                 "Tumor--T1 Coded",                
                 "Node",               
                 "Node-Coded",               
                 "Metastasis",              
                 "Metastasis-Coded",             
                 "AJCC Stage",            
                 "Converted Stage",           
                 "Survival Data Form",          
                 "Vital Status",         
                 "Days to Date of Last Contact",       
                 "Days to date of Death",       
                 "OS event",      
                 "OS Time",     
                 "PAM50 mRNA",    
                 "SigClust Unsupervised mRNA",   
                 "SigClust Intrinsic mRNA",  
                 "miRNA Clusters", 
                 "methylation Clusters",
                 "RPPA Clusters",
                 "CN Clusters",
                 "Integrated Clusters (with PAM50)",
                 "Integrated Clusters (no exp)",
                 "Integrated Clusters (unsup exp)",
                 "NP_958782")

BC_var_table <- tibble(BC_variable,BC_variable_explinaion)
print(BC_var_table, n=31)





# Wrangle data ------------------------------------------------------------
my_data_clean_aug %>% ...


# Model data
my_data_clean_aug %>% ...


# Visualise data ----------------------------------------------------------

# PLOT: Cancer type, vital status and Age
p_AVC1 <- ggplot(data = BC_data_clean_aug,
                 mapping = aes(x = `PAM50 mRNA`,
                               y = `Age at Initial Pathologic Diagnosis`,
                               fill = `Vital Status`)) +
  geom_boxplot(alpha = 0.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = "Breast Cancer Type")

p_AVC2 <- ggplot(data = BC_data_clean_aug,
       mapping = aes(x = `PAM50 mRNA`)) +
  geom_bar(mapping = aes(fill = `Vital Status`),
           width = 0.6,
           alpha = 0.5) +
  theme(legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  labs(x = "Breast Cancer Type")

Plot_AgeVitalCtype <- (p_AVC1 + p_AVC2) +
  plot_annotation(title = "Distribution of Age and Vital Status for each observed Breast Cancer Type",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

ggsave(file = "results/04_plot_AgeVitalCancerType.png",
       width = 8.96, 
       height = 5.42, 
       dpi = 150)

# Write data --------------------------------------------------------------
write_tsv(...)





readPNG(system.file("results/04_plot_AgeVitalCancerType.rda", package="png"), TRUE)
t1 <- num_table
