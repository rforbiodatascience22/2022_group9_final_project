# Load libraries ----------------------------------------------------------
library(tidyverse)
library(readr)
library(ggplot2)
library(patchwork)
library(ggthemes)
rm(list=ls())

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
patients  <- read_csv(file = "data/01_patients.csv")


# visualise data ----------------------------------------------------------

# Box plot of the age of diagnosis for the different subtypes 
p_AVC1 <- ggplot(data = patients,
                 mapping = aes(x = `PAM50 mRNA`,
                               y = `Age at initial pathologic diagnosis`,
                               fill = `Vital Status`)) +
  geom_boxplot(alpha = 0.5) +
  expand_limits(y = c(0, NA)) + 
  new_theme + 
  theme(legend.position = "none") +
  labs(x = element_blank())
  

# Bar plot of vital state of breast cancer patients for different subtypes
p_AVC2 <- ggplot(data = patients,
       mapping = aes(x = `PAM50 mRNA`)) +
  geom_bar(mapping = aes(fill = `Vital Status`),
           width = 0.6,
           alpha = 0.5) +
  labs(y = "Number of cases") + 
  new_theme +
  theme(legend.position="right") + 
  labs(x = element_blank())
  

# Plot of the two combined
(p_AVC1 + p_AVC2) +
  plot_annotation(title = "Breast cancer data for different subtypes",
                  caption = "SOURCE DATA",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))
#NB: ADD SOURCE DATA

# Save plot ---------------------------------------------------------------
ggsave(file = "results/04_plot_AgeVitalCancerType.png",
       width = 12, 
       height = 5, 
       dpi = 150)
