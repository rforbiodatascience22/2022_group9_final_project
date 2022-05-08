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
p_AVC1 <- ggplot(data = patients,
                 mapping = aes(x = `PAM50 mRNA`,
                               y = `Age at Initial Pathologic Diagnosis`,
                               fill = `Vital Status`)) +
  geom_boxplot(alpha = 0.5) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = element_blank())

p_AVC2 <- ggplot(data = patients,
       mapping = aes(x = `PAM50 mRNA`)) +
  geom_bar(mapping = aes(fill = `Vital Status`),
           width = 0.6,
           alpha = 0.5) +
  theme(legend.position="right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  labs(x = element_blank())

(p_AVC1 + p_AVC2) +
  plot_annotation(title = "Publicly available proteomic Breast Cancer data",
                  theme = theme(plot.title = element_text(hjust = 0.5,
                                                          size = 16)))

# Save plot ---------------------------------------------------------------
ggsave(file = "results/04_plot_AgeVitalCancerType.png",
       width = 9, 
       height = 4, 
       dpi = 150)
