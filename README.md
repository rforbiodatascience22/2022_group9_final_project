# 2022 R Final Project
This repository is for the final project of course 22100 by Group 9

## Introduction 
This project aims to investigate whether 4 Breast Cancer subtypes can be detected and classified correctly by the gene expression Using the proteomic data. Analyzes used are GLM, PCA and K-means.


## Project Structure
**data:**
All data files, which are derivatives from the raw data generated by project scripts\
patients            <- read_csv(file = "data/01_patients.csv") \n
PAM50               <- read_csv(file = "data/01_PAM50.csv")             # not final file\
PAM50_clean         <- read_csv(file = "data/02_PAM50_clean.csv")
proteomes           <- read_csv(file = "data/01_proteomes.csv")         # not final file
proteomes_clean     <- read_csv(file = "data/02_proteomes_clean.csv")
BC_data_clean       <- read_csv(file = "data/02_BC_data_clean.csv")     # not final file
BC_data_clean_aug   <- read_csv(file = "data/03_BC_data_clean_aug.csv")
BC_data_PAM50_clean <- read_csv(file = "data/02_BC_data_PAM50_clean.csv")

**R:**
Include all project scripts\
file.edit("R/01_load.R")
file.edit("R/02_clean.R")
file.edit("R/03_augment.R")
file.edit("R/04_PlotData.R")
file.edit("R/05_analysis_glm.R")
file.edit("R/06_analysis_PickGenes.R")
file.edit("R/07_analysis_PAM50_clustering.R")

**results:**
library(magick)
image_read("results/04_plot_AgeVitalCancerType.png")
read_csv(file = "results/05_Basal_glm.csv")
read_csv(file = "results/05_Her2_glm.csv")
read_csv(file = "results/05_LumA_glm.csv")
read_csv(file = "results/05_LumB_glm.csv")
image_read("results/05_subtype_heatmap.png")
image_read("results/05_venndiagram.png")
read_csv(file = "results/06_BC_overlap_genes.csv")
image_read("results/06_BC_overlap_PCA.png")
image_read("results/07_BC_data_cumulative_kmeans.png")
image_read("results/07_BC_data_PAM50_cumulative_kmeans.png")

**doc:**
rmarkdown::render(input = "doc/presentation.Rmd",
                  output_file="../doc/group9_presentation.html", 
                  knit_root_dir = "..")


## Running the Scripts
source(file = "R/00_doit.R")
