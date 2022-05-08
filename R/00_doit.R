# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_PlotData.R")
source(file = "R/05_analysis_PickGenes.R")
source(file = "R/06_analysis_CommonGene_clustering.R")
source(file = "R/07_analysis_PAM50_clustering.R")

# Knit presentation
rmarkdown::render(input = "doc/presentation.Rmd",
                  output_file="../doc/group9_presentation.html", 
                  knit_root_dir = "..")
