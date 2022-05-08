# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/01_PlotData.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/05_analysis_glm.R")
source(file = "R/06_analysis_PickGenes.R")

# Add all files
# Make sure to include a programmatic call to knitr at the end of your doit-script
