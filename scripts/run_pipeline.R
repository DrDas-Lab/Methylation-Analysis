############################################################
## AML MULTIOMICS TRAJECTORY PIPELINE
## Author: Your Name
## Description:
## Bulk RNA + DNA methylation integration in AML
############################################################

rm(list = ls())

# -----------------------------
# Load package manager
# -----------------------------
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::activate()

source("setup.R")

# -----------------------------
# Run analysis modules
# -----------------------------
source("scripts/01_preprocess_counts.R")
source("scripts/02_CIBERSORT_pipeline.R")
source("scripts/03_ModuleScore_trajectory.R")
source("scripts/04_RNA_heatmap.R")
source("scripts/05_Methylation_processing.R")
source("scripts/06_Correlation_plots.R")

cat("Pipeline completed successfully\n")
