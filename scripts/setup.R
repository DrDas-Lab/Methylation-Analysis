packages <- c(
  "data.table",
  "tidyverse",
  "pheatmap",
  "ggplot2",
  "ggridges",
  "Seurat",
  "e1071",
  "preprocessCore",
  "matrixStats",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}