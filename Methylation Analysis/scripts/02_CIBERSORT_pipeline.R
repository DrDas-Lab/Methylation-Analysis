############################################################
## 02_CIBERSORT_pipeline.R
## Purpose:
## Run AML cell-state deconvolution using CIBERSORT
############################################################

message("Running CIBERSORT AML pipeline")
graphics.off()

# ---------------------------------------------------------
# 1. Load required libraries
# ---------------------------------------------------------
# library(tidyverse)
# library(pheatmap)
# library(ggplot2)

# ---------------------------------------------------------
# 2. Source CIBERSORT function
# ---------------------------------------------------------
source("CIBERSORT.R")

# ---------------------------------------------------------
# 3. Load log2RPM expression
# ---------------------------------------------------------
expr_log <- read.table(
  "input/BAML_Normalized.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

expr_log <- as.matrix(expr_log)

# ---------------------------------------------------------
# 4. Convert log2RPM → linear RPM
# ---------------------------------------------------------
expr_linear <- (2^expr_log) - 1

# ---------------------------------------------------------
# 5. Collapse duplicated gene symbols
# ---------------------------------------------------------
expr_linear <- as.data.frame(expr_linear)
expr_linear$Gene <- rownames(expr_linear)

expr_clean <- expr_linear %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), max)) %>%
  ungroup()

expr_clean <- as.data.frame(expr_clean)
rownames(expr_clean) <- expr_clean$Gene
expr_clean$Gene <- NULL

expr_clean <- as.matrix(expr_clean)

stopifnot(!any(duplicated(rownames(expr_clean))))

# ---------------------------------------------------------
# 6. Write mixture file
# ---------------------------------------------------------
write.table(
  expr_clean,
  "output/CIBERSORT_mixture_RPM.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ---------------------------------------------------------
# 7. Load AML signature matrix
# ---------------------------------------------------------
sig <- read.table(
  "input/AML_signature_matrix_M.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

sig <- as.matrix(sig)

# ---------------------------------------------------------
# 8. Gene overlap check
# ---------------------------------------------------------
overlap <- length(intersect(rownames(sig), rownames(expr_clean)))
cat("Gene overlap:", overlap, "\n")

if (overlap < 500) {
  stop("Low gene overlap between mixture and signature matrix")
}

# ---------------------------------------------------------
# 9. Run CIBERSORT
# ---------------------------------------------------------
cibersort_res <- CIBERSORT(
  sig_matrix   = "input/AML_signature_matrix_M.txt",
  mixture_file = "output/CIBERSORT_mixture_RPM.txt",
  perm         = 100,
  QN           = FALSE
)

# Save results
write.csv(
  cibersort_res,
  "output/CIBERSORT_AML_results.csv"
)

# ---------------------------------------------------------
# 10. Extract malignant AML states
# ---------------------------------------------------------
malignant_states <- c(
  "LSPC-Quiescent",
  "LSPC-Primed",
  "LSPC-Cycle",
  "GMP-like",
  "ProMono-like",
  "Mono-like",
  "cDC-like"
)

malignant <- cibersort_res[, malignant_states]

# ---------------------------------------------------------
# 11. Normalize malignant fractions
# ---------------------------------------------------------
malignant_norm <- malignant / rowSums(malignant)

write.csv(
  malignant_norm,
  "output/CIBERSORT_AML_malignant_normalized.csv"
)

# ---------------------------------------------------------
# 12. Dominant AML state per sample
# ---------------------------------------------------------
dominant_state <- apply(malignant_norm, 1, which.max)
dominant_state <- malignant_states[dominant_state]

dominant_df <- data.frame(
  Sample = rownames(malignant_norm),
  Dominant_Cell_State = dominant_state
)

write.csv(
  dominant_df,
  "output/AML_dominant_cell_state.csv",
  row.names = FALSE
)

message("Generating CIBERSORT Heatmap & PCA visualization")
# ---------------------------------------------------------
# 1. Heatmap of AML cell states
# ---------------------------------------------------------

pdf("plots/CIBERSORT_AML_heatmap.pdf", width = 8, height = 6)

pheatmap(
  t(malignant_norm),
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  fontsize = 10,
  main = "AML Cell-State Composition (CIBERSORT)"
)

dev.off()

# ---------------------------------------------------------
# 2. PCA visualization - from 02_CIBERSORT_pipeline.R
# ---------------------------------------------------------
pca <- prcomp(malignant_norm, scale. = TRUE)

pca_df <- data.frame(
  pca$x,
  Cell_State = dominant_state
)

pdf("plots/CIBERSORT_AML_PCA.pdf", width = 7, height = 5)

ggplot(pca_df, aes(PC1, PC2, color = Cell_State)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic() +
  labs(title = "PCA of AML Cell States")

dev.off()

message("Plots Saved")

gc()
message("CIBERSORT pipeline complete")
