############################################################
## 04_ModuleScore_trajectory.R
## Purpose:
## Calculate AML differentiation module scores
## and generate ordered trajectory samples
############################################################

message("Running ModuleScore trajectory pipeline")

# ---------------------------------------------------------
# 1. Load libraries
# ---------------------------------------------------------
# library(Seurat)
# library(tidyverse)

# ---------------------------------------------------------
# 2. Load cleaned RPM expression (mixture-ready)
# ---------------------------------------------------------
expr_clean <- read.table(
  "output/CIBERSORT_mixture_RPM.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

expr_clean <- as.matrix(expr_clean)

# ---------------------------------------------------------
# 3. Load dominant AML cell states (from CIBERSORT step)
# ---------------------------------------------------------
dominant_df <- read.csv("output/AML_dominant_cell_state.csv")

rownames(dominant_df) <- dominant_df$Sample
dominant_states <- dominant_df$Dominant_Cell_State

# ---------------------------------------------------------
# 4. Load AML cell-state gene sets (GMT)
# ---------------------------------------------------------
gmt_file <- "input/AMLCellType_Genesets.gmt"

gmt_raw <- readLines(gmt_file)
gmt_list <- strsplit(gmt_raw, "\t")

genesets <- lapply(gmt_list, function(x) x[-c(1,2)])
names(genesets) <- sapply(gmt_list, `[`, 1)

# ---------------------------------------------------------
# 5. Merge gene sets by biological state
# ---------------------------------------------------------
merge_sets <- function(pattern) {
  unique(unlist(genesets[grep(pattern, names(genesets))]))
}

final_genesets <- list(
  HSC          = merge_sets("^LSPC"),
  GMP_like     = merge_sets("^GMP-like"),
  ProMono_like = merge_sets("^ProMono-like"),
  Mono_like    = merge_sets("^Mono-like"),
  cDC_like     = merge_sets("^cDC-like")
)

# ---------------------------------------------------------
# 6. Intersect gene sets with bulk genes
# ---------------------------------------------------------
bulk_genes <- rownames(expr_clean)

final_genesets <- lapply(final_genesets, function(g) {
  intersect(g, bulk_genes)
})

gene_counts <- sapply(final_genesets, length)
print(gene_counts)

if (any(gene_counts < 15)) {
  stop("Some gene sets contain <15 genes")
}

# ---------------------------------------------------------
# 7. Create Seurat object
# ---------------------------------------------------------
aml_seu <- CreateSeuratObject(counts = expr_clean)
aml_seu <- NormalizeData(aml_seu)

# ---------------------------------------------------------
# 8. Add Module Scores
# ---------------------------------------------------------
aml_seu <- AddModuleScore(
  aml_seu,
  features = final_genesets,
  name = "CellState"
)

score_cols <- grep("^CellState", colnames(aml_seu@meta.data), value = TRUE)

module_scores <- aml_seu@meta.data[, score_cols, drop = FALSE]
colnames(module_scores) <- names(final_genesets)

# Save module scores
write.csv(
  module_scores,
  "output/AML_ModuleScores.csv"
)
# 
# head(rownames(module_scores))
# head(names(dominant_states))

# ---------------------------------------------------------
# 9. Build PCA dataframe
# ---------------------------------------------------------
pca <- prcomp(module_scores, scale. = TRUE)

pca_df <- data.frame(
  pca$x,
  Cell_State = dominant_df$Dominant_Cell_State[match(rownames(module_scores), 
                                                     dominant_df$Sample)]
)
print(pca_df)

write.csv(
  pca_df,
  "output/AML_ModuleScore_PCA.csv"
)

# ---------------------------------------------------------
# 10. Define trajectory states
# ---------------------------------------------------------
trajectory_states <- list(
  
  HSC_like = list(
    states = c("LSPC-Quiescent","LSPC-Primed","LSPC-Cycle"),
    score  = "HSC",
    order  = "decreasing"
  ),
  
  GMP_like = list(
    states = c("GMP-like"),
    score  = "GMP_like",
    order  = "decreasing"
  ),
  
  ProMono_like = list(
    states = c("ProMono-like"),
    score  = "ProMono_like",
    order  = "increasing"
  ),
  
  Mono_like = list(
    states = c("Mono-like"),
    score  = "Mono_like",
    order  = "increasing"
  ),
  
  cDC_like = list(
    states = c("cDC-like"),
    score  = "cDC_like",
    order  = "increasing"
  )
)

# ---------------------------------------------------------
# 11. Build ordered trajectory samples
# ---------------------------------------------------------
trajectory_samples <- character()
trajectory_state_labels <- character()
trajectory_scores <- numeric()

# print(trajectory_scores)
for (nm in names(trajectory_states)) {
  
  st <- trajectory_states[[nm]]
  
  samples <- rownames(pca_df)[
    pca_df$Cell_State %in% st$states
  ]
  
  if (length(samples) == 0) next
  
  scores <- module_scores[samples, st$score]
  names(scores) <- samples
  scores <- scores[!is.na(scores)]
  
  if (st$order == "decreasing") {
    scores <- sort(scores, decreasing = TRUE)
  } else {
    scores <- sort(scores, decreasing = FALSE)
  }
  
  trajectory_samples <- c(trajectory_samples, names(scores))
  trajectory_state_labels <- c(
    trajectory_state_labels,
    rep(nm, length(scores))
  )
  trajectory_scores <- c(trajectory_scores, scores)
}

trajectory_df <- data.frame(
  Sample = trajectory_samples,
  State  = trajectory_state_labels,
  Score  = trajectory_scores
)

write.csv(
  trajectory_df,
  "output/AML_Trajectory_Ordering.csv",
  row.names = FALSE
)

message("ModuleScore trajectory pipeline complete")
gc()