############################################################
## 06_Methylation_processing.R
## Promoter methylation processing + heatmap
############################################################
message("Processing methylation data")

# library(tidyverse)
# library(pheatmap)

graphics.off()

# ---------------------------------------------------------
# 1. Load methylation Expression Data
# ---------------------------------------------------------
df <- read.table('output/Meth_Expression.txt',
                 header = TRUE, sep = '\t')

expr_clean <- as.matrix(df)
dim(expr_clean)

# write.csv(expr_clean, 'Methylation_expression.csv', row.names = TRUE)

traj_df <- read.csv("output/AML_Trajectory_Ordering.csv")

common_samples <- intersect(traj_df$Sample, colnames(expr_clean))
length(common_samples)
traj_df_common <- traj_df[traj_df$Sample %in% common_samples, ]
length(rownames(traj_df_common))

trajectory_samples <- traj_df_common$Sample
trajectory_states  <- traj_df_common$State
trajectory_scores  <- traj_df_common$Score

# ---------------------------------------------------------
# 2. Define program genes
# ---------------------------------------------------------
program_list <- list(
  "mature myeloid" = c("MPO","ELANE","AZU1","CTSG","SELL","CD36","SRGN","MNDA","CX3CR1","FCGR3A","MAFB","CD14","FCN1"),
  "immune regulation" = c("CD74","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","TYROBP","LILRA5","LILRB2","FCER1G","CTSS"),
  "proinflammatory signaling" = c("MAP3K8","CFD","CXCL8","NCF2","S100A12","SERPINA1","CCL5","S100A8","S100A9"),
  "AML core" = c("CDK6","MYB","HOXA9"),
  "stress / metabolism" = c("MGST1","CLU","FOS","HSPA5","CYP1B1")
)

genes_all <- unique(unlist(program_list))
genes_use <- intersect(genes_all, rownames(expr_clean))
samples_use <- intersect(trajectory_samples, colnames(expr_clean))

expr_plot <- expr_clean[genes_use, samples_use]
# ---------------------------------------------------------
# 4. Z-score
# ---------------------------------------------------------

meth_z <- t(scale(t(expr_plot)))
meth_z[!is.finite(meth_z)] <- 0

# ---------------------------------------------------------
# 5. Row annotation
# ---------------------------------------------------------
gene_program <- sapply(rownames(meth_z), function(g){
  names(Filter(function(x) g %in% x, program_list))[1]
})

annotation_row <- data.frame(
  Program = factor(gene_program, levels = names(program_list)),
  row.names = rownames(meth_z)
)

gaps_row <- cumsum(table(annotation_row$Program))

# ---------------------------------------------------------
# 6. Column annotation
# ---------------------------------------------------------
annotation_col <- data.frame(
  State = factor(trajectory_states, levels = unique(trajectory_states)),
  Score = trajectory_scores,
  row.names = trajectory_samples
)

annotation_col$Score[annotation_col$Score > 1] <- 1
annotation_col$Score[annotation_col$Score < 0.4] <- 0.4

gaps_col <- cumsum(table(annotation_col$State))

# ---------------------------------------------------------
# 7. Plot
# ---------------------------------------------------------
ann_colors <- list(
  
  State = c(
    HSC_like = "#bd3030",
    GMP_like = "#3b3b3b",
    ProMono_like = "#0073c2",
    Mono_like = "#003c67",
    cDC_like = "#efc000"
  ),
  
  Score = colorRampPalette(c("#e5f5e0", "#3311aa")) (100)
)

pdf("plots/Methylation_Heatmap.pdf", width = 20, height = 12)

pheatmap(
  mat = meth_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,   # preserves trajectory
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 12,
  fontsize_col = 8,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  main = "Promoter DNA Methylation",
  
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cellwidth  = 1.5,
  cellheight = 15,
  border_color = NA
)


dev.off()
gc()
