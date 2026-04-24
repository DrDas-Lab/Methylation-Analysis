############################################################
## 05_RNA_heatmap.R
## RNA differentiation trajectory heatmap
############################################################

message("Generating RNA trajectory heatmap")
graphics.off()

# library(pheatmap)
# library(tidyverse)

# ---------------------------------------------------------
# 1. Load expression
# ---------------------------------------------------------
expr_clean <- read.table('output/RNA_Expression.txt',
                               header = TRUE, sep = '\t')


# ---------------------------------------------------------
# 2. Load trajectory ordering
# ---------------------------------------------------------
traj_df <- read.csv("output/AML_Trajectory_Ordering.csv")

common_genes <- as.matrix(expr_clean)
common_samples <- intersect(traj_df$Sample, colnames(common_genes))
traj_df_common <- traj_df[traj_df$Sample %in% common_samples, ]
# length(row.names(traj_df_common))

trajectory_samples <- traj_df_common$Sample
trajectory_states  <- traj_df_common$State
trajectory_scores  <- traj_df_common$Score

# ---------------------------------------------------------
# 3. Define program genes - From fig 6F
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

expr_plot <- expr_clean[genes_use, trajectory_samples]

# ---------------------------------------------------------
# 4. Z-score
# ---------------------------------------------------------
expr_z <- t(scale(t(expr_plot)))
expr_z[!is.finite(expr_z)] <- 0

# ---------------------------------------------------------
# 5. Row annotation
# ---------------------------------------------------------
gene_program <- sapply(rownames(expr_z), function(g){
  names(Filter(function(x) g %in% x, program_list))[1]
})

annotation_row <- data.frame(
  Program = factor(gene_program, levels = names(program_list)),
  row.names = rownames(expr_z)
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
# 7. Colors
# ---------------------------------------------------------
ann_colors <- list(
  
  State = c(
    HSC_like = "#bd3030",
    GMP_like = "#3b3b3b",
    ProMono_like = "#0073c2",
    Mono_like = "#003c67",
    cDC_like = "#efc000"
    ),

  Score = colorRampPalette(c("#e5f5e0", "#31a")) (100)
)

# ---------------------------------------------------------
# 8. Plot
# ---------------------------------------------------------
pdf("plots/RNA_Heatmap.pdf", width = 20, height = 12)


pheatmap(
  mat = expr_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,   # preserves trajectory
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 12,
  fontsize_col = 8,
  # annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  main = "RNA Expression",
  
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cellwidth  = 1.5,
  cellheight = 15,
  border_color = NA
)

dev.off()
gc()
