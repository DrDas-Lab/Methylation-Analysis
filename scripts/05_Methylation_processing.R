############################################################
## 05_Methylation_processing.R
## Promoter methylation processing + heatmap
############################################################

message("Processing methylation data")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
getwd()

library(tidyverse)
library(pheatmap)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ---------------------------------------------------------
# 1. Load methylation beta values
# ---------------------------------------------------------
meth_raw <- read.csv(
  "input/GSE159907_BetaValues_BA_mapped_clean.csv",
  check.names = FALSE
)

colnames(meth_raw)[1] <- "CpG_ID"

# ---------------------------------------------------------
# 2. Annotation
# ---------------------------------------------------------
probe_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
probe_anno <- as.data.frame(probe_anno)
probe_anno$CpG_ID <- rownames(probe_anno)
# print(probe_anno)

meth_full <- merge(probe_anno, meth_raw, by="CpG_ID")

# ---------------------------------------------------------
# 3. Promoter filtering
# ---------------------------------------------------------
meth_prom <- meth_full %>%
  filter(grepl("TSS2000|TSS1500", UCSC_RefGene_Group))

# meth_prom %>% filter(UCSC_RefGene_Name == "ELANE") %>% head()

# write.csv(meth_prom, 'Methylation_TSS_new.csv', row.names = FALSE)
meth_prom <- read.csv('Methylation_TSS_new.csv', check.names = FALSE)

anno_cols <- colnames(meth_prom)[1:46]
meth_samples <- setdiff(colnames(meth_prom), anno_cols)

length(meth_samples)

# ---------------------------------------------------------
# 4. Collapse to gene level
# ---------------------------------------------------------
meth_gene <- meth_prom %>%
  select(UCSC_RefGene_Name, all_of(meth_samples)) %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(across(all_of(meth_samples), mean, na.rm=TRUE))

meth_gene <- as.data.frame(meth_gene)
rownames(meth_gene) <- meth_gene$UCSC_RefGene_Name
meth_gene$UCSC_RefGene_Name <- NULL

# ---------------------------------------------------------
# 5. Load RNA heatmap gene order
# ---------------------------------------------------------
rna_heatmap <- read.table(
  "temp_output/CIBERSORT_mixture_RPM.txt",
  header = TRUE,
  row.names = 1
)

df <- read.table('input/Methylation_GeneLevel_TSS2000_TSS1500.txt',
                 header = TRUE)
list_genes <- list(df[[1]])
print(list_genes)

rna_filtered <- rna_heatmap[rownames(rna_heatmap) %in% unlist(list_genes), ]

common_samples <- intersect(colnames(rna_heatmap), colnames(meth_gene))

common_genes <- intersect(rownames(rna_heatmap), rownames(meth_gene))
meth_gene <- meth_gene[common_genes,common_samples]

expr_clean <- as.matrix(df)
dim(expr_clean)

write.csv(expr_clean, 'Methylation_expression.csv', row.names = TRUE)

traj_df <- read.csv("temp_output/AML_Trajectory_Ordering.csv")

common_samples <- intersect(traj_df$Sample, colnames(expr_clean))
# length(common_samples)
traj_df_common <- traj_df[traj_df$Sample %in% common_samples, ]

trajectory_samples <- traj_df_common$Sample
trajectory_states  <- traj_df_common$State
trajectory_scores  <- traj_df_common$Score

# ---------------------------------------------------------
# 3. Define program genes
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
  
  Score = colorRampPalette(c("#e5f5e0", "#31a")) (100)
)

pdf("temp_output/Methylation_Trajectory_Heatmap.pdf", width = 20, height = 12)

pheatmap(
  mat = meth_z,
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
  main = "Promoter DNA Methylation",
  
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cellwidth  = 1.2,
  cellheight = 15,
  border_color = NA
)


dev.off()

