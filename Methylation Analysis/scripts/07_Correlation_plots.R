############################################################
## 07_Correlation_plots.R
## RNA vs methylation correlation plots
############################################################

message("Generating correlation plots")

graphics.off()

# library(tidyverse)
# library(ggplot2)
# library(ggridges)

# cor_df <- read.csv("temp_output/Gene_RNA_Methylation_Correlations.csv")
expr_z <- read.table("output/RNA_Expression.txt",
                     sep = '\t',
                     check.names = T, 
                     row.names = 1)
dim(expr_z)
print(colnames(expr_z)) 
print(rownames(expr_z))

meth_z <- read.table("output/Meth_Expression.txt", 
                     sep = '\t',
                     check.names = T, 
                     row.names = 1)

dim(meth_z)
print(colnames(meth_z)) 
print(rownames(meth_z))

expr_z_sorted <- expr_z[rownames(meth_z), colnames(meth_z)]
common_genes <- intersect(rownames(expr_z), rownames(meth_z))

program_list <- list(
  "mature myeloid" = c("MPO","ELANE","AZU1","CTSG","SELL","CD36","SRGN","MNDA","CX3CR1","FCGR3A","MAFB","CD14","FCN1"),
  "immune regulation" = c("CD74","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","TYROBP","LILRA5","LILRB2","FCER1G","CTSS"),
  "proinflammatory signaling" = c("MAP3K8","CFD","CXCL8","NCF2","S100A12","SERPINA1","CCL5","S100A8","S100A9"),
  "AML core" = c("CDK6","MYB","HOXA9"),
  "stress / metabolism" = c("MGST1","CLU","FOS","HSPA5","CYP1B1")
)

genes_all <- unique(unlist(program_list))
genes_use <- intersect(genes_all, rownames(expr_z_sorted))

# expr_plot <- expr_clean[genes_use, trajectory_samples]

# ---------------------------------------------------------
# 4. Z-score
# # ---------------------------------------------------------
# expr_z <- t(scale(t(expr_plot)))
# expr_z[!is.finite(expr_z)] <- 0

# ---------------------------------------------------------
# 5. Row annotation
# ---------------------------------------------------------
gene_program <- sapply(rownames(expr_z_sorted), function(g){
  names(Filter(function(x) g %in% x, program_list))[1]
})

annotation_row <- data.frame(
  Program = factor(gene_program, levels = names(program_list)),
  row.names = rownames(expr_z)
)

cor_list <- lapply(common_genes, function(g) {
  
  rna  <- as.numeric(expr_z_sorted[g, ])
  meth <- as.numeric(meth_z[g, ])
  
  ok <- is.finite(rna) & is.finite(meth)
  if (sum(ok) < 5) return(NULL)
  
  prog <- annotation_row[g, "Program", drop = TRUE]
  if (is.na(prog)) return(NULL)
  
  ct <- cor.test(rna[ok], meth[ok], method = "spearman", exact = FALSE)
  
  data.frame(
    Gene        = g,
    Correlation = unname(ct$estimate),
    p_value     = ct$p.value,
    Program     = as.character(prog),
    stringsAsFactors = FALSE
  )
})

cor_df <- bind_rows(cor_list)

stopifnot(nrow(cor_df) > 0)

cor_df <- cor_df %>%
  mutate(
    logP = -log10(p_value),
    logP = pmin(logP, 20),
    Program = factor(
      Program,
      levels = c(
        "mature myeloid",
        "immune regulation",
        "proinflammatory signaling",
        "AML core",
        "stress / metabolism"
      )
    )
  )


custom_colors <- c(
  "mature myeloid" = "#fc8d62",
  "immune regulation" = "#8da0cb",
  "proinflammatory signaling" = "#e78ac3",
  "AML core" = "#66c2a5",
  "stress / metabolism" = "#ffd92f"
)

# ---------------------------------------------------------
# Ridge plot
# ---------------------------------------------------------
# cor_df$Program <- factor(
#   cor_df$Program,
#   levels = c(
#     "mature myeloid",
#     "immune regulation",
#     "proinflammatory signaling",
#     "AML core",
#     "stress / metabolism"
#   )
# )

pdf("plots/Correlation_RidgePlot.pdf", width=8, height=6)

ggplot(cor_df, aes(Correlation, Program, fill=Program)) +
  geom_density_ridges(alpha=.8, color="white") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_y_discrete(limits = rev(levels(cor_df$Program))) +
  scale_fill_manual(values = custom_colors) +
  # coord_cartesian(xlim = c(-0.50, 0.25)) +
  theme_classic()

dev.off()

# ---------------------------------------------------------
# Mirror plot
# ---------------------------------------------------------
cor_df <- cor_df %>%
  mutate(logP = -log10(p_value))

pdf("plots/Correlation_MirrorPlot.pdf", width=10, height=15)

ggplot(cor_df,
       aes(Correlation, reorder(Gene, Correlation), fill=Program)) +
  geom_col(alpha=.9, width=0.5) +
  geom_point(aes(size=logP), color="black") +
  geom_vline(xintercept=0) +
  facet_grid(Program~., scales="free_y", space="free") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

dev.off()
gc()