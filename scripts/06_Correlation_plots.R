############################################################
## 06_Correlation_plots.R
## RNA vs methylation correlation plots
############################################################

message("Generating correlation plots")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
getwd()

library(tidyverse)
library(ggplot2)
library(ggridges)

# cor_df <- read.csv("temp_output/Gene_RNA_Methylation_Correlations.csv")
expr_z <- read.csv("RNA_expression.csv", check.names = FALSE, row.names = 1)
meth_z <- read.csv("Methylation_expression.csv", check.names = FALSE, row.names = 1)
expr_z_sorted <- expr_z[rownames(meth_z), colnames(meth_z)]

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
print(cor_list)

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
cor_df$Program <- factor(
  cor_df$Program,
  levels = c(
    "mature myeloid",
    "immune regulation",
    "proinflammatory signaling",
    "AML core",
    "stress / metabolism"
  )
)

pdf("temp_output/Correlation_RidgePlot.pdf", width=8, height=6)

ggplot(cor_df, aes(Correlation, Program, fill=Program)) +
  geom_density_ridges(alpha=.8, color="white") +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_y_discrete(limits = rev(levels(cor_df$Program))) +
  scale_fill_manual(values = custom_colors) +
  coord_cartesian(xlim = c(-0.50, 0.25)) +
  theme_classic()

dev.off()

# ---------------------------------------------------------
# Mirror plot
# ---------------------------------------------------------
cor_df <- cor_df %>%
  mutate(logP = -log10(p_value))

pdf("temp_output/Correlation_MirrorPlot.pdf", width=10, height=15)

ggplot(cor_df,
       aes(Correlation, reorder(Gene, Correlation), fill=Program)) +
  geom_col(alpha=.9, width=0.5) +
  geom_point(aes(size=logP), color="black") +
  geom_vline(xintercept=0) +
  facet_grid(Program~., scales="free_y", space="free") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal()

dev.off()
