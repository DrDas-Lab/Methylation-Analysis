############################################################
## 03_Processing_RNA_Methylation_datasets.R
############################################################

library(data.table)
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ---------------------------------------------------------
# 1. Load RNA CIBERSORT Mixture file
# ---------------------------------------------------------
rna <- read.table("output/CIBERSORT_mixture_RPM.txt",
                  header = TRUE,
                  sep = "\t",
                  row.names = 1,
                  check.names = FALSE)

# ---------------------------------------------------------
# 2. Load methylation beta values
# ---------------------------------------------------------
meth_raw <- read.csv("input/MethylationBetaValues_clean.csv",
                     check.names = FALSE)

colnames(meth_raw)[1] <- "CpG_ID"

# ---------------------------------------------------------
# 3. Annotation of Methylation data
# ---------------------------------------------------------
probe_anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
probe_anno <- as.data.frame(probe_anno)
probe_anno$CpG_ID <- rownames(probe_anno)
# print(probe_anno)

meth_full <- merge(probe_anno, meth_raw, by="CpG_ID")

# ---------------------------------------------------------
# 4. Promoter filtering
# ---------------------------------------------------------
meth_prom <- meth_full %>%
  filter(grepl("TSS2000|TSS1500", UCSC_RefGene_Group))

# write.csv(meth_prom, 'Methylation_TSS_new.csv', row.names = FALSE)
# meth_prom <- read.csv('Methylation_TSS_new.csv', check.names = FALSE)

anno_cols <- colnames(meth_prom)[1:46]
meth_samples <- setdiff(colnames(meth_prom), anno_cols)

# ---------------------------------------------------------
# 4. Collapse to gene level
# ---------------------------------------------------------
meth_gene <- meth_prom %>%
  dplyr::select(UCSC_RefGene_Name, all_of(meth_samples)) %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(across(all_of(meth_samples), mean, na.rm=TRUE))

meth_gene <- as.data.frame(meth_gene)
rownames(meth_gene) <- meth_gene$UCSC_RefGene_Name
meth_gene$UCSC_RefGene_Name <- NULL

# ---------------------------------------------------------
# 4. RNA and Methylation filtering - Based on Figure 6D
# ---------------------------------------------------------
program_list <- list(
  "mature myeloid" = c("MPO","ELANE","AZU1","CTSG","SELL","CD36","SRGN","MNDA","CX3CR1","FCGR3A","MAFB","CD14","FCN1"),
  "immune regulation" = c("CD74","HLA-DPA1","HLA-DRB1","HLA-DPB1","HLA-DQA1","TYROBP","LILRA5","LILRB2","FCER1G","CTSS"),
  "proinflammatory signaling" = c("MAP3K8","CFD","CXCL8","NCF2","S100A12","SERPINA1","CCL5","S100A8","S100A9"),
  "AML core" = c("CDK6","MYB","HOXA9"),
  "stress / metabolism" = c("MGST1","CLU","FOS","HSPA5","CYP1B1")
)

valid_genes <- Reduce(intersect, list(unlist(program_list), rownames(rna), rownames(meth_gene)))

common_samples <- intersect(colnames(rna), colnames(meth_gene))

rna_targeted <- rna[valid_genes, common_samples, drop = FALSE]
meth_targeted <- meth_gene[valid_genes, common_samples, drop = FALSE]

# ---------------------------------------------------------
# 4. Export the datasets
# ---------------------------------------------------------
write.table(rna_targeted, 
            file = "output/RNA_Expression.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

write.table(meth_targeted, 
            file = "output/Meth_Expression.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE)

message(sprintf("RNA Dims: %d x %d | Meth Dims: %d x %d", 
                nrow(rna_targeted), ncol(rna_targeted), 
                nrow(meth_targeted), ncol(meth_targeted)))

gc()