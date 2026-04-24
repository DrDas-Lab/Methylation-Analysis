############################################################
## 01_preprocess_counts.R
############################################################
# library(data.table)
# library(org.Hs.eg.db)

output_file <- "output/BAML_Normalized.txt"
req_file <- "input/BAML_Normalized.txt"

if (file.exists(req_file)) {
  
  message(sprintf("File '%s' already exists. Skipping preprocessing.", output_file))
  
} else {
  
  message("Running count preprocessing...")
  

  counts <- fread("input/beataml_waves1to4_counts_dbgap.txt", data.table = FALSE)
  gene_ids <- counts[, 1]
  count_matrix <- counts[, -1]
  counts_matrix <- as.matrix(counts[, sapply(counts, is.numeric)])
  rownames(count_matrix) <- gene_ids
  
  lib_size <- colSums(counts_matrix)
  rpm <- t(t(counts_matrix) / lib_size) * 1e6
  log2rpm <- log2(rpm + 1)
  

  write.table(
    log2rpm,
    output_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  message(sprintf("Preprocessing complete. Processed %d genes.", nrow(log2rpm)))
}