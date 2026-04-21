library(data.table)
message("Running count preprocessing")

counts <- fread("data/raw/raw_counts.txt", data.table = FALSE)

rownames(counts) <- counts[,1]
counts <- counts[,-1]

counts <- as.matrix(counts)
mode(counts) <- "numeric"

lib_size <- colSums(counts)

rpm <- t(t(counts) / lib_size) * 1e6
log2rpm <- log2(rpm + 1)

write.table(
  log2rpm,
  "data/processed/log2RPM_expression.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)