# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Fits the DESeq2 model to the Tn5 insertions in the consensus peaks of all
# samples, so that all replicates contribute to the dispersion estimates

library(BiocParallel)
library(DESeq2)

# Load Snakemake parameters
samples <- unlist(snakemake@params[["samples"]])
conditions <- unlist(snakemake@params[["conditions"]])
# Reference condition first
levels <- unlist(snakemake@params[["levels"]])

# Tn5 insertions per peak; the samples are the columns after the annotation
counts <- read.delim(
  snakemake@input[["counts"]],
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)
counts <- as.matrix(counts[, -(1:5), drop = FALSE])
colnames(counts) <- sub(".nodup.bam$", "", basename(colnames(counts)))
counts <- counts[, samples, drop = FALSE]

sample_info <- data.frame(
  condition = factor(conditions, levels = levels),
  row.names = samples
)
dds <- DESeqDataSetFromMatrix(counts, sample_info, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(
  dds,
  parallel = TRUE,
  BPPARAM = MulticoreParam(snakemake@threads)
)
saveRDS(dds, snakemake@output[["dds"]])
