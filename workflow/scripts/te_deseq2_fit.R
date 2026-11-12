# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Fits the DESeq2 model (~ condition) to the Tn5 insertions per repeat
# subfamily of all samples (TE module), so that all replicates contribute to
# the dispersion estimates

library(BiocParallel)
library(DESeq2)

# Load Snakemake parameters
samples <- unlist(snakemake@params[["samples"]])
conditions <- unlist(snakemake@params[["conditions"]])
# Reference condition first
levels <- unlist(snakemake@params[["levels"]])
normalisation <- snakemake@params[["normalisation"]]

# Insertions per subfamily (all read pairs, multimappers at one of their best
# alignments); the samples are the columns after the annotation
counts <- read.delim(
  snakemake@input[["counts"]],
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)
counts <- as.matrix(counts[, -(1:5), drop = FALSE])
colnames(counts) <- sub(".filtered.bam$", "", basename(colnames(counts)))
counts <- counts[, samples, drop = FALSE]

sample_info <- data.frame(
  condition = factor(conditions, levels = levels),
  row.names = samples
)
dds <- DESeqDataSetFromMatrix(counts, sample_info, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

if (normalisation == "insertions") {
  # Library size: all kept insertions of a sample, in and outside repeats
  stats <- read.delim(snakemake@input[["stats"]], row.names = 1)
  insertions <- stats[samples, "kept_reads"]
  sizeFactors(dds) <- insertions / exp(mean(log(insertions)))
} else {
  # Median-of-ratios over the subfamilies (DESeq2 default)
  dds <- estimateSizeFactors(dds)
}
print(sizeFactors(dds))

dds <- DESeq(
  dds,
  parallel = TRUE,
  BPPARAM = MulticoreParam(snakemake@threads)
)
saveRDS(dds, snakemake@output[["dds"]])
