# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# PCA and correlation of all samples, from the Tn5 insertions in the consensus
# peaks after the variance stabilising transformation of DESeq2

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Load Snakemake parameters
samples <- unlist(snakemake@params[["samples"]])
conditions <- unlist(snakemake@params[["conditions"]])
# Reference condition first
levels <- unlist(snakemake@params[["levels"]])
ntop <- snakemake@params[["ntop"]]

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

# design ~ 1: the transformation does not use the conditions (blind)
sample_info <- data.frame(
  condition = factor(conditions, levels = levels),
  row.names = samples
)
dds <- DESeqDataSetFromMatrix(counts, sample_info, design = ~1)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- estimateSizeFactors(dds)
print(sizeFactors(dds))

# vst() needs at least 1000 peaks with a mean count above 5
if (sum(rowMeans(counts(dds, normalized = TRUE)) > 5) >= 1000) {
  vsd <- vst(dds)
} else {
  vsd <- varianceStabilizingTransformation(dds)
}
vst_counts <- assay(vsd)
write.table(
  data.frame(peak = rownames(vst_counts), vst_counts, check.names = FALSE),
  snakemake@output[["vst"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# PCA of the most variable peaks
pca <- plotPCA(vsd, intgroup = "condition", ntop = ntop, returnData = TRUE)
percent_var <- round(100 * attr(pca, "percentVar"))
write.table(
  pca[, c("name", "condition", "PC1", "PC2")],
  snakemake@output[["pca_data"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
plot <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE, size = 3) +
  labs(
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 12)
ggsave(snakemake@output[["pca"]], plot, width = 6, height = 4.5, dpi = 300)

# Pearson correlation of the samples over all peaks
correlation <- cor(vst_counts)
write.table(
  data.frame(sample = rownames(correlation), correlation, check.names = FALSE),
  snakemake@output[["correlation_data"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
size <- 2.5 + 0.5 * length(samples)
pheatmap(
  correlation,
  annotation_col = sample_info,
  display_numbers = TRUE,
  number_format = "%.3f",
  filename = snakemake@output[["correlation"]],
  width = size + 1.5,
  height = size
)
