# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# DESeq2 size factors (median-of-ratios) of the TE consensus peaks, next to
# the library sizes (kept insertions) that CPM divides by. scale_factor scales
# the insertions of a sample to insertions per million of a typical library:
# 1e6 / (size factor x geometric mean of the library sizes). Samples whose
# size factor and library factor are the same get the same signal with both
# scalings; a ratio above 1 means more insertions in the peaks (e.g. a higher
# FRiP), so the DESeq2-scaled signal is lower than the CPM signal

library(DESeq2)

samples <- unlist(snakemake@params[["samples"]])

# Tn5 insertions per peak; the samples are the columns after the annotation.
# Same peak filter as the DESeq2 fit (deseq2_fit.R), so the size factors are
# those of the differential test; they do not depend on the design, so a ~1
# design is used and no fit is needed
counts <- read.delim(
  snakemake@input[["counts"]],
  comment.char = "#",
  row.names = 1,
  check.names = FALSE
)
counts <- as.matrix(counts[, -(1:5), drop = FALSE])
colnames(counts) <- sub("\\.filtered\\.bam$", "", basename(colnames(counts)))
counts <- counts[, samples, drop = FALSE]
counts <- counts[rowSums(counts) >= snakemake@params[["min_count"]], ]

dds <- DESeqDataSetFromMatrix(
  counts,
  data.frame(row.names = samples),
  design = ~1
)
dds <- estimateSizeFactors(dds)
stats <- read.delim(snakemake@input[["stats"]], check.names = FALSE)

size_factors <- sizeFactors(dds)
samples <- names(size_factors)
insertions <- stats$kept_reads[match(samples, stats$sample)]
if (anyNA(insertions)) {
  stop(
    "No kept_reads in te_library_stats.tsv for some samples of the DESeq2 object"
  )
}
geo_mean <- exp(mean(log(insertions)))

result <- data.frame(
  sample = samples,
  insertions = insertions,
  library_factor = insertions / geo_mean,
  size_factor = unname(size_factors),
  ratio = unname(size_factors) / (insertions / geo_mean),
  scale_factor = 1e6 / (unname(size_factors) * geo_mean)
)
print(result)
write.table(
  result,
  snakemake@output[[1]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
