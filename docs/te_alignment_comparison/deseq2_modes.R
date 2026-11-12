# Differential accessibility of TE subfamilies (treatment vs reference
# condition) with DESeq2 for each bowtie2 alignment mode, to check that the mode
# does not change the result
#
# Usage: Rscript deseq2_modes.R <analysis directory>
# Writes te_alignment_comparison/comparison/deseq2_*.tsv

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
c_dir <- file.path(args[1], "te_alignment_comparison")
out_dir <- file.path(c_dir, "comparison")
dir.create(out_dir, showWarnings = FALSE)

modes <- c("default", "very_sensitive", "k10", "k20", "k100")
modes <- modes[file.exists(file.path(
  c_dir,
  "counts",
  paste0(modes, ".multi.txt")
))]
young <- c(
  "L1MdTf_I",
  "L1MdTf_II",
  "L1MdTf_III",
  "L1MdA_I",
  "L1MdA_II",
  "L1MdA_III",
  "L1MdGf_I",
  "L1MdGf_II",
  "L1MdN_I",
  "L1MdF_I",
  "L1MdF_II",
  "IAPEz-int",
  "IAPLTR1_Mm",
  "IAPLTR1a_Mm"
)

sample_sheet <- read.csv(file.path(args[1], "config", "samples.csv"))
reference <- unique(sample_sheet$condition[sample_sheet$reference == "yes"])
treatment <- setdiff(unique(sample_sheet$condition), reference)
stopifnot(length(reference) == 1, length(treatment) == 1)
sample_sheet$condition <- factor(
  sample_sheet$condition,
  levels = c(reference, treatment)
)
rownames(sample_sheet) <- sample_sheet$sample

results_all <- list()
for (mode in modes) {
  counts <- read.delim(
    file.path(c_dir, "counts", paste0(mode, ".multi.txt")),
    comment.char = "#",
    check.names = FALSE
  )
  rownames(counts) <- counts$Geneid
  counts <- counts[, -(1:6)]
  colnames(counts) <- sub("\\..*", "", basename(colnames(counts)))
  # featureCounts -M --fraction gives fractional counts
  counts <- round(as.matrix(counts[, sample_sheet$sample]))
  counts <- counts[rowSums(counts) >= 10, ]

  dds <- DESeqDataSetFromMatrix(counts, sample_sheet, design = ~condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", treatment, reference))
  res <- data.frame(mode = mode, subfamily = rownames(res), as.data.frame(res))
  results_all[[mode]] <- res
}
res <- do.call(rbind, results_all)
rownames(res) <- NULL
write.table(
  res,
  file.path(out_dir, "deseq2_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

res_young <- res[res$subfamily %in% young, ]
write.table(
  res_young,
  file.path(out_dir, "deseq2_young.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("## log2 fold change %s vs %s\n", treatment, reference))
print(round(
  xtabs(log2FoldChange ~ subfamily + mode, res_young)[young, modes],
  2
))
cat("\n## Adjusted p-value\n")
print(signif(xtabs(padj ~ subfamily + mode, res_young)[young, modes], 2))

# Agreement of log2 fold changes and significant calls with -k 100
if ("k100" %in% modes) {
  ref <- results_all[["k100"]]
  cat("\n## Agreement with k100 (subfamilies with baseMean >= 50)\n")
  for (mode in setdiff(modes, "k100")) {
    x <- results_all[[mode]]
    shared <- intersect(
      rownames(x)[x$baseMean >= 50],
      rownames(ref)[ref$baseMean >= 50]
    )
    sig_x <- shared[which(x[shared, "padj"] < 0.05)]
    sig_ref <- shared[which(ref[shared, "padj"] < 0.05)]
    cat(sprintf(
      "%s: n = %d, Pearson LFC = %.4f, padj < 0.05: %d (k100 %d, shared %d)\n",
      mode,
      length(shared),
      cor(x[shared, "log2FoldChange"], ref[shared, "log2FoldChange"]),
      length(sig_x),
      length(sig_ref),
      length(intersect(sig_x, sig_ref))
    ))
  }
}
