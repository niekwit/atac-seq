# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Differential accessibility of repeat subfamilies between two conditions (TE
# module): DESeq2 Wald test, with log2 fold changes shrunken by ashr for
# ranking and plots

library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load Snakemake parameters
treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
fdr <- snakemake@params[["fdr"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
n_labels <- 15

dds <- readRDS(snakemake@input[["dds"]])
contrast <- c("condition", treatment, control)

# The test (and thus the p-values) uses the unshrunken fold changes
res <- results(
  dds,
  contrast = contrast,
  alpha = fdr,
  lfcThreshold = lfc_threshold
)
summary(res)
shrunken <- lfcShrink(dds, contrast = contrast, res = res, type = "ashr")

table <- data.frame(
  subfamily = rownames(res),
  baseMean = res$baseMean,
  log2FoldChange = res$log2FoldChange,
  log2FoldChange_shrunken = shrunken$log2FoldChange,
  lfcSE = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj
)
table$change <- "unchanged"
significant <- !is.na(table$padj) & table$padj < fdr
table$change[significant & table$log2FoldChange > 0] <- "more accessible"
table$change[significant & table$log2FoldChange < 0] <- "less accessible"
print(table(table$change))

# Add family, class, number of copies and genome coverage of each subfamily
subfamilies <- read.delim(
  snakemake@input[["subfamilies"]],
  check.names = FALSE
)
output <- merge(subfamilies, table, by = "subfamily")
output <- output[order(output$padj, output$pvalue, na.last = TRUE), ]
write.table(
  output,
  snakemake@output[["table"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

colours <- c(
  "more accessible" = "#D55E00",
  "less accessible" = "#0072B2",
  "unchanged" = "grey70"
)
title <- paste(treatment, "vs", control)
# Label the most significant subfamilies
labelled <- head(output$subfamily[output$change != "unchanged"], n_labels)
table$label <- ifelse(table$subfamily %in% labelled, table$subfamily, "")
table$change <- factor(table$change, levels = rev(names(colours)))
table <- table[order(table$change), ]

# MA plot
ma <- ggplot(
  table,
  aes(baseMean, log2FoldChange_shrunken, colour = change, label = label)
) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = 0) +
  geom_text_repel(
    size = 2.5,
    colour = "black",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_x_log10() +
  scale_colour_manual(values = colours, drop = FALSE) +
  labs(
    title = title,
    x = "Mean normalised insertions",
    y = "Shrunken log2 fold change",
    colour = NULL
  ) +
  theme_bw(base_size = 12)
ggsave(snakemake@output[["ma"]], ma, width = 6, height = 4.5, dpi = 300)

# Volcano plot
volcano <- ggplot(
  table[!is.na(table$pvalue), ],
  aes(log2FoldChange_shrunken, -log10(pvalue), colour = change, label = label)
) +
  geom_point(size = 0.8) +
  geom_text_repel(
    size = 2.5,
    colour = "black",
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_colour_manual(values = colours, drop = FALSE) +
  labs(
    title = title,
    x = "Shrunken log2 fold change",
    y = "-log10(p-value)",
    colour = NULL
  ) +
  theme_bw(base_size = 12)
ggsave(
  snakemake@output[["volcano"]],
  volcano,
  width = 6,
  height = 4.5,
  dpi = 300
)
