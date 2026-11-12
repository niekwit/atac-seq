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

# Heatmap faceted by class: row z-scores of variance-stabilised insertions in
# the samples of the two conditions. Classes and subfamilies are lists from the
# config, where "all" selects everything. With all subfamilies, the significant
# ones are plotted (at most max_per_class per class); named subfamilies are
# plotted whether significant or not, with significant ones marked by *
heatmap_classes <- unlist(snakemake@params[["heatmap"]][["classes"]])
heatmap_subfamilies <- unlist(snakemake@params[["heatmap"]][["subfamilies"]])
max_per_class <- snakemake@params[["heatmap"]][["max_per_class"]]
all_classes <- "all" %in% heatmap_classes
all_subfamilies <- "all" %in% heatmap_subfamilies

if (all_subfamilies) {
  selected <- output[output$change != "unchanged", ]
} else {
  missing <- setdiff(heatmap_subfamilies, output$subfamily)
  if (length(missing) > 0) {
    message(
      "Heatmap: subfamilies not tested (unknown or too few insertions): ",
      paste(missing, collapse = ", ")
    )
  }
  selected <- output[output$subfamily %in% heatmap_subfamilies, ]
}
if (!all_classes) {
  outside <- selected$subfamily[!selected$class %in% heatmap_classes]
  if (!all_subfamilies && length(outside) > 0) {
    message(
      "Heatmap: subfamilies outside the selected classes: ",
      paste(outside, collapse = ", ")
    )
  }
  selected <- selected[selected$class %in% heatmap_classes, ]
}
if (all_subfamilies && nrow(selected) > 0) {
  # output is sorted by significance
  selected <- do.call(
    rbind,
    lapply(split(selected, selected$class), head, max_per_class)
  )
}

if (nrow(selected) == 0) {
  heatmap <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = if (all_subfamilies) {
        "No significant subfamilies in the selected classes"
      } else {
        "None of the selected subfamilies were tested"
      }
    ) +
    labs(title = title) +
    theme_void(base_size = 12)
  width <- 5
  height <- 3
} else {
  in_contrast <- dds$condition %in% c(control, treatment)
  vsd <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  vsd <- vsd[selected$subfamily, in_contrast, drop = FALSE]
  z <- t(scale(t(vsd)))
  # Classes in config order (or by number of rows with all classes); rows
  # within a class by log2 fold change, most increased at the top
  class_levels <- if (all_classes) {
    names(sort(table(selected$class), decreasing = TRUE))
  } else {
    intersect(heatmap_classes, selected$class)
  }
  selected <- selected[
    order(match(selected$class, class_levels), selected$log2FoldChange),
  ]
  selected$label <- selected$subfamily
  if (!all_subfamilies) {
    significant <- selected$change != "unchanged"
    selected$label[significant] <- paste(selected$label[significant], "*")
  }
  long <- data.frame(
    subfamily = rep(rownames(z), times = ncol(z)),
    sample = rep(colnames(z), each = nrow(z)),
    z = as.vector(z)
  )
  rows <- match(long$subfamily, selected$subfamily)
  long$class <- factor(selected$class[rows], levels = class_levels)
  long$label <- factor(selected$label[rows], levels = selected$label)
  conditions <- as.character(dds$condition[in_contrast])
  long$condition <- factor(
    conditions[match(long$sample, colnames(z))],
    levels = c(control, treatment)
  )
  limit <- max(abs(long$z), na.rm = TRUE)
  heatmap <- ggplot(long, aes(sample, label, fill = z)) +
    geom_tile() +
    facet_grid(class ~ condition, scales = "free", space = "free") +
    scale_fill_gradient2(
      low = "#0072B2",
      mid = "white",
      high = "#D55E00",
      limits = c(-limit, limit)
    ) +
    labs(
      title = title,
      caption = if (!all_subfamilies) "* adjusted p < FDR" else NULL,
      x = NULL,
      y = NULL,
      fill = "z-score"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0),
      panel.grid = element_blank()
    )

  # Square cells: give each panel an absolute size of its number of samples
  # and subfamilies times the cell size, then size the figure to fit
  cell <- grid::unit(0.2, "in")
  heatmap <- ggplotGrob(heatmap)
  panels <- heatmap$layout[grepl("^panel", heatmap$layout$name), ]
  panel_columns <- sort(unique(panels$l))
  panel_rows <- sort(unique(panels$t))
  samples_per_condition <- table(long$condition[!duplicated(long$sample)])
  subfamilies_per_class <- table(selected$class)[class_levels]
  for (i in seq_along(panel_columns)) {
    heatmap$widths[panel_columns[i]] <- samples_per_condition[[i]] * cell
  }
  for (i in seq_along(panel_rows)) {
    heatmap$heights[panel_rows[i]] <- subfamilies_per_class[[i]] * cell
  }
  # Unit conversion needs an open device
  pdf(NULL)
  width <- grid::convertWidth(sum(heatmap$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(heatmap$heights), "in", valueOnly = TRUE)
  invisible(dev.off())
}
ggsave(
  snakemake@output[["heatmap"]],
  heatmap,
  width = width,
  height = height,
  dpi = 300,
  limitsize = FALSE
)
