# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Integration of the differential accessibility of the peaks of a contrast
# with the differential regions of an epigenomic mark (H3K9me3, DNA
# methylation) of the same contrast (epi_annotate_peaks: the region each peak
# overlaps most).
#
# - enrichment.tsv, enrichment.pdf: percentage of the more, less and
#   unchanged accessible peaks in regions where the mark decreases (or
#   increases), with Fisher's exact tests of the more and less accessible
#   against the unchanged peaks
# - distribution.pdf, distribution_tests.tsv: cumulative distribution of the
#   ATAC-seq log2 fold changes of the peaks by the change of the mark in their
#   region, Wilcoxon rank-sum test of each group against the peaks outside
#   the regions

library(ggplot2)
library(cowplot)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
peak_set <- snakemake@params[["peak_set"]]
# Mark (h3k9me3, methylation): prefix of its columns in the input and output
prefix <- snakemake@params[["prefix"]]
# Change codes of a decrease and an increase in the treatment (lost/gained,
# hypo/hyper), and plot labels
down <- snakemake@params[["down"]]
up <- snakemake@params[["up"]]
labels <- snakemake@params[["labels"]]

peaks <- read.delim(
  snakemake@input[["peaks"]],
  colClasses = c(chr = "character")
)
# Generic names of the columns of the mark (mark_change, ...)
names(peaks) <- sub(paste0("^", prefix, "_"), "mark_", names(peaks))
cat(sprintf("%d peaks tested by DESeq2\n", nrow(peaks)))
print(table(peaks$change, peaks$mark_change))

# Writes a table with the column names of the mark (e.g. mark_change)
write_tsv <- function(data, path) {
  names(data) <- sub("^mark_", paste0(prefix, "_"), names(data))
  write.table(data, path, quote = FALSE, sep = "\t", row.names = FALSE)
}

message_plot <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message) +
    theme_void(base_size = 12)
}

stars <- function(p) {
  ifelse(
    is.na(p) | p >= 0.05,
    "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  )
}

atac_levels <- c("more accessible", "less accessible", "unchanged")
atac_colours <- c(
  "more accessible" = "#dd3b3b",
  "less accessible" = "#377EB8",
  "unchanged" = "#cccccc"
)
mark_levels <- c(down, up, "unchanged", "outside")
mark_labels <- setNames(
  c(
    labels[["down"]],
    labels[["up"]],
    labels[["unchanged"]],
    labels[["outside"]]
  ),
  mark_levels
)
mark_colours <- setNames(
  c("#377EB8", "#dd3b3b", "grey30", "#cccccc"),
  mark_levels
)
# Changes of the mark in the data (BED file: one direction)
changes <- intersect(c(down, up), snakemake@params[["changes"]])

# Enrichment of differential peaks in regions where the mark changes
# -----------------------------------------------------
reference <- peaks$change == "unchanged"
enrichment <- do.call(
  rbind,
  lapply(changes, function(mark_change) {
    in_regions <- peaks$mark_change == mark_change
    do.call(
      rbind,
      lapply(atac_levels, function(atac_change) {
        group <- peaks$change == atac_change
        testable <- atac_change != "unchanged" && any(group) && any(reference)
        test <- list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
        if (testable) {
          keep <- group | reference
          table_2x2 <- table(
            factor(group[keep], levels = c(TRUE, FALSE)),
            factor(in_regions[keep], levels = c(TRUE, FALSE))
          )
          if (all(colSums(table_2x2) > 0)) {
            test <- fisher.test(table_2x2)
          }
        }
        data.frame(
          mark_change = mark_change,
          atac_change = atac_change,
          n_peaks = sum(group),
          n_in_regions = sum(group & in_regions),
          percent_in_regions = if (any(group)) {
            100 * mean(in_regions[group])
          } else {
            NA
          },
          reference = if (atac_change == "unchanged") NA else "unchanged",
          odds_ratio = unname(test$estimate),
          ci_low = test$conf.int[1],
          ci_high = test$conf.int[2],
          p = test$p.value
        )
      })
    )
  })
)
enrichment$padj <- p.adjust(enrichment$p, method = "BH")
print(enrichment)
write_tsv(enrichment, snakemake@output[["enrichment"]])

plot_data <- enrichment[enrichment$n_peaks > 0, ]
if (nrow(plot_data) == 0) {
  enrichment_plot <- message_plot("No peaks tested by DESeq2")
} else {
  plot_data$atac_change <- factor(plot_data$atac_change, levels = atac_levels)
  plot_data$mark_change <- factor(
    mark_labels[plot_data$mark_change],
    levels = mark_labels[changes]
  )
  plot_data$label <- paste0(
    plot_data$n_in_regions,
    "/",
    plot_data$n_peaks,
    "\n",
    stars(plot_data$padj)
  )
  enrichment_plot <- ggplot(
    plot_data,
    aes(atac_change, percent_in_regions, fill = atac_change)
  ) +
    geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = label), vjust = -0.2, size = 3.5, lineheight = 0.9) +
    facet_wrap(~mark_change) +
    scale_fill_manual(values = atac_colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      x = paste("Accessibility in", treatment, "vs", control),
      y = paste0("Peaks in the ", labels[["regions"]], " (%)"),
      caption = paste0(
        peak_set,
        "\nFisher's exact test against the unchanged peaks",
        "\nBH-adjusted p: * < 0.05, ** < 0.01, *** < 0.001"
      )
    ) +
    theme_cowplot(12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0, size = 10)
    )
}
ggsave(
  snakemake@output[["enrichment_plot"]],
  enrichment_plot,
  width = 3 + 2.5 * max(1, length(changes)),
  height = 5,
  dpi = 300,
  bg = "white"
)

# Distribution of the ATAC-seq log2 fold changes by change of the mark
# -----------------------------------------------------
outside <- peaks$log2FoldChange_shrunken[peaks$mark_change == "outside"]
tests <- do.call(
  rbind,
  lapply(mark_levels, function(level) {
    values <- peaks$log2FoldChange_shrunken[peaks$mark_change == level]
    testable <- level != "outside" &&
      length(values) > 0 &&
      length(outside) > 0
    data.frame(
      mark_change = level,
      n_peaks = length(values),
      median_log2FoldChange_shrunken = if (length(values)) {
        median(values)
      } else {
        NA
      },
      reference = if (level == "outside") NA else "outside",
      median_difference = if (testable) {
        median(values) - median(outside)
      } else {
        NA
      },
      wilcoxon_p = if (testable) {
        suppressWarnings(wilcox.test(values, outside)$p.value)
      } else {
        NA
      }
    )
  })
)
tests <- tests[tests$n_peaks > 0, ]
tests$wilcoxon_padj <- p.adjust(tests$wilcoxon_p, method = "BH")
print(tests)
write_tsv(tests, snakemake@output[["distribution_tests"]])

if (nrow(peaks) == 0) {
  distribution <- message_plot("No peaks tested by DESeq2")
} else {
  cdf_data <- peaks[, c("mark_change", "log2FoldChange_shrunken")]
  present <- mark_levels[mark_levels %in% cdf_data$mark_change]
  legend_labels <- sapply(present, function(level) {
    row <- tests[tests$mark_change == level, ]
    sprintf(
      "%s (n = %d)%s",
      mark_labels[[level]],
      row$n_peaks,
      if (is.na(row$wilcoxon_padj)) {
        ""
      } else {
        # Below the smallest double, p is 0
        sprintf(
          ", p %s",
          if (row$wilcoxon_padj == 0) {
            "< 1e-300"
          } else {
            sprintf("= %.2g", row$wilcoxon_padj)
          }
        )
      }
    )
  })
  cdf_data$mark_change <- factor(cdf_data$mark_change, levels = present)
  distribution <- ggplot(
    cdf_data,
    aes(log2FoldChange_shrunken, colour = mark_change)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    scale_colour_manual(
      values = mark_colours[present],
      labels = legend_labels
    ) +
    labs(
      x = paste(
        "ATAC-seq shrunken log2 fold change,",
        treatment,
        "vs",
        control
      ),
      y = "Cumulative fraction of peaks",
      colour = NULL,
      caption = paste0(
        peak_set,
        "\nWilcoxon rank-sum test against the peaks outside the ",
        labels[["regions"]],
        " (BH-adjusted p)"
      )
    ) +
    theme_cowplot(12) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      plot.caption = element_text(hjust = 0, size = 10)
    )
}
ggsave(
  snakemake@output[["distribution"]],
  distribution,
  width = 7,
  height = 6.5,
  dpi = 300,
  bg = "white"
)
