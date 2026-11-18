# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Integration of the differential accessibility of repeat subfamilies (TE
# module) with their differential expression (TE rows of the RNA-seq DESeq2
# results, e.g. from TEtranscripts) of the same contrast. Accessibility of a
# subfamily is measured in two ways: the DESeq2 test of the insertions in all
# its copies (differential_subfamilies.tsv), and the differential TE peaks
# whose largest repeat overlap is a copy of the subfamily.
#
# - subfamilies.tsv: per subfamily tested in both, the ATAC-seq and RNA-seq
#   results and the number of more and less accessible TE peaks
# - correlation.tsv, scatter.pdf: Spearman correlation of the ATAC-seq and
#   RNA-seq log2 fold changes, over all subfamilies and per class
# - cdf.pdf, cdf_tests.tsv: cumulative distribution of the RNA log2 fold
#   changes of subfamilies grouped by their accessibility, Wilcoxon rank-sum
#   test of each group against the unchanged subfamilies
# - enrichment.tsv: Fisher's exact tests of up- and downregulated subfamilies
#   against more and less accessible subfamilies

library(ggplot2)
library(ggrepel)
library(cowplot)

# Load Snakemake parameters
treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
n_labels <- 20

atac <- read.delim(snakemake@input[["atac"]], check.names = FALSE)
te_peaks <- read.delim(
  snakemake@input[["peaks"]],
  check.names = FALSE,
  colClasses = c(chr = "character")
)
rna <- read.delim(
  snakemake@input[["rnaseq"]],
  colClasses = c(subfamily = "character")
)
rna <- rna[rna$type == "te" & rna$rna_change != "not tested", ]
cat(sprintf("%d subfamilies tested in the RNA-seq data\n", nrow(rna)))

# Subfamilies and TE peaks tested by DESeq2 (padj not NA), as for the RNA-seq
# data; otherwise they would count as unchanged
cat(sprintf(
  "%d subfamilies and %d TE peaks not tested by DESeq2 (padj NA) left out\n",
  sum(is.na(atac$padj)),
  sum(is.na(te_peaks$padj))
))
atac <- atac[!is.na(atac$padj), ]
te_peaks <- te_peaks[!is.na(te_peaks$padj), ]

# Subfamilies tested in both data sets
subfamilies <- merge(
  atac[, c(
    "subfamily",
    "family",
    "class",
    "baseMean",
    "log2FoldChange",
    "log2FoldChange_shrunken",
    "padj",
    "change"
  )],
  data.frame(
    subfamily = rna$subfamily,
    rna_id = rna$id,
    rna_baseMean = rna$baseMean,
    rna_log2FoldChange = rna$log2FoldChange,
    rna_padj = rna$padj,
    rna_change = rna$rna_change
  ),
  by = "subfamily"
)
atac_columns <- c(
  "baseMean",
  "log2FoldChange",
  "log2FoldChange_shrunken",
  "padj",
  "change"
)
renamed <- names(subfamilies) %in% atac_columns
names(subfamilies)[renamed] <- paste0("atac_", names(subfamilies)[renamed])
cat(sprintf(
  "%d of %d subfamilies tested in the ATAC-seq data are in the RNA-seq data\n",
  nrow(subfamilies),
  nrow(atac)
))

# Differential TE peaks per subfamily of their largest repeat overlap
count_peaks <- function(change) {
  keep <- !is.na(te_peaks$repeat_subfamily) & te_peaks$change == change
  counts <- table(te_peaks$repeat_subfamily[keep])
  out <- as.integer(counts[subfamilies$subfamily])
  out[is.na(out)] <- 0L
  out
}
in_subfamily <- table(te_peaks$repeat_subfamily[
  !is.na(te_peaks$repeat_subfamily)
])
subfamilies$n_te_peaks <- as.integer(in_subfamily[subfamilies$subfamily])
subfamilies$n_te_peaks[is.na(subfamilies$n_te_peaks)] <- 0L
subfamilies$n_te_peaks_more_accessible <- count_peaks("more accessible")
subfamilies$n_te_peaks_less_accessible <- count_peaks("less accessible")
more <- subfamilies$n_te_peaks_more_accessible > 0
less <- subfamilies$n_te_peaks_less_accessible > 0
subfamilies$te_peak_status <- ifelse(
  subfamilies$n_te_peaks == 0,
  "no peak",
  ifelse(
    more & less,
    "both",
    ifelse(
      more,
      "more accessible",
      ifelse(less, "less accessible", "unchanged")
    )
  )
)
up <- subfamilies$rna_change == "up"
down <- subfamilies$rna_change == "down"
subfamilies$concordance <- ifelse(
  (up & subfamilies$atac_change == "more accessible") |
    (down & subfamilies$atac_change == "less accessible"),
  "concordant",
  ifelse(
    (up & subfamilies$atac_change == "less accessible") |
      (down & subfamilies$atac_change == "more accessible"),
    "discordant",
    ""
  )
)
subfamilies <- subfamilies[
  order(subfamilies$rna_padj, subfamilies$atac_padj, na.last = TRUE),
]
write.table(
  subfamilies,
  snakemake@output[["subfamilies"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
if (nrow(subfamilies)) {
  cat("Subfamilies by RNA change and accessibility (subfamily counts):\n")
  print(table(subfamilies$rna_change, subfamilies$atac_change))
  cat("Subfamilies by RNA change and accessibility (TE peaks):\n")
  print(table(subfamilies$rna_change, subfamilies$te_peak_status))
}

empty_plot <- function(path, message) {
  plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = message) +
    theme_void(base_size = 12)
  ggsave(path, plot, width = 6, height = 4, bg = "white")
}

# Correlation of the ATAC-seq and RNA-seq log2 fold changes
# -----------------------------------------------------
spearman <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  if (sum(keep) < 3) {
    return(c(rho = NA, p = NA, n = sum(keep)))
  }
  test <- suppressWarnings(
    cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
  )
  c(rho = unname(test$estimate), p = test$p.value, n = sum(keep))
}
groups <- c("all", sort(unique(subfamilies$class)))
correlation <- do.call(
  rbind,
  lapply(groups, function(group) {
    keep <- if (group == "all") {
      rep(TRUE, nrow(subfamilies))
    } else {
      subfamilies$class == group
    }
    result <- spearman(
      subfamilies$atac_log2FoldChange_shrunken[keep],
      subfamilies$rna_log2FoldChange[keep]
    )
    data.frame(
      class = group,
      n_subfamilies = result[["n"]],
      spearman_rho = result[["rho"]],
      p = result[["p"]]
    )
  })
)
print(correlation)
write.table(
  correlation,
  snakemake@output[["correlation"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

if (nrow(subfamilies) == 0) {
  message <- "No subfamilies tested in both the ATAC-seq and RNA-seq data"
  empty_plot(snakemake@output[["scatter"]], message)
} else {
  plot_data <- subfamilies
  plot_data$concordance[plot_data$concordance == ""] <- "other"
  plot_data$concordance <- factor(
    plot_data$concordance,
    levels = c("other", "discordant", "concordant")
  )
  plot_data <- plot_data[order(plot_data$concordance), ]
  labelled <- head(
    subfamilies$subfamily[subfamilies$concordance != ""],
    n_labels
  )
  plot_data$label <- ifelse(
    plot_data$subfamily %in% labelled,
    plot_data$subfamily,
    ""
  )
  overall <- correlation[correlation$class == "all", ]
  scatter <- ggplot(
    plot_data,
    aes(atac_log2FoldChange_shrunken, rna_log2FoldChange, label = label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      aes(fill = class, shape = concordance, alpha = concordance),
      size = 2.5,
      colour = "black",
      stroke = 0.3
    ) +
    geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
    scale_fill_brewer(palette = "Set1") +
    # Shape and transparency share one legend, with the groups present only
    scale_shape_manual(
      values = c("other" = 21, "discordant" = 24, "concordant" = 21),
      drop = TRUE
    ) +
    scale_alpha_manual(
      values = c("other" = 0.35, "discordant" = 1, "concordant" = 1),
      drop = TRUE
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, alpha = 1)
      ),
      shape = guide_legend(order = 2, override.aes = list(fill = "grey50")),
      alpha = guide_legend(order = 2)
    ) +
    labs(
      x = "ATAC-seq shrunken log2 fold change (all copies)",
      y = "RNA-seq log2 fold change",
      fill = "Class",
      shape = "Significant in both",
      alpha = "Significant in both",
      caption = sprintf(
        "Repeat subfamilies; Spearman rho = %.2f (p = %.2g, n = %d)",
        overall$spearman_rho,
        overall$p,
        overall$n_subfamilies
      )
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
  ggsave(
    snakemake@output[["scatter"]],
    scatter,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

# Cumulative distribution of RNA log2 fold changes per accessibility group
# -----------------------------------------------------
status_levels <- c(
  "more accessible",
  "less accessible",
  "both",
  "unchanged",
  "no peak"
)
status_colours <- c(
  "more accessible" = "#dd3b3b",
  "less accessible" = "#377EB8",
  "both" = "#984EA3",
  "unchanged" = "grey30",
  "no peak" = "#cccccc"
)
measures <- c(
  "Subfamily insertions (all copies)" = "atac_change",
  "TE peaks in the subfamily" = "te_peak_status"
)
tests <- do.call(
  rbind,
  lapply(names(measures), function(measure) {
    status <- subfamilies[[measures[[measure]]]]
    reference <- subfamilies$rna_log2FoldChange[status == "unchanged"]
    do.call(
      rbind,
      lapply(status_levels, function(level) {
        values <- subfamilies$rna_log2FoldChange[status == level]
        testable <- level != "unchanged" &&
          length(values) > 0 &&
          length(reference) > 0
        data.frame(
          accessibility = measure,
          atac_status = level,
          n_subfamilies = length(values),
          median_log2FoldChange = if (length(values)) median(values) else NA,
          reference = "unchanged",
          median_difference = if (testable) {
            median(values) - median(reference)
          } else {
            NA
          },
          wilcoxon_p = if (testable) {
            suppressWarnings(wilcox.test(values, reference)$p.value)
          } else {
            NA
          }
        )
      })
    )
  })
)
tests <- tests[tests$n_subfamilies > 0 | tests$atac_status != "no peak", ]
tests$wilcoxon_padj <- p.adjust(tests$wilcoxon_p, method = "BH")
print(tests)
write.table(
  tests,
  snakemake@output[["cdf_tests"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

if (nrow(subfamilies) == 0) {
  empty_plot(snakemake@output[["cdf"]], message)
} else {
  cdf_data <- do.call(
    rbind,
    lapply(names(measures), function(measure) {
      data.frame(
        accessibility = measure,
        atac_status = subfamilies[[measures[[measure]]]],
        rna_log2FoldChange = subfamilies$rna_log2FoldChange
      )
    })
  )
  cdf_data$atac_status <- factor(cdf_data$atac_status, levels = status_levels)
  cdf_data$accessibility <- factor(
    cdf_data$accessibility,
    levels = names(measures)
  )
  counts <- aggregate(
    rna_log2FoldChange ~ accessibility + atac_status,
    data = cdf_data,
    FUN = length
  )
  counts$label <- paste0(counts$atac_status, ": ", counts$rna_log2FoldChange)
  count_text <- aggregate(
    label ~ accessibility,
    data = counts,
    FUN = paste,
    collapse = "\n"
  )
  cdf <- ggplot(cdf_data, aes(rna_log2FoldChange, colour = atac_status)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    geom_text(
      data = count_text,
      aes(x = Inf, y = -Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.2,
      size = 3
    ) +
    facet_wrap(~accessibility) +
    scale_colour_manual(values = status_colours, drop = TRUE) +
    labs(
      x = "RNA-seq log2 fold change",
      y = "Cumulative fraction of subfamilies",
      colour = "Accessibility",
      caption = "Repeat subfamilies; Wilcoxon tests in cdf_tests.tsv"
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
  ggsave(
    snakemake@output[["cdf"]],
    cdf,
    width = 10,
    height = 5,
    dpi = 300,
    bg = "white"
  )
}

# Enrichment of differentially accessible among differentially expressed
# subfamilies
# -----------------------------------------------------
fisher_row <- function(measure, in_rna, in_atac, rna_change, atac_change) {
  table_2x2 <- table(
    factor(in_rna, levels = c(TRUE, FALSE)),
    factor(in_atac, levels = c(TRUE, FALSE))
  )
  # NA when a margin is empty, where the test has no information
  informative <- all(rowSums(table_2x2) > 0) && all(colSums(table_2x2) > 0)
  test <- if (informative) {
    fisher.test(table_2x2)
  } else {
    list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  }
  data.frame(
    accessibility = measure,
    rna_change = rna_change,
    atac_change = atac_change,
    n_subfamilies = length(in_rna),
    n_rna = sum(in_rna),
    n_atac = sum(in_atac),
    n_both = sum(in_rna & in_atac),
    odds_ratio = unname(test$estimate),
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    p = test$p.value
  )
}
enrichment <- do.call(
  rbind,
  lapply(names(measures), function(measure) {
    status <- subfamilies[[measures[[measure]]]]
    do.call(
      rbind,
      lapply(c("up", "down"), function(rna_change) {
        do.call(
          rbind,
          lapply(
            c("more accessible", "less accessible"),
            function(atac_change) {
              fisher_row(
                measure,
                subfamilies$rna_change == rna_change,
                status %in% c(atac_change, "both"),
                rna_change,
                atac_change
              )
            }
          )
        )
      })
    )
  })
)
enrichment$padj <- p.adjust(enrichment$p, method = "BH")
print(enrichment)
write.table(
  enrichment,
  snakemake@output[["enrichment"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
