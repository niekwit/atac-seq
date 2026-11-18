# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Integration of differential accessibility (peaks of the ENCODE part or the
# TE module) with differential expression (RNA-seq DESeq2 results) of the
# same contrast. Each peak is assigned to the gene of its nearest TSS
# (ChIPseeker annotation), if that TSS is within the largest distance bin.
#
# - genes.tsv: per tested gene, the RNA-seq results and the number of more
#   and less accessible peaks, with the most significant peak
# - enrichment.tsv, enrichment.pdf: Fisher's exact tests of up- and
#   downregulated genes against genes with more and less accessible peaks,
#   per distance bin
# - correlation.tsv: Spearman correlation of the peak and gene log2 fold
#   changes, per distance bin, of all and of the differential peaks
# - pairs.tsv, scatter.pdf, volcano.pdf, distance.pdf: the peak-gene pairs
#   with a significant change in both (padj below the fdr of each data set
#   and |log2 fold change| > pair_lfc_threshold), by category (concordant
#   up or down, more accessible but down, less accessible but up): scatter
#   plot of the two log2 fold changes, "volcano" plot of their product
#   against the larger of the two adjusted p-values (a p-value for a change
#   in both, intersection-union test; Berger 1982), and distances of the
#   peaks to the TSS

library(ggplot2)
library(cowplot)
library(ggrepel)

# Load Snakemake parameters
treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
features <- snakemake@params[["features"]]
# Named in the captions, e.g. ENCODE consensus peaks or TE module peaks
peak_set <- snakemake@params[["peak_set"]]
bins <- sort(as.numeric(unlist(snakemake@params[["distance_bins"]])))
# Minimum |log2 fold change| of the peak and the gene for the pairs
pair_lfc <- as.numeric(snakemake@params[["pair_lfc_threshold"]])
n_labels <- 20

atac <- read.delim(
  snakemake@input[["atac"]],
  check.names = FALSE,
  colClasses = c(chr = "character", geneId = "character")
)
rna <- read.delim(
  snakemake@input[["rnaseq"]],
  colClasses = c(gene_id = "character", gene_name = "character")
)

# Genes tested in the RNA-seq data are the universe of all tests
rna <- rna[rna$type == "gene" & rna$rna_change != "not tested", ]
rna <- rna[, c(
  "gene_id",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "padj",
  "rna_change"
)]
cat(sprintf("%d genes tested in the RNA-seq data\n", nrow(rna)))
print(table(rna$rna_change))

# Distance bins: (previous bin, bin], the first from 0
kb <- function(x) {
  ifelse(x >= 1000, paste0(x / 1000, " kb"), paste0(x, " bp"))
}
bin_labels <- paste0(kb(c(0, head(bins, -1))), "-", kb(bins))
bin_labels[1] <- paste0("<= ", kb(bins[1]))
all_label <- paste0("all (<= ", kb(max(bins)), ")")

# Plot with a message only, when there is nothing to plot
message_plot <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message) +
    theme_void(base_size = 12)
}

# Peaks tested by DESeq2 (padj not NA: not removed by independent filtering
# or as outliers), of tested genes, within the largest distance bin of their
# TSS
n_untested <- sum(is.na(atac$padj))
cat(sprintf(
  "%d %s not tested by DESeq2 (padj NA) left out\n",
  n_untested,
  features
))
peaks <- atac[
  !is.na(atac$padj) & !is.na(atac$geneId) & !is.na(atac$distanceToTSS),
  c(
    "peak_id",
    "geneId",
    "distanceToTSS",
    "log2FoldChange_shrunken",
    "padj",
    "change"
  )
]
names(peaks)[names(peaks) == "geneId"] <- "gene_id"
peaks$distance <- abs(peaks$distanceToTSS)
n_annotated <- nrow(peaks)
peaks <- peaks[peaks$distance <= max(bins) & peaks$gene_id %in% rna$gene_id, ]
peaks$bin <- cut(
  peaks$distance,
  breaks = c(0, bins),
  labels = bin_labels,
  include.lowest = TRUE
)
cat(sprintf(
  "%d of %d %s assigned to a tested gene within %s of its TSS\n",
  nrow(peaks),
  nrow(atac),
  features,
  kb(max(bins))
))
cat(sprintf(
  "%d %s without gene annotation\n",
  nrow(atac) - n_annotated,
  features
))
print(table(peaks$bin, peaks$change))

# Gene-level summary
# -----------------------------------------------------
count_per_gene <- function(keep) {
  counts <- table(factor(peaks$gene_id[keep], levels = rna$gene_id))
  as.integer(counts[rna$gene_id])
}
genes <- rna
genes$n_peaks <- count_per_gene(rep(TRUE, nrow(peaks)))
genes$n_more_accessible <- count_per_gene(peaks$change == "more accessible")
genes$n_less_accessible <- count_per_gene(peaks$change == "less accessible")
genes$atac_status <- ifelse(
  genes$n_peaks == 0,
  "no peak",
  ifelse(
    genes$n_more_accessible > 0 & genes$n_less_accessible > 0,
    "both",
    ifelse(
      genes$n_more_accessible > 0,
      "more accessible",
      ifelse(genes$n_less_accessible > 0, "less accessible", "unchanged")
    )
  )
)

# Most significant peak of each gene (lowest padj, then largest change)
ordered <- peaks[
  order(peaks$padj, -abs(peaks$log2FoldChange_shrunken), na.last = TRUE),
]
best <- ordered[!duplicated(ordered$gene_id), ]
best <- data.frame(
  gene_id = best$gene_id,
  peak_id = best$peak_id,
  peak_distanceToTSS = best$distanceToTSS,
  peak_log2FoldChange_shrunken = best$log2FoldChange_shrunken,
  peak_padj = best$padj,
  peak_change = best$change
)
genes <- merge(genes, best, by = "gene_id", all.x = TRUE)

up <- genes$rna_change == "up"
down <- genes$rna_change == "down"
concordant <- (up & genes$n_more_accessible > 0) |
  (down & genes$n_less_accessible > 0)
discordant <- (up & genes$n_less_accessible > 0) |
  (down & genes$n_more_accessible > 0)
genes$concordance <- ifelse(
  concordant,
  "concordant",
  ifelse(discordant, "discordant", "")
)
genes <- genes[order(genes$padj, na.last = TRUE), ]
write.table(
  genes,
  snakemake@output[["genes"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
write.table(
  genes[genes$concordance == "concordant", ],
  snakemake@output[["concordant"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
cat("Genes by RNA change and accessibility of their peaks:\n")
print(table(genes$rna_change, genes$atac_status))

# Colours of the accessibility of the peaks
status_colours <- c(
  "more accessible" = "#dd3b3b",
  "less accessible" = "#377EB8"
)

# Enrichment of differential peaks at differentially expressed genes
# -----------------------------------------------------
fisher_row <- function(in_rna, in_atac, bin, rna_change, atac_change) {
  table_2x2 <- table(
    factor(in_rna, levels = c(TRUE, FALSE)),
    factor(in_atac, levels = c(TRUE, FALSE))
  )
  informative <- all(rowSums(table_2x2) > 0) && all(colSums(table_2x2) > 0)
  test <- if (informative) {
    fisher.test(table_2x2)
  } else {
    list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  }
  data.frame(
    distance = bin,
    rna_change = rna_change,
    atac_change = atac_change,
    n_genes = length(in_rna),
    n_rna = sum(in_rna),
    n_atac = sum(in_atac),
    n_both = sum(in_rna & in_atac),
    fraction_rna_with_atac = if (sum(in_rna)) {
      sum(in_rna & in_atac) / sum(in_rna)
    } else {
      NA
    },
    fraction_other_with_atac = if (sum(!in_rna)) {
      sum(!in_rna & in_atac) / sum(!in_rna)
    } else {
      NA
    },
    odds_ratio = unname(test$estimate),
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    p = test$p.value
  )
}
bin_sets <- c(
  setNames(lapply(bin_labels, function(b) peaks$bin == b), bin_labels),
  setNames(list(rep(TRUE, nrow(peaks))), all_label)
)
enrichment <- do.call(
  rbind,
  lapply(names(bin_sets), function(bin) {
    in_bin <- bin_sets[[bin]]
    do.call(
      rbind,
      lapply(c("up", "down"), function(rna_change) {
        do.call(
          rbind,
          lapply(
            c("more accessible", "less accessible"),
            function(atac_change) {
              with_peak <- unique(peaks$gene_id[
                in_bin & peaks$change == atac_change
              ])
              fisher_row(
                genes$rna_change == rna_change,
                genes$gene_id %in% with_peak,
                bin,
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

atac_colours <- status_colours
feature <- sub("s$", "", features)
plot_data <- enrichment[!is.na(enrichment$odds_ratio), ]
if (nrow(plot_data) == 0) {
  enrichment_plot <- message_plot(paste(
    "No test possible: no up- or downregulated genes,\nor no genes with a",
    "differential",
    feature
  ))
} else {
  plot_data$distance <- factor(plot_data$distance, levels = names(bin_sets))
  plot_data$rna_change <- factor(
    paste(plot_data$rna_change, "in RNA-seq"),
    levels = c("up in RNA-seq", "down in RNA-seq")
  )
  # Odds ratios of infinity or 0 (no genes in one cell of the table) are
  # drawn at the edge of the plot as triangles
  log2_or <- log2(plot_data$odds_ratio)
  finite <- c(log2_or, log2(plot_data$ci_low), log2(plot_data$ci_high))
  finite <- finite[is.finite(finite)]
  cap <- min(10, ceiling(max(c(abs(finite), 2))) + 1)
  plot_data$estimate <- factor(
    ifelse(is.finite(log2_or), "finite", ifelse(log2_or > 0, "infinite", "0")),
    levels = c("finite", "infinite", "0")
  )
  plot_data$log2_or <- pmin(pmax(log2_or, -cap), cap)
  plot_data$log2_low <- pmax(log2(plot_data$ci_low), -cap)
  plot_data$log2_high <- pmin(log2(plot_data$ci_high), cap)
  plot_data$label <- ifelse(
    !is.na(plot_data$padj) & plot_data$padj < 0.05,
    ifelse(
      plot_data$padj < 0.001,
      "***",
      ifelse(plot_data$padj < 0.01, "**", "*")
    ),
    ""
  )
  plot_data$label_y <- plot_data$log2_high + 0.04 * cap
  dodge <- position_dodge(width = 0.6)
  shapes <- c("finite" = 21, "infinite" = 24, "0" = 25)
  shape_labels <- c(
    "finite" = "estimate",
    "infinite" = "infinite (at the top)",
    "0" = "0 (at the bottom)"
  )
  present <- levels(droplevels(plot_data$estimate))
  enrichment_plot <- ggplot(
    plot_data,
    aes(distance, log2_or, group = atac_change)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_linerange(
      aes(ymin = log2_low, ymax = log2_high, colour = atac_change),
      position = dodge,
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = atac_change, shape = estimate),
      position = dodge,
      size = 3,
      colour = "black"
    ) +
    geom_text(
      aes(y = label_y, label = label),
      position = dodge,
      vjust = 0,
      size = 4.5
    ) +
    facet_wrap(~rna_change) +
    scale_colour_manual(values = atac_colours) +
    scale_fill_manual(values = atac_colours) +
    scale_shape_manual(
      values = shapes[present],
      labels = shape_labels[present],
      drop = TRUE
    ) +
    scale_y_continuous(limits = c(-cap, cap * 1.2)) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(shape = 21)),
      # Explains the triangles; not needed with finite estimates only
      shape = if (!identical(present, "finite")) {
        guide_legend(order = 2, override.aes = list(fill = "grey50"))
      } else {
        "none"
      }
    ) +
    labs(
      x = paste("Distance of the", feature, "to the TSS"),
      y = "log2 odds ratio (95% CI)",
      fill = paste("Genes with a", feature),
      shape = "Odds ratio",
      caption = paste0(
        peak_set,
        "\nFisher's exact test, all tested genes; ",
        "BH-adjusted p: * < 0.05, ** < 0.01, *** < 0.001"
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
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)

# Correlation of peak and gene log2 fold changes
# -----------------------------------------------------
pairs <- merge(
  peaks[!is.na(peaks$log2FoldChange_shrunken), ],
  data.frame(
    gene_id = rna$gene_id,
    gene_name = rna$gene_name,
    rna_log2FoldChange = rna$log2FoldChange,
    rna_padj = rna$padj,
    rna_change = rna$rna_change
  ),
  by = "gene_id"
)
spearman <- function(x, y) {
  if (length(x) < 3) {
    return(c(rho = NA, p = NA))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p = test$p.value)
}
correlation <- do.call(
  rbind,
  lapply(c("all", "differential"), function(which_peaks) {
    do.call(
      rbind,
      lapply(c(bin_labels, all_label), function(bin) {
        keep <- if (bin == all_label) {
          rep(TRUE, nrow(pairs))
        } else {
          pairs$bin == bin
        }
        if (which_peaks == "differential") {
          keep <- keep & pairs$change != "unchanged"
        }
        result <- spearman(
          pairs$log2FoldChange_shrunken[keep],
          pairs$rna_log2FoldChange[keep]
        )
        data.frame(
          peaks = which_peaks,
          distance = bin,
          n_peaks = sum(keep),
          n_genes = length(unique(pairs$gene_id[keep])),
          spearman_rho = result[["rho"]],
          p = result[["p"]]
        )
      })
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

# Peak-gene pairs with a significant change in both
# -----------------------------------------------------
sig <- pairs[
  pairs$change != "unchanged" &
    abs(pairs$log2FoldChange_shrunken) > pair_lfc &
    pairs$rna_change %in% c("up", "down") &
    abs(pairs$rna_log2FoldChange) > pair_lfc,
]
category_colours <- c(
  "concordant up" = "#dd3b3b",
  "concordant down" = "#377EB8",
  "more accessible, down" = "#FF7F00",
  "less accessible, up" = "#984EA3"
)
sig$category <- factor(
  ifelse(
    sig$change == "more accessible",
    ifelse(sig$rna_change == "up", "concordant up", "more accessible, down"),
    ifelse(sig$rna_change == "down", "concordant down", "less accessible, up")
  ),
  levels = names(category_colours)
)
sig$combined_score <- sig$log2FoldChange_shrunken * sig$rna_log2FoldChange
# A p-value for a change in both: the larger of the two adjusted p-values
# (intersection-union test, Berger 1982)
sig$padj_both <- pmax(sig$padj, sig$rna_padj)
sig <- sig[order(sig$padj_both, -abs(sig$combined_score)), ]
pairs_out <- data.frame(
  peak_id = sig$peak_id,
  gene_id = sig$gene_id,
  gene_name = sig$gene_name,
  distanceToTSS = sig$distanceToTSS,
  distance_bin = sig$bin,
  atac_log2FoldChange_shrunken = sig$log2FoldChange_shrunken,
  atac_padj = sig$padj,
  rna_log2FoldChange = sig$rna_log2FoldChange,
  rna_padj = sig$rna_padj,
  category = sig$category,
  combined_score = sig$combined_score,
  padj_both = sig$padj_both
)
write.table(
  pairs_out,
  snakemake@output[["pairs"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
counts <- table(sig$category)
cat(sprintf(
  "%d %s-gene pairs significant in both (|log2 fold change| > %g):\n",
  nrow(sig),
  feature,
  pair_lfc
))
print(counts)

# Legend labels with the numbers of pairs; only the categories present
present <- names(counts)[counts > 0]
category_labels <- setNames(
  sprintf("%s (n = %d)", present, as.integer(counts[present])),
  present
)
category_scale <- function(aesthetic) {
  scale_fill_manual(
    values = category_colours[present],
    breaks = present,
    labels = unname(category_labels),
    aesthetics = aesthetic
  )
}
empty_message <- paste0(
  "No ",
  feature,
  "-gene pairs with a significant change in both\n(within ",
  kb(max(bins)),
  " of the TSS, |log2 fold change| > ",
  pair_lfc,
  ")"
)
save_plot <- function(name, plot, width, height) {
  ggsave(
    snakemake@output[[name]],
    plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}
# The genes with the smallest padj_both (one label per gene)
labelled <- head(unique(sig$gene_name[!is.na(sig$gene_name)]), n_labels)
label_data <- sig[!duplicated(sig$gene_name) & sig$gene_name %in% labelled, ]

# Scatter plot of the two log2 fold changes
if (nrow(sig) == 0) {
  scatter <- message_plot(empty_message)
} else {
  rho <- spearman(sig$log2FoldChange_shrunken, sig$rna_log2FoldChange)
  scatter <- ggplot(
    sig,
    aes(log2FoldChange_shrunken, rna_log2FoldChange)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      aes(fill = category),
      shape = 21,
      size = 1.8,
      colour = "black",
      stroke = 0.2,
      alpha = 0.8
    ) +
    category_scale("fill") +
    labs(
      x = paste("ATAC-seq shrunken log2 fold change of the", feature),
      y = "RNA-seq log2 fold change of the gene",
      fill = NULL,
      caption = sprintf(
        paste0(
          "%s\n%s-gene pairs (nearest TSS within %s) with a significant ",
          "change in both; Spearman rho = %.2f (p = %.2g)"
        ),
        peak_set,
        sub("^(.)", "\\U\\1", feature, perl = TRUE),
        kb(max(bins)),
        rho[["rho"]],
        rho[["p"]]
      )
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
}
save_plot("scatter", scatter, 8, 6)

# "Volcano" plot: product of the log2 fold changes (positive: concordant)
# against the larger of the two adjusted p-values
if (nrow(sig) == 0) {
  volcano <- message_plot(empty_message)
} else {
  # padj of 0 (below the smallest double) at the top of the plot
  floor_p <- function(p) pmax(p, .Machine$double.xmin)
  volcano <- ggplot(
    sig,
    aes(combined_score, -log10(floor_p(padj_both)))
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      aes(fill = category),
      shape = 21,
      size = 1.8,
      colour = "black",
      stroke = 0.2,
      alpha = 0.8
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = gene_name),
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0
    ) +
    category_scale("fill") +
    labs(
      x = "ATAC-seq × RNA-seq log2 fold change (> 0: concordant)",
      y = "-log10 adjusted p-value (larger of the two)",
      fill = NULL,
      caption = paste0(
        peak_set,
        "\nPairs with a significant change in both; y: the larger of the ",
        "two adjusted p-values (intersection-union test)"
      )
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
}
save_plot("volcano", volcano, 8, 6)

# Distances of the peaks to the TSS, per category
if (nrow(sig) == 0) {
  distance_plot <- message_plot(empty_message)
} else {
  sig$distance_kb <- sig$distance / 1000
  distance_plot <- ggplot(sig, aes(distance_kb, fill = category)) +
    geom_histogram(
      binwidth = max(bins) / 1000 / 40,
      boundary = 0,
      colour = "black",
      linewidth = 0.2
    ) +
    facet_wrap(~category, scales = "free_y", drop = TRUE) +
    category_scale("fill") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = paste("Distance of the", feature, "to the TSS (kb)"),
      y = "Pairs",
      caption = paste0(
        peak_set,
        "\nPairs with a significant change in both"
      )
    ) +
    theme_cowplot(12) +
    theme(
      legend.position = "none",
      plot.caption = element_text(hjust = 0, size = 10)
    )
}
save_plot("distance", distance_plot, 8, 6)
