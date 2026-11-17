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
# - cdf.pdf, cdf_tests.tsv: cumulative distribution of the RNA log2 fold
#   changes of genes grouped by the accessibility of their peaks, Wilcoxon
#   rank-sum test of each group against the genes with only unchanged peaks
# - enrichment.tsv, enrichment.pdf: Fisher's exact tests of up- and
#   downregulated genes against genes with more and less accessible peaks,
#   per distance bin
# - correlation.tsv, scatter.pdf: Spearman correlation of the peak and gene
#   log2 fold changes, per distance bin

library(ggplot2)

# Load Snakemake parameters
treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
features <- snakemake@params[["features"]]
bins <- sort(as.numeric(unlist(snakemake@params[["distance_bins"]])))
title <- paste(treatment, "vs", control)

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

# Peaks of tested genes, within the largest distance bin of their TSS
peaks <- atac[
  !is.na(atac$geneId) & !is.na(atac$distanceToTSS),
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
cat(sprintf("%d %s without gene annotation\n", nrow(atac) - n_annotated, features))
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
  "more accessible" = "#D55E00",
  "less accessible" = "#0072B2",
  "both" = "#CC79A7",
  "unchanged" = "grey40",
  "no peak" = "grey75"
)
reference <- genes$log2FoldChange[genes$atac_status == "unchanged"]
tests <- do.call(
  rbind,
  lapply(status_levels, function(status) {
    values <- genes$log2FoldChange[genes$atac_status == status]
    testable <- status != "unchanged" && length(values) > 0 &&
      length(reference) > 0
    data.frame(
      atac_status = status,
      n_genes = length(values),
      median_log2FoldChange = if (length(values)) median(values) else NA,
      reference = "unchanged",
      median_difference = if (testable) {
        median(values) - median(reference)
      } else {
        NA
      },
      wilcoxon_p = if (testable) {
        wilcox.test(values, reference)$p.value
      } else {
        NA
      }
    )
  })
)
tests$wilcoxon_padj <- p.adjust(tests$wilcoxon_p, method = "BH")
print(tests)
write.table(
  tests,
  snakemake@output[["cdf_tests"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

present <- status_levels[status_levels %in% genes$atac_status]
legend_labels <- setNames(
  sprintf(
    "%s (n = %d%s)",
    present,
    tests$n_genes[match(present, tests$atac_status)],
    ifelse(
      is.na(tests$wilcoxon_p[match(present, tests$atac_status)]),
      "",
      sprintf(", p = %.2g", tests$wilcoxon_p[match(present, tests$atac_status)])
    )
  ),
  present
)
cdf_data <- genes[genes$atac_status %in% present, ]
cdf_data$atac_status <- factor(cdf_data$atac_status, levels = present)
limits <- quantile(cdf_data$log2FoldChange, c(0.005, 0.995), na.rm = TRUE)
cdf <- ggplot(cdf_data, aes(log2FoldChange, colour = atac_status)) +
  stat_ecdf(geom = "step") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  coord_cartesian(xlim = limits) +
  scale_colour_manual(
    values = status_colours,
    breaks = present,
    labels = unname(legend_labels[present])
  ) +
  labs(
    title = title,
    subtitle = paste0(
      "Genes by the accessibility of their ",
      features,
      " (nearest TSS within ",
      kb(max(bins)),
      "); p: Wilcoxon vs unchanged"
    ),
    x = "RNA-seq log2 fold change",
    y = "Cumulative fraction of genes",
    colour = paste0(features, " of the gene")
  ) +
  theme_bw(base_size = 12)
ggsave(snakemake@output[["cdf"]], cdf, width = 8, height = 5, dpi = 300)

# Enrichment of differential peaks at differentially expressed genes
# -----------------------------------------------------
fisher_row <- function(in_rna, in_atac, bin, rna_change, atac_change) {
  table_2x2 <- table(
    factor(in_rna, levels = c(TRUE, FALSE)),
    factor(in_atac, levels = c(TRUE, FALSE))
  )
  test <- tryCatch(
    fisher.test(table_2x2),
    error = function(e) list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
  )
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
          lapply(c("more accessible", "less accessible"), function(atac_change) {
            with_peak <- unique(peaks$gene_id[in_bin & peaks$change == atac_change])
            fisher_row(
              genes$rna_change == rna_change,
              genes$gene_id %in% with_peak,
              bin,
              rna_change,
              atac_change
            )
          })
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

atac_colours <- c("more accessible" = "#D55E00", "less accessible" = "#0072B2")
plot_data <- enrichment
plot_data$distance <- factor(plot_data$distance, levels = names(bin_sets))
plot_data$rna_change <- factor(
  paste(plot_data$rna_change, "in RNA-seq"),
  levels = c("up in RNA-seq", "down in RNA-seq")
)
# Odds ratios of 0 or infinity (no or all genes with a peak) are not shown
finite <- is.finite(log2(plot_data$odds_ratio))
plot_data$log2_or <- ifelse(finite, log2(plot_data$odds_ratio), NA)
plot_data$log2_low <- ifelse(finite, pmax(log2(plot_data$ci_low), -10), NA)
plot_data$log2_high <- ifelse(finite, pmin(log2(plot_data$ci_high), 10), NA)
plot_data$label <- ifelse(
  !is.na(plot_data$padj) & plot_data$padj < 0.05,
  ifelse(plot_data$padj < 0.001, "***", ifelse(plot_data$padj < 0.01, "**", "*")),
  ""
)
dodge <- position_dodge(width = 0.5)
enrichment_plot <- ggplot(
  plot_data,
  aes(distance, log2_or, colour = atac_change, label = label)
) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_pointrange(
    aes(ymin = log2_low, ymax = log2_high),
    position = dodge,
    na.rm = TRUE
  ) +
  geom_text(
    aes(y = log2_high),
    position = dodge,
    vjust = -0.3,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  facet_wrap(~rna_change) +
  scale_colour_manual(values = atac_colours) +
  labs(
    title = title,
    subtitle = paste0(
      "Fisher's exact test, all tested genes; * BH-adjusted p < 0.05, ",
      "** < 0.01, *** < 0.001"
    ),
    x = paste("Distance of the", sub("s$", "", features), "to the TSS"),
    y = "log2 odds ratio (95% CI)",
    colour = paste("Gene with a", sub("s$", "", features))
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  snakemake@output[["enrichment_plot"]],
  enrichment_plot,
  width = 9,
  height = 5,
  dpi = 300
)

# Correlation of peak and gene log2 fold changes
# -----------------------------------------------------
pairs <- merge(
  peaks[!is.na(peaks$log2FoldChange_shrunken), ],
  rna[, c("gene_id", "log2FoldChange", "rna_change")],
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
          pairs$log2FoldChange[keep]
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

pairs$class <- "other"
pairs$class[
  (pairs$change == "more accessible" & pairs$rna_change == "up") |
    (pairs$change == "less accessible" & pairs$rna_change == "down")
] <- "concordant"
pairs$class[
  (pairs$change == "more accessible" & pairs$rna_change == "down") |
    (pairs$change == "less accessible" & pairs$rna_change == "up")
] <- "discordant"
pairs$class <- factor(pairs$class, levels = c("other", "discordant", "concordant"))
pairs <- pairs[order(pairs$class), ]
class_colours <- c(
  "concordant" = "#009E73",
  "discordant" = "#E69F00",
  "other" = "grey70"
)
facet_text <- correlation[correlation$peaks == "all" & correlation$distance != all_label, ]
facet_text$bin <- factor(facet_text$distance, levels = bin_labels)
facet_text$label <- sprintf(
  "rho = %.2f\nn = %d",
  facet_text$spearman_rho,
  facet_text$n_peaks
)
scatter <- ggplot(
  pairs,
  aes(log2FoldChange_shrunken, log2FoldChange, colour = class)
) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_text(
    data = facet_text,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 3.5
  ) +
  facet_wrap(~bin, drop = FALSE) +
  scale_colour_manual(values = class_colours, drop = FALSE) +
  labs(
    title = title,
    subtitle = paste(
      "Each", sub("s$", "", features),
      "against its nearest gene, by distance to the TSS; Spearman rho"
    ),
    x = paste("ATAC-seq shrunken log2 fold change of the", sub("s$", "", features)),
    y = "RNA-seq log2 fold change of the gene",
    colour = "Significant in both"
  ) +
  theme_bw(base_size = 12)
ggsave(snakemake@output[["scatter"]], scatter, width = 9, height = 7, dpi = 300)
