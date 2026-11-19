# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Integration of all data sets of a contrast per TE peak (TE module): each
# TE peak tested by DESeq2 with its change in accessibility, the change of
# H3K9me3 and DNA methylation in the region it overlaps most, and the
# differential expression of the gene of its nearest TSS (within the largest
# rnaseq: distance_bins, genes tested in the RNA-seq data only).
#
# - peaks.tsv: the peaks with these annotations
# - upset.pdf: UpSet plot (ComplexHeatmap) of the peaks in the sets more and
#   less accessible, the changes of each mark, and nearest gene up- and
#   downregulated: the max_combinations largest combinations (each peak
#   counted in exactly one combination)
# - enrichment.tsv, enrichment.pdf: percentage of the peaks that are more or
#   less accessible, and of the peaks with a tested gene whose gene is up- or
#   downregulated, by the combination of the first change of each mark (e.g.
#   H3K9me3 lost and hypomethylated, either one, neither); each combination
#   against neither (Fisher's exact test, Benjamini-Hochberg adjusted over
#   all tests)

library(ComplexHeatmap)
library(ggplot2)
library(cowplot)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
marks <- snakemake@params[["marks"]]
max_distance <- snakemake@params[["max_distance"]]
max_combinations <- snakemake@params[["max_combinations"]]

atac <- read.delim(
  snakemake@input[["atac"]],
  check.names = FALSE,
  colClasses = c(chr = "character", geneId = "character")
)
cat(sprintf(
  "%d TE peaks not tested by DESeq2 (padj NA) left out\n",
  sum(is.na(atac$padj))
))
peaks <- atac[
  !is.na(atac$padj),
  c(
    "peak_id",
    "chr",
    "start",
    "end",
    "log2FoldChange_shrunken",
    "padj",
    "change",
    "geneId",
    "distanceToTSS"
  )
]
for (mark in names(marks)) {
  annotated <- read.delim(snakemake@input[[mark]], check.names = FALSE)
  column <- paste0(mark, "_change")
  peaks[[column]] <- annotated[[column]][match(
    peaks$peak_id,
    annotated$peak_id
  )]
}

# Gene of the nearest TSS, if tested in the RNA-seq data and close enough
rna <- read.delim(snakemake@input[["rnaseq"]])
rna <- rna[rna$type == "gene" & rna$rna_change != "not tested", ]
near <- !is.na(peaks$geneId) &
  !is.na(peaks$distanceToTSS) &
  abs(peaks$distanceToTSS) <= max_distance &
  peaks$geneId %in% rna$gene_id
peaks$gene_id <- ifelse(near, peaks$geneId, NA)
peaks$gene_name <- rna$gene_name[match(peaks$gene_id, rna$gene_id)]
peaks$rna_log2FoldChange <- rna$log2FoldChange[match(
  peaks$gene_id,
  rna$gene_id
)]
peaks$rna_change <- rna$rna_change[match(peaks$gene_id, rna$gene_id)]
peaks$geneId <- NULL
cat(sprintf(
  "%d tested TE peaks, %d with a tested gene within %g kb of its TSS\n",
  nrow(peaks),
  sum(near),
  max_distance / 1000
))

# Sets of the UpSet plot
mark_label <- function(mark, change) {
  key <- if (change == marks[[mark]][["down"]]) "down_short" else "up_short"
  marks[[mark]][["labels"]][[key]]
}
sets <- list(
  "More accessible" = peaks$change == "more accessible",
  "Less accessible" = peaks$change == "less accessible"
)
for (mark in names(marks)) {
  for (change in unlist(marks[[mark]][["changes"]])) {
    sets[[mark_label(mark, change)]] <- peaks[[paste0(mark, "_change")]] %in%
      change
  }
}
sets[["Nearest gene up"]] <- peaks$rna_change %in% "up"
sets[["Nearest gene down"]] <- peaks$rna_change %in% "down"
for (name in names(sets)) {
  cat(sprintf("%s: %d peaks\n", name, sum(sets[[name]])))
}

# Combination of the first change of each mark (the decrease if both are
# tested), for the enrichment
first <- lapply(names(marks), function(mark) {
  change <- unlist(marks[[mark]][["changes"]])[1]
  list(
    label = mark_label(mark, change),
    member = peaks[[paste0(mark, "_change")]] %in% change
  )
})
in_first <- first[[1]]$member
in_second <- first[[2]]$member
group_levels <- c(
  sprintf("%s + %s", first[[1]]$label, tolower(first[[2]]$label)),
  sprintf("%s only", first[[1]]$label),
  sprintf("%s only", first[[2]]$label),
  "Neither"
)
peaks$mark_group <- factor(
  ifelse(
    in_first & in_second,
    group_levels[1],
    ifelse(
      in_first,
      group_levels[2],
      ifelse(in_second, group_levels[3], group_levels[4])
    )
  ),
  levels = group_levels
)
print(table(peaks$mark_group))

write.table(
  peaks,
  snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# UpSet plot
# -----------------------------------------------------
binary <- as.data.frame(lapply(sets, as.integer), check.names = FALSE)
binary <- binary[, colSums(binary) > 0, drop = FALSE]
# Peaks in no set left out (make_comb_mat would add them as the complement
# combination)
binary <- binary[rowSums(binary) > 0, , drop = FALSE]
if (ncol(binary) == 0 || nrow(binary) == 0) {
  pdf(snakemake@output[["upset"]], width = 6, height = 4)
  grid::grid.text("No TE peak in any set")
  dev.off()
} else {
  combinations <- make_comb_mat(binary, mode = "distinct")
  sizes <- comb_size(combinations)
  keep <- order(sizes, decreasing = TRUE)[seq_len(
    min(max_combinations, length(sizes))
  )]
  combinations <- combinations[keep]
  cat(sprintf(
    "UpSet plot: %d of %d combinations (%d of %d peaks in a set)\n",
    length(keep),
    length(sizes),
    sum(comb_size(combinations)),
    nrow(binary)
  ))
  # Combinations with more accessible peaks in red
  colours <- ifelse(
    grepl("^1", comb_name(combinations)) &
      "More accessible" %in% set_name(combinations),
    "#dd3b3b",
    "grey20"
  )
  upset <- UpSet(
    combinations,
    set_order = intersect(names(sets), set_name(combinations)),
    comb_order = order(comb_size(combinations), decreasing = TRUE),
    comb_col = colours,
    top_annotation = upset_top_annotation(
      combinations,
      add_numbers = TRUE,
      numbers_rot = 90,
      numbers_gp = gpar(fontsize = 7),
      annotation_name_gp = gpar(fontsize = 9),
      axis_param = list(gp = gpar(fontsize = 8))
    ),
    right_annotation = upset_right_annotation(
      combinations,
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 7),
      annotation_name_gp = gpar(fontsize = 9),
      axis_param = list(gp = gpar(fontsize = 8))
    ),
    row_names_gp = gpar(fontsize = 10),
    column_title = sprintf(
      "TE peaks (%s vs %s; %d tested)",
      treatment,
      control,
      nrow(peaks)
    ),
    column_title_gp = gpar(fontsize = 10)
  )
  pdf(
    snakemake@output[["upset"]],
    width = max(6, 2.5 + 0.25 * length(keep)),
    height = 3 + 0.25 * ncol(binary)
  )
  draw(upset)
  dev.off()
}

# Enrichment by combination of the marks
# -----------------------------------------------------
outcomes <- list(
  "More accessible" = list(
    value = peaks$change == "more accessible",
    tested = rep(TRUE, nrow(peaks))
  ),
  "Less accessible" = list(
    value = peaks$change == "less accessible",
    tested = rep(TRUE, nrow(peaks))
  ),
  "Nearest gene up" = list(
    value = peaks$rna_change %in% "up",
    tested = !is.na(peaks$rna_change)
  ),
  "Nearest gene down" = list(
    value = peaks$rna_change %in% "down",
    tested = !is.na(peaks$rna_change)
  )
)
rows <- list()
for (outcome in names(outcomes)) {
  value <- outcomes[[outcome]]$value
  tested <- outcomes[[outcome]]$tested
  reference <- tested & peaks$mark_group == "Neither"
  for (group in group_levels) {
    selected <- tested & peaks$mark_group == group
    row <- data.frame(
      outcome = outcome,
      group = group,
      n = sum(selected),
      n_outcome = sum(value[selected]),
      percent = if (any(selected)) 100 * mean(value[selected]) else NA,
      odds_ratio = NA,
      conf_low = NA,
      conf_high = NA,
      pvalue = NA
    )
    if (group != "Neither" && any(selected) && any(reference)) {
      counts <- matrix(
        c(
          sum(value[selected]),
          sum(!value[selected]),
          sum(value[reference]),
          sum(!value[reference])
        ),
        nrow = 2
      )
      test <- fisher.test(counts)
      row$odds_ratio <- unname(test$estimate)
      row$conf_low <- test$conf.int[1]
      row$conf_high <- test$conf.int[2]
      row$pvalue <- test$p.value
    }
    rows[[length(rows) + 1]] <- row
  }
}
enrichment <- do.call(rbind, rows)
enrichment$padj <- p.adjust(enrichment$pvalue, method = "BH")
write.table(
  enrichment,
  snakemake@output[["enrichment"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
print(enrichment)

stars <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "***",
      ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
    )
  )
}
enrichment$group <- factor(enrichment$group, levels = group_levels)
enrichment$outcome <- factor(enrichment$outcome, levels = names(outcomes))
enrichment$label <- stars(enrichment$padj)
group_colours <- setNames(
  c("#dd3b3b", "#377EB8", "#4DAF4A", "#cccccc"),
  group_levels
)
plotted <- enrichment[!is.na(enrichment$percent), ]
if (nrow(plotted) == 0) {
  pdf(snakemake@output[["enrichment_plot"]], width = 6, height = 4)
  grid::grid.text("No tested TE peak")
  dev.off()
} else {
  plot <- ggplot(
    plotted,
    aes(group, percent, fill = group)
  ) +
    geom_col(colour = "black", width = 0.7) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.5) +
    geom_text(aes(y = 0, label = paste0("n=", n)), vjust = 1.4, size = 2.8) +
    facet_wrap(~outcome, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colours, guide = "none") +
    scale_x_discrete(labels = function(x) sub(" (only|\\+ )", "\n\\1", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(
      x = NULL,
      y = sprintf("TE peaks (%%; %s vs %s)", treatment, control)
    ) +
    theme_cowplot(12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 30)
    )
  ggsave(
    snakemake@output[["enrichment_plot"]],
    plot,
    width = 10,
    height = 4.5,
    bg = "white"
  )
}
cat("Done!\n")
