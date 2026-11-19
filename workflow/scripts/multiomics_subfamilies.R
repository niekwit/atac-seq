# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Integration of all data sets of a contrast per repeat subfamily (TE
# module): differential accessibility (DESeq2 of the insertions in all
# copies), differential expression (TE rows of the RNA-seq DESeq2 results,
# e.g. from TEtranscripts), and the change of H3K9me3 and DNA methylation
# over all copies (from the BigWig files of the control and the treatment).
#
# - subfamilies.tsv: every subfamily tested in the ATAC-seq and RNA-seq data,
#   with the four changes, the percentage of its copies in the regions where
#   each mark changes, and whether it is in the heatmap of the significant
#   subfamilies
# - heatmap.pdf: ComplexHeatmap of the four changes (one colour scale each,
#   blue: decrease, red: increase in the treatment) for the subfamilies with
#   a significant change in accessibility or expression (at most
#   max_per_class per class, those with the smallest adjusted p-value),
#   split by class (profile subfamilies in bold). Bars: the percentage of
#   copies in the regions where the marks change.
# - profile_heatmap.pdf: the same for the subfamilies of te: profile:
#   subfamilies (not of family: and class: entries; also if not tested in a
#   data set, grey)
#
# H3K9me3: log2 of the ratio of the mean signal per base of the copies in the
# treatment and the control, with a pseudocount of 1% of the median signal of
# the subfamilies. DNA methylation: the mean methylation of the covered CpGs
# of the copies in the treatment minus the control.

library(ComplexHeatmap)
library(circlize)
library(grid)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
marks <- snakemake@params[["marks"]]
min_copies <- snakemake@params[["min_copies"]]
max_per_class <- snakemake@params[["max_per_class"]]

profile_entries <- unlist(snakemake@params[["profile_entries"]])
# All counted subfamilies (resources/te/subfamilies.tsv)
annotation <- read.delim(
  snakemake@input[["annotation"]],
  check.names = FALSE,
  colClasses = c(subfamily = "character")
)
atac <- read.delim(snakemake@input[["atac"]], check.names = FALSE)
rna <- read.delim(
  snakemake@input[["rnaseq"]],
  colClasses = c(subfamily = "character")
)
rna <- rna[rna$type == "te" & rna$rna_change != "not tested", ]
signal <- read.delim(snakemake@input[["signal"]], check.names = FALSE)

# Subfamilies of the profiles (te: profile: subfamilies), own heatmap: a
# subfamily, or the subfamilies <entry>_* if no subfamily has that name (as
# for the profiles); not family: and class: entries, which would be too
# many. Also those without full-length copies or data in a data set
# (missing values are grey)
profile_subfamilies <- character()
for (entry in profile_entries) {
  if (grepl("^(family|class):", entry)) {
    next
  }
  members <- if (entry %in% annotation$subfamily) {
    entry
  } else {
    annotation$subfamily[startsWith(annotation$subfamily, paste0(entry, "_"))]
  }
  if (length(members) == 0) {
    cat(sprintf("te: profile: subfamilies: %s: no subfamily\n", entry))
  }
  profile_subfamilies <- union(profile_subfamilies, members)
}

# Subfamilies tested in both data sets (ATAC-seq padj not NA), and those of
# the profiles
atac_tested <- atac$subfamily[!is.na(atac$padj)]
keep <- (annotation$subfamily %in%
  atac_tested &
  annotation$subfamily %in% rna$subfamily) |
  annotation$subfamily %in% profile_subfamilies
subfamilies <- annotation[keep, c("subfamily", "family", "class", "n_copies")]
atac <- atac[!is.na(atac$padj), ]
match_atac <- match(subfamilies$subfamily, atac$subfamily)
subfamilies$atac_log2FoldChange_shrunken <- atac$log2FoldChange_shrunken[
  match_atac
]
subfamilies$atac_padj <- atac$padj[match_atac]
subfamilies$atac_change <- atac$change[match_atac]
match_rna <- match(subfamilies$subfamily, rna$subfamily)
subfamilies$rna_log2FoldChange <- rna$log2FoldChange[match_rna]
subfamilies$rna_padj <- rna$padj[match_rna]
subfamilies$rna_change <- rna$rna_change[match_rna]
subfamilies$profile <- subfamilies$subfamily %in% profile_subfamilies
cat(sprintf(
  paste(
    "%d subfamilies tested in the ATAC-seq and RNA-seq data; %d profile",
    "subfamilies (%d not tested in both)\n"
  ),
  sum(!is.na(match_atac) & !is.na(match_rna)),
  sum(subfamilies$profile),
  sum(subfamilies$profile & (is.na(match_atac) | is.na(match_rna)))
))

# H3K9me3 and DNA methylation over the copies
pseudocount <- function(values) {
  positive <- values[is.finite(values) & values > 0]
  if (length(positive) == 0) 1 else 0.01 * median(positive)
}
h3 <- c(signal$h3k9me3_control, signal$h3k9me3_treatment)
p <- pseudocount(h3)
cat(sprintf("H3K9me3 pseudocount: %g\n", p))
signal$h3k9me3_log2ratio <- log2(
  (signal$h3k9me3_treatment + p) / (signal$h3k9me3_control + p)
)
signal$methylation_difference <- signal$methylation_treatment -
  signal$methylation_control
subfamilies <- merge(
  subfamilies,
  signal[, c(
    "subfamily",
    "h3k9me3_control",
    "h3k9me3_treatment",
    "h3k9me3_log2ratio",
    "methylation_control",
    "methylation_treatment",
    "methylation_difference",
    "methylation_control_covered_bp",
    "methylation_treatment_covered_bp"
  )],
  by = "subfamily",
  all.x = TRUE
)

# Percentage of the copies in the regions where each mark changes
bars <- list()
for (mark in names(marks)) {
  copies <- read.delim(snakemake@input[[mark]], check.names = FALSE)
  for (change in unlist(marks[[mark]][["changes"]])) {
    column <- sprintf("percent_copies_%s_%s", mark, change)
    n <- copies[[paste0("n_copies_", change)]]
    percent <- setNames(100 * n / copies$n_copies, copies$subfamily)
    subfamilies[[column]] <- unname(percent[subfamilies$subfamily])
    label <- marks[[mark]][["labels"]][[paste0(
      if (change == marks[[mark]][["down"]]) "down" else "up",
      "_short"
    )]]
    bars[[column]] <- label
  }
}


# Heatmap rows: significant in ATAC-seq or RNA-seq, at most max_per_class
# per class, most significant first
atac_significant <- subfamilies$atac_change %in%
  c("more accessible", "less accessible")
significant <- subfamilies$n_copies >= min_copies &
  (atac_significant | subfamilies$rna_change %in% c("up", "down"))
best_padj <- pmin(
  ifelse(atac_significant, subfamilies$atac_padj, 1),
  ifelse(
    subfamilies$rna_change %in% c("up", "down"),
    subfamilies$rna_padj,
    1
  ),
  na.rm = TRUE
)
ranked <- which(significant)
ranked <- ranked[order(best_padj[ranked])]
selected <- unlist(lapply(
  split(ranked, subfamilies$class[ranked]),
  head,
  max_per_class
))
subfamilies$in_heatmap <- seq_len(nrow(subfamilies)) %in% selected
cat(sprintf(
  paste(
    "%d subfamilies (at least %d copies) significant in ATAC-seq or",
    "RNA-seq; %d in the heatmap (at most %d per class)\n"
  ),
  sum(significant),
  min_copies,
  length(selected),
  max_per_class
))

subfamilies <- subfamilies[
  order(!subfamilies$in_heatmap, -subfamilies$atac_log2FoldChange_shrunken),
]
write.table(
  subfamilies,
  snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# Heatmaps
# -----------------------------------------------------
# Diverging colours (blue: decrease, red: increase), limit: the largest
# absolute value up to the 95th percentile
diverging <- function(values) {
  limit <- quantile(abs(values), 0.95, na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) {
    limit <- 1
  }
  colorRamp2(c(-limit, 0, limit), c("#377EB8", "white", "#dd3b3b"))
}
cell_height <- unit(4.5, "mm")
# Square cells
cell_width <- cell_height
star <- function(significant) {
  function(j, i, x, y, width, height, fill) {
    if (significant[i]) {
      grid.text("*", x, y, vjust = 0.75, gp = gpar(fontsize = 11))
    }
  }
}
significance <- Legend(
  labels = "padj < FDR (DESeq2)",
  type = "points",
  pch = "*",
  legend_gp = gpar(fontsize = 11),
  background = NA,
  title = sprintf("%s vs %s", treatment, control),
  title_gp = gpar(fontsize = 9)
)

# Heatmap of the four changes of the subfamilies in data, split by class
# (classes with the most subfamilies first), sorted by the ATAC-seq fold
# change; profile subfamilies in bold
plot_heatmap <- function(data, file, empty) {
  if (nrow(data) == 0) {
    pdf(file, width = 6, height = 4)
    grid.text(empty)
    dev.off()
    return(invisible())
  }
  class_order <- names(sort(table(data$class), decreasing = TRUE))
  data$class <- factor(data$class, levels = class_order)
  data <- data[order(data$class, -data$atac_log2FoldChange_shrunken), ]

  column <- function(values, name, label, title, significant = NULL) {
    Heatmap(
      matrix(values, ncol = 1, dimnames = list(data$subfamily, label)),
      name = name,
      col = diverging(values),
      na_col = "grey90",
      width = cell_width,
      height = cell_height * nrow(data),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      column_names_rot = 90,
      column_names_gp = gpar(fontsize = 10),
      rect_gp = gpar(col = "black", lwd = 0.5),
      cell_fun = if (is.null(significant)) NULL else star(significant),
      heatmap_legend_param = list(title = title, title_gp = gpar(fontsize = 9))
    )
  }

  bar_annotations <- lapply(names(bars), function(column) {
    anno_barplot(
      data[[column]],
      which = "row",
      gp = gpar(fill = "grey40", col = NA),
      border = FALSE,
      width = unit(18, "mm"),
      axis_param = list(gp = gpar(fontsize = 7), side = "bottom")
    )
  })
  names(bar_annotations) <- names(bars)
  right <- do.call(
    rowAnnotation,
    c(
      bar_annotations,
      list(
        annotation_label = sprintf("%s\n(%% copies)", unlist(bars)),
        annotation_name_side = "top",
        annotation_name_rot = 0,
        annotation_name_gp = gpar(fontsize = 8),
        gap = unit(4, "mm")
      )
    )
  )

  # The first heatmap splits the rows by class; row names after the last one
  heatmaps <- Heatmap(
    matrix(
      data$atac_log2FoldChange_shrunken,
      ncol = 1,
      dimnames = list(data$subfamily, "ATAC-seq")
    ),
    name = "atac",
    col = diverging(data$atac_log2FoldChange_shrunken),
    na_col = "grey90",
    width = cell_width,
    height = cell_height * nrow(data),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = data$class,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 10),
    row_gap = unit(2, "mm"),
    show_row_names = FALSE,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "black", lwd = 0.5),
    cell_fun = star(
      data$atac_change %in% c("more accessible", "less accessible")
    ),
    heatmap_legend_param = list(
      title = "ATAC-seq\nlog2 FC",
      title_gp = gpar(fontsize = 9)
    )
  ) +
    column(
      data$rna_log2FoldChange,
      "rna",
      "RNA-seq",
      "RNA-seq\nlog2 FC",
      data$rna_change %in% c("up", "down")
    ) +
    column(
      data$h3k9me3_log2ratio,
      "h3k9me3",
      "H3K9me3",
      "H3K9me3\nlog2 ratio"
    ) +
    Heatmap(
      matrix(
        data$methylation_difference,
        ncol = 1,
        dimnames = list(data$subfamily, "CpG methylation")
      ),
      name = "methylation",
      col = diverging(data$methylation_difference),
      na_col = "grey90",
      width = cell_width,
      height = cell_height * nrow(data),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_rot = 90,
      column_names_gp = gpar(fontsize = 10),
      rect_gp = gpar(col = "black", lwd = 0.5),
      right_annotation = right,
      show_row_names = TRUE,
      row_names_gp = gpar(
        fontsize = 9,
        fontface = ifelse(data$profile, "bold", "plain")
      ),
      heatmap_legend_param = list(
        title = "Methylation\ndifference",
        title_gp = gpar(fontsize = 9)
      )
    )

  height <- 1.5 + 0.18 * nrow(data) + 0.08 * nlevels(droplevels(data$class))
  pdf(file, width = 6 + 0.9 * length(bars), height = max(4, height))
  draw(
    heatmaps,
    annotation_legend_list = list(significance),
    merge_legend = TRUE,
    ht_gap = unit(1, "mm")
  )
  dev.off()
}

plot_heatmap(
  subfamilies[subfamilies$in_heatmap, ],
  snakemake@output[["heatmap"]],
  "No subfamily with a significant change in accessibility or expression"
)
plot_heatmap(
  subfamilies[subfamilies$profile, ],
  snakemake@output[["profile_heatmap"]],
  "No subfamily of te: profile: subfamilies"
)
cat("Done!\n")
