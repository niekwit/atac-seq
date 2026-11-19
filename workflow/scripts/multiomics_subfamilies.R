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
#   each mark changes, and whether it is in the heatmap
# - heatmap.pdf: ComplexHeatmap of the four changes (one colour scale each,
#   blue: decrease, red: increase in the treatment) for the subfamilies with
#   a significant change in accessibility or expression (at most
#   max_subfamilies, those with the smallest adjusted p-value), and the
#   subfamilies of te: profile: subfamilies (bold; not of family: and class:
#   entries), split by class. Bars: the
#   percentage of copies in the regions where the marks change.
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
max_subfamilies <- snakemake@params[["max_subfamilies"]]

atac <- read.delim(snakemake@input[["atac"]], check.names = FALSE)
atac <- atac[!is.na(atac$padj), ]
rna <- read.delim(
  snakemake@input[["rnaseq"]],
  colClasses = c(subfamily = "character")
)
rna <- rna[rna$type == "te" & rna$rna_change != "not tested", ]
signal <- read.delim(snakemake@input[["signal"]], check.names = FALSE)

# Subfamilies tested in both data sets
subfamilies <- merge(
  atac[, c(
    "subfamily",
    "family",
    "class",
    "n_copies",
    "log2FoldChange_shrunken",
    "padj",
    "change"
  )],
  data.frame(
    subfamily = rna$subfamily,
    rna_log2FoldChange = rna$log2FoldChange,
    rna_padj = rna$padj,
    rna_change = rna$rna_change
  ),
  by = "subfamily"
)
names(subfamilies)[names(subfamilies) == "log2FoldChange_shrunken"] <-
  "atac_log2FoldChange_shrunken"
names(subfamilies)[names(subfamilies) == "padj"] <- "atac_padj"
names(subfamilies)[names(subfamilies) == "change"] <- "atac_change"
cat(sprintf(
  "%d subfamilies tested in the ATAC-seq and RNA-seq data\n",
  nrow(subfamilies)
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

# Subfamilies of the profiles (te: profile: subfamilies), always plotted:
# those of subfamily and prefix entries (e.g. L1MdTf), not of family: and
# class: entries, which would fill the heatmap. In the BED file, the copies of
# an entry are followed by "#<entry> (n=...)"
profile_subfamilies <- character()
full_length <- snakemake@input[["full_length"]]
if (length(full_length) > 0) {
  pending <- character()
  for (line in readLines(full_length[[1]])) {
    if (startsWith(line, "#")) {
      entry <- sub(" \\(n=.*$", "", sub("^#", "", line))
      if (!grepl("^(family|class):", entry)) {
        profile_subfamilies <- union(profile_subfamilies, pending)
      }
      pending <- character()
    } else if (nzchar(line)) {
      pending <- union(pending, strsplit(line, "\t")[[1]][4])
    }
  }
}
subfamilies$profile <- subfamilies$subfamily %in% profile_subfamilies

# Heatmap rows: significant in ATAC-seq or RNA-seq, most significant first
significant <- subfamilies$n_copies >= min_copies &
  (subfamilies$atac_change != "unchanged" |
    subfamilies$rna_change %in% c("up", "down"))
best_padj <- pmin(
  ifelse(subfamilies$atac_change != "unchanged", subfamilies$atac_padj, 1),
  ifelse(
    subfamilies$rna_change %in% c("up", "down"),
    subfamilies$rna_padj,
    1
  ),
  na.rm = TRUE
)
ranked <- which(significant & !subfamilies$profile)
ranked <- ranked[order(best_padj[ranked])]
n_free <- max(0, max_subfamilies - sum(subfamilies$profile))
selected <- c(which(subfamilies$profile), head(ranked, n_free))
subfamilies$in_heatmap <- seq_len(nrow(subfamilies)) %in% selected
cat(sprintf(
  paste(
    "%d subfamilies (at least %d copies) significant in ATAC-seq or",
    "RNA-seq; %d in the heatmap, %d of them from te: profile: subfamilies\n"
  ),
  sum(significant),
  min_copies,
  length(selected),
  sum(subfamilies$profile)
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

# Heatmap
# -----------------------------------------------------
data <- subfamilies[subfamilies$in_heatmap, ]
if (nrow(data) == 0) {
  pdf(snakemake@output[["heatmap"]], width = 6, height = 4)
  grid.text(
    "No subfamily with a significant change in accessibility or expression"
  )
  dev.off()
  quit(save = "no")
}
# Classes in order of their number of subfamilies
class_order <- names(sort(table(data$class), decreasing = TRUE))
data$class <- factor(data$class, levels = class_order)
data <- data[order(data$class, -data$atac_log2FoldChange_shrunken), ]

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
cell_width <- unit(10, "mm")
star <- function(significant) {
  function(j, i, x, y, width, height, fill) {
    if (significant[i]) {
      grid.text("*", x, y, vjust = 0.75, gp = gpar(fontsize = 11))
    }
  }
}
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
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "white", lwd = 0.5),
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
  width = cell_width,
  height = cell_height * nrow(data),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = data$class,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  row_gap = unit(2, "mm"),
  show_row_names = FALSE,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "white", lwd = 0.5),
  cell_fun = star(data$atac_change != "unchanged"),
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
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = "white", lwd = 0.5),
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

significance <- Legend(
  labels = "padj < FDR (DESeq2)",
  type = "points",
  pch = "*",
  legend_gp = gpar(fontsize = 11),
  background = NA,
  title = sprintf("%s vs %s", treatment, control),
  title_gp = gpar(fontsize = 9)
)
height <- 1.5 + 0.18 * nrow(data) + 0.08 * nlevels(droplevels(data$class))
pdf(
  snakemake@output[["heatmap"]],
  width = 7 + 0.9 * length(bars),
  height = max(4, height)
)
draw(
  heatmaps,
  annotation_legend_list = list(significance),
  merge_legend = TRUE,
  ht_gap = unit(1, "mm")
)
dev.off()
cat("Done!\n")
