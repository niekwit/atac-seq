# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Area-proportional Venn (Euler) diagrams (eulerr) of the differential TE
# peaks of a contrast and the regions of each data set in te: correlation,
# with the GAT enrichment of the overlap. One row per data set, one column
# for the peaks more and less accessible in the treatment.
#
# Areas: peaks without overlap, regions without overlap, and the peaks that
# overlap at least one region. GAT tests the overlap in bp, against random
# placements of the peaks within the TE consensus peaks.

library(eulerr)
library(grid)
library(gridExtra)

results <- read.delim(snakemake@input[["table"]], check.names = FALSE)
title <- paste(
  snakemake@params[["treatment"]],
  "vs",
  snakemake@params[["control"]]
)
datasets <- unique(results$dataset)
directions <- c("more accessible", "less accessible")
colours <- c("more accessible" = "#D55E00", "less accessible" = "#0072B2")

gat_label <- function(r) {
  if (is.na(r$fold)) {
    return("GAT: not run")
  }
  sprintf(
    "GAT fold %.2f (observed %s bp, expected %s bp)\nq = %s",
    r$fold,
    format(r$observed_overlap_bp, big.mark = ",", scientific = FALSE),
    format(round(r$expected_overlap_bp), big.mark = ",", scientific = FALSE),
    format(signif(r$qvalue, 2))
  )
}

panel <- function(dataset, direction) {
  r <- results[results$dataset == dataset & results$peaks == direction, ]
  heading <- textGrob(
    paste0(dataset, "\n", direction, " peaks"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  if (r$n_peaks == 0) {
    body <- textGrob(paste("No", direction, "peaks"), gp = gpar(fontsize = 9))
  } else {
    sizes <- c(
      A = r$n_peaks - r$n_peaks_overlap,
      B = r$n_regions - r$n_regions_overlap,
      "A&B" = r$n_peaks_overlap
    )
    fit <- euler(sizes)
    # Set names and totals in the legend: labels next to the circles are
    # clipped when one set is much smaller than the other
    legend_labels <- c(
      sprintf(
        "%s peaks: %s (%s overlap a region)",
        direction,
        format(r$n_peaks, big.mark = ","),
        format(r$n_peaks_overlap, big.mark = ",")
      ),
      sprintf(
        "Regions: %s (%s overlap a peak)",
        format(r$n_regions, big.mark = ","),
        format(r$n_regions_overlap, big.mark = ",")
      )
    )
    body <- plot(
      fit,
      fills = list(fill = c(colours[[direction]], "grey70"), alpha = 0.6),
      edges = list(col = "grey30"),
      labels = FALSE,
      # Counts in each area: peaks without overlap, regions without overlap,
      # peaks overlapping a region
      quantities = list(
        labels = format(sizes, big.mark = ",", trim = TRUE),
        fontsize = 8
      ),
      legend = list(labels = legend_labels, side = "bottom", fontsize = 8)
    )
  }
  arrangeGrob(
    heading,
    body,
    textGrob(gat_label(r), gp = gpar(fontsize = 8)),
    heights = unit(c(0.5, 3.6, 0.45), "in")
  )
}

panels <- list()
for (dataset in datasets) {
  for (direction in directions) {
    panels[[length(panels) + 1]] <- panel(dataset, direction)
  }
}
figure <- arrangeGrob(
  grobs = panels,
  ncol = length(directions),
  top = textGrob(title, gp = gpar(fontsize = 13, fontface = "bold")),
  padding = unit(0.3, "in")
)

png(
  snakemake@output[["venn"]],
  width = 4 * length(directions),
  height = 0.5 + 4.65 * length(datasets),
  units = "in",
  res = 300
)
grid.draw(figure)
invisible(dev.off())
