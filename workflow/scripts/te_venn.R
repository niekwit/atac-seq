# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Area-proportional Venn (Euler) diagrams (eulerr) of the differential TE
# peaks of a contrast and the regions of each data set in te: correlation,
# with the GAT enrichment of the overlap. One row per data set, one column
# for the peaks more and less accessible in the treatment.
#
# Circles: the peaks and the regions of the data set, each merged within its
# own set (overlapping intervals only), labelled with their totals; overlap:
# the intersections (overlapping peak-region pairs; counted by te_gat.py).
# All panels share one scale, so a data set has the same size in each panel.
# GAT tests the overlap in bp, against random placements of the peaks within
# the TE consensus peaks.

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

# Euler fit per panel: areas of peaks only, regions only and intersections.
# Intersections are pairs, so they can exceed the smaller set (one peak
# overlapping many regions): then capped at the smaller set for the drawing.
fit_panel <- function(r) {
  both <- min(r$n_intersections, r$n_peaks, r$n_regions)
  euler(c(A = r$n_peaks - both, B = r$n_regions - both, "A&B" = both))
}

# Width or height of a fitted diagram, whichever is larger: plot.euler scales
# this to the viewport, so viewports scaled by it share one scale
extent <- function(fit) {
  e <- fit$ellipses
  max(
    diff(range(c(e$h - e$a, e$h + e$a))),
    diff(range(c(e$k - e$b, e$k + e$b)))
  )
}

rows <- expand.grid(
  direction = directions,
  dataset = datasets,
  stringsAsFactors = FALSE
)
fits <- lapply(seq_len(nrow(rows)), function(i) {
  r <- results[
    results$dataset == rows$dataset[i] & results$peaks == rows$direction[i],
  ]
  if (r$n_peaks == 0) NULL else fit_panel(r)
})
extents <- vapply(fits, function(f) if (is.null(f)) 0 else extent(f), 0)
max_extent <- max(extents)

panel <- function(i) {
  dataset <- rows$dataset[i]
  direction <- rows$direction[i]
  r <- results[results$dataset == dataset & results$peaks == direction, ]
  heading <- textGrob(
    paste0(dataset, "\n", direction, " peaks"),
    gp = gpar(fontsize = 10, fontface = "bold")
  )
  if (r$n_peaks == 0) {
    body <- textGrob(paste("No", direction, "peaks"), gp = gpar(fontsize = 9))
  } else {
    fit <- fits[[i]]
    # Totals in the circles, intersections in the overlap (quantities are
    # placed in the areas A only, B only and A&B)
    labels <- format(
      c(r$n_peaks, r$n_regions, r$n_intersections),
      big.mark = ",",
      trim = TRUE
    )
    diagram <- plot(
      fit,
      fills = list(fill = c(colours[[direction]], "grey70"), alpha = 0.6),
      edges = list(col = "grey30"),
      labels = FALSE,
      quantities = list(labels = labels, fontsize = 8)
    )
    size <- extents[i] / max_extent
    diagram <- grobTree(
      diagram,
      vp = viewport(width = unit(size, "snpc"), height = unit(size, "snpc"))
    )
    legend <- legendGrob(
      c(
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
        ),
        sprintf(
          "Intersections (overlapping pairs): %s",
          format(r$n_intersections, big.mark = ",")
        )
      ),
      pch = c(21, 21, NA),
      gp = gpar(
        fill = c(colours[[direction]], "grey70", NA),
        col = "grey30",
        fontsize = 8
      )
    )
    body <- arrangeGrob(diagram, legend, heights = unit(c(2.9, 0.7), "in"))
  }
  arrangeGrob(
    heading,
    body,
    textGrob(gat_label(r), gp = gpar(fontsize = 8)),
    heights = unit(c(0.5, 3.6, 0.45), "in")
  )
}

panels <- lapply(seq_len(nrow(rows)), panel)
figure <- arrangeGrob(
  grobs = panels,
  ncol = length(directions),
  top = textGrob(title, gp = gpar(fontsize = 13, fontface = "bold")),
  padding = unit(0.3, "in")
)

pdf(
  snakemake@output[["venn"]],
  width = 4 * length(directions),
  height = 0.5 + 4.65 * length(datasets)
)
grid.draw(figure)
invisible(dev.off())
