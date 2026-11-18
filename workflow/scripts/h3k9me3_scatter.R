# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Correlation of the ATAC-seq log2 fold changes of the peaks with the H3K9me3
# log2 fold changes (DiffBind) of the region each peak overlaps most
# (h3k9me3_annotate_peaks), for the peaks in a DiffBind site. A site with
# several peaks contributes one point per peak.
#
# - correlation.tsv: Spearman correlation of the peak-site pairs
#   (descriptive: peaks in the same site are not independent, so no p-value)
#   and of the sites, with the mean shrunken log2 fold change of their peaks
#   (one value per site, with p-value), over all pairs or sites and over those
#   with a significant change in both (sites: significant H3K9me3 change and
#   at least one differential peak), with the number of pairs per quadrant
#   (significant in both)
# - scatter.pdf: the pairs, coloured when significant in both

library(ggplot2)
library(cowplot)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
peak_set <- snakemake@params[["peak_set"]]

peaks <- read.delim(
  snakemake@input[["peaks"]],
  colClasses = c(chr = "character")
)
pairs <- peaks[
  peaks$h3k9me3_change != "outside" &
    !is.na(peaks$h3k9me3_log2FoldChange),
]
pairs$significant <- pairs$change != "unchanged" &
  pairs$h3k9me3_change != "unchanged"
cat(sprintf(
  "%d peaks in a DiffBind site, %d significant in both\n",
  nrow(pairs),
  sum(pairs$significant)
))

spearman <- function(x, y) {
  if (length(x) < 3) {
    return(c(rho = NA, p = NA, n = length(x)))
  }
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p = test$p.value, n = length(x))
}
quadrant <- function(atac, h3k9me3) {
  sum(
    pairs$significant &
      pairs$change == atac &
      pairs$h3k9me3_change == h3k9me3
  )
}
# One value per site: the mean shrunken log2 fold change of its peaks
sites <- if (nrow(pairs)) {
  merge(
    aggregate(
      cbind(log2FoldChange_shrunken, significant) ~
        h3k9me3_region +
          h3k9me3_log2FoldChange,
      data = pairs,
      FUN = mean
    ),
    aggregate(significant ~ h3k9me3_region, data = pairs, FUN = any),
    by = "h3k9me3_region",
    suffixes = c("_fraction", "")
  )
} else {
  data.frame(
    h3k9me3_log2FoldChange = numeric(),
    log2FoldChange_shrunken = numeric(),
    significant = logical()
  )
}
correlation <- do.call(
  rbind,
  lapply(c("peaks", "sites"), function(level) {
    data <- if (level == "peaks") pairs else sites
    do.call(
      rbind,
      lapply(c("all", "significant in both"), function(set) {
        keep <- if (set == "all") rep(TRUE, nrow(data)) else data$significant
        result <- spearman(
          data$h3k9me3_log2FoldChange[keep],
          data$log2FoldChange_shrunken[keep]
        )
        data.frame(
          level = level,
          pairs = set,
          n = result[["n"]],
          spearman_rho = result[["rho"]],
          # Peaks of a site are not independent: no p-value
          p = if (level == "peaks") NA else result[["p"]],
          n_lost_more_accessible = quadrant("more accessible", "lost"),
          n_lost_less_accessible = quadrant("less accessible", "lost"),
          n_gained_more_accessible = quadrant("more accessible", "gained"),
          n_gained_less_accessible = quadrant("less accessible", "gained")
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

if (nrow(pairs) == 0) {
  scatter <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No peaks in DiffBind sites") +
    theme_void(base_size = 12)
} else {
  pairs$colour <- ifelse(pairs$significant, pairs$change, "other")
  pairs$colour <- factor(
    pairs$colour,
    levels = c("other", "more accessible", "less accessible")
  )
  pairs <- pairs[order(pairs$colour), ]
  # Pairs significant in both per quadrant, in the corners
  corners <- data.frame(
    x = c(-Inf, -Inf, Inf, Inf),
    y = c(Inf, -Inf, Inf, -Inf),
    hjust = c(-0.2, -0.2, 1.2, 1.2),
    vjust = c(1.5, -0.5, 1.5, -0.5),
    label = c(
      quadrant("more accessible", "lost"),
      quadrant("less accessible", "lost"),
      quadrant("more accessible", "gained"),
      quadrant("less accessible", "gained")
    )
  )
  by_peak <- correlation[
    correlation$level == "peaks" & correlation$pairs == "all",
  ]
  by_site <- correlation[
    correlation$level == "sites" & correlation$pairs == "all",
  ]
  scatter <- ggplot(
    pairs,
    aes(h3k9me3_log2FoldChange, log2FoldChange_shrunken)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      aes(fill = colour, alpha = colour),
      shape = 21,
      size = 2,
      colour = "black",
      stroke = 0.2
    ) +
    geom_text(
      data = corners,
      aes(x, y, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_fill_manual(
      values = c(
        "other" = "#cccccc",
        "more accessible" = "#dd3b3b",
        "less accessible" = "#377EB8"
      ),
      drop = TRUE
    ) +
    scale_alpha_manual(
      values = c(
        "other" = 0.4,
        "more accessible" = 1,
        "less accessible" = 1
      ),
      drop = TRUE
    ) +
    labs(
      x = paste("H3K9me3 log2 fold change,", treatment, "vs", control),
      y = paste(
        "ATAC-seq shrunken log2 fold change,",
        treatment,
        "vs",
        control
      ),
      fill = "Significant in both",
      alpha = "Significant in both",
      caption = sprintf(
        paste0(
          "%s in DiffBind sites; corners: pairs significant in both\n",
          "Spearman rho = %.2f (%d peaks); per site (mean of its peaks): ",
          "rho = %.2f (p = %.2g, n = %d)"
        ),
        peak_set,
        by_peak$spearman_rho,
        by_peak$n,
        by_site$spearman_rho,
        by_site$p,
        by_site$n
      )
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
}
ggsave(
  snakemake@output[["scatter"]],
  scatter,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
