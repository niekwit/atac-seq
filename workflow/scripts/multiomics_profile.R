# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Profile of all data sets of a contrast over the full-length copies of the
# te: profile: subfamilies entries (integration of all data sets, TE
# module): mean signal of the control (grey) and the treatment (red) per bin,
# one row per data set (ATAC-seq, H3K9me3, CpG methylation) and one column per
# entry, copies scaled from 5' to 3' with flanks (from profile.tsv of
# multiomics_copies.py)

library(ggplot2)
library(cowplot)

control <- snakemake@params[["control"]]
treatment <- snakemake@params[["treatment"]]
body_length <- snakemake@params[["body_length"]]
flank <- snakemake@params[["flank"]]

profile <- read.delim(snakemake@input[[1]], check.names = FALSE)
profile <- profile[!is.na(profile$signal), ]
groups <- unique(profile$group)
n_copies <- tapply(profile$n_copies, profile$group, max)
profile$group <- factor(
  profile$group,
  levels = groups,
  labels = sprintf("%s (n = %d)", groups, n_copies[groups])
)
units <- unlist(snakemake@params[["units"]])
datasets <- unique(profile$dataset)
# Units on a second line of the y axis titles (strips)
profile$dataset <- factor(
  profile$dataset,
  levels = datasets,
  labels = sprintf("%s\n(%s)", datasets, units[datasets])
)
profile$condition <- factor(profile$condition, levels = c(control, treatment))

plot <- ggplot(profile, aes(position, signal, colour = condition)) +
  annotate(
    "rect",
    xmin = 0,
    xmax = body_length,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey95"
  ) +
  geom_line(linewidth = 0.7) +
  facet_grid(dataset ~ group, scales = "free_y", switch = "y") +
  scale_colour_manual(
    values = setNames(c("#999999", "#dd3b3b"), c(control, treatment)),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(-flank, 0, body_length, body_length + flank),
    labels = c(
      sprintf("-%g kb", flank / 1000),
      "5'",
      "3'",
      sprintf("+%g kb", flank / 1000)
    )
  ) +
  labs(x = NULL, y = NULL) +
  theme_cowplot(12) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size = 9),
    legend.position = "top"
  )
ggsave(
  snakemake@output[[1]],
  plot,
  width = 2.2 + 2.6 * length(groups),
  height = 1.2 + 2.2 * nlevels(profile$dataset),
  bg = "white",
  # Wide with many entries of te: profile: subfamilies
  limitsize = FALSE
)
cat("Done!\n")
