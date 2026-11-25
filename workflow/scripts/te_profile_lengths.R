# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# QC of the full-length copies of the te: profile: subfamilies entries: length
# distribution of the full-length copies (red) and of all copies (grey) of each
# entry, one panel per entry, with the median length of both (dashed). The
# histograms are scaled to the same area (density), as an entry has many more
# copies than full-length copies.

library(ggplot2)
library(cowplot)

copies <- read.delim(
  snakemake@input[["copies"]],
  colClasses = c(entry = "character", chrom = "character")
)
copies$length <- copies$length / 1000
copies$full_length <- as.logical(copies$full_length)

# Panel labels with the number of full-length copies out of all copies
entries <- unique(copies$entry)
n_all <- table(factor(copies$entry, levels = entries))
n_full <- tapply(
  copies$full_length,
  factor(copies$entry, levels = entries),
  sum
)
panel_labels <- sprintf(
  "%s (%s of %s full-length)",
  entries,
  format(n_full, big.mark = ",", trim = TRUE),
  format(as.vector(n_all), big.mark = ",", trim = TRUE)
)
names(panel_labels) <- entries

# All copies and the full-length copies as two sets
all_label <- "All copies"
full_label <- "Full-length copies"
sets <- rbind(
  data.frame(entry = copies$entry, length = copies$length, set = all_label),
  data.frame(
    entry = copies$entry[copies$full_length],
    length = copies$length[copies$full_length],
    set = full_label
  )
)
sets$entry <- factor(panel_labels[sets$entry], levels = panel_labels)
sets$set <- factor(sets$set, levels = c(all_label, full_label))

medians <- aggregate(length ~ entry + set, data = sets, FUN = median)
for (i in seq_len(nrow(medians))) {
  message(sprintf(
    "%s, %s: median length %.2f kb",
    medians$entry[i],
    medians$set[i],
    medians$length[i]
  ))
}

colours <- setNames(c("#cccccc", "#dd3b3b"), c(all_label, full_label))
p <- ggplot(sets, aes(x = length, fill = set)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    position = "identity",
    alpha = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    data = medians,
    aes(xintercept = length, colour = set),
    linetype = "dashed",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = colours) +
  # Darker grey than the bars for the median of all copies
  scale_colour_manual(
    values = setNames(c("grey40", "#dd3b3b"), c(all_label, full_label))
  ) +
  facet_wrap(~entry, scales = "free", ncol = 3) +
  labs(
    x = "Copy length (kb)",
    y = "Density",
    fill = NULL,
    caption = "Dashed lines: median length"
  ) +
  theme_cowplot(12) +
  theme(strip.background = element_blank(), legend.position = "bottom")

n_rows <- ceiling(length(entries) / 3)
ggsave(
  snakemake@output[[1]],
  p,
  width = 11,
  height = 1 + 2.5 * n_rows,
  bg = "white",
  # Height grows with the number of entries, beyond ggsave's 50 inch limit
  # from about 58 entries
  limitsize = FALSE
)
