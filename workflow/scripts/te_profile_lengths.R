# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# QC of the full-length copies of the te: profile: subfamilies entries
# (resources/te/full_length_copies.bed): histogram of their genomic lengths,
# one panel per entry, with the median length (dashed) and the length that
# computeMatrix scale-regions scales each copy to (te: profile: body_length,
# dotted), explained in the legend. Copies of an entry are followed by a line
# "#<entry> (n=...)".

library(ggplot2)
library(cowplot)

body_length <- snakemake@params[["body_length"]]

# Copies with the label of their entry
lines <- readLines(snakemake@input[["bed"]])
is_label <- startsWith(lines, "#")
labels <- sub("^#", "", lines[is_label])
# Entry of each line: the next label line (the label follows its copies)
entry <- labels[length(labels) - rev(cumsum(rev(is_label))) + 1]
copies <- read.delim(
  text = lines[!is_label],
  header = FALSE,
  colClasses = c(
    "character",
    "integer",
    "integer",
    "character",
    "numeric",
    "character"
  )
)
copies$entry <- factor(entry[!is_label], levels = labels)
copies$length <- (copies$V3 - copies$V2) / 1000

medians <- aggregate(length ~ entry, data = copies, FUN = median)
for (i in seq_len(nrow(medians))) {
  message(sprintf(
    "%s: median length %.2f kb",
    medians$entry[i],
    medians$length[i]
  ))
}

# Vertical lines with a legend: median per entry and the scaled length
median_label <- "Median length"
body_label <- sprintf(
  "Scaled length in the profile (%g kb)",
  body_length / 1000
)
lines_df <- rbind(
  data.frame(entry = medians$entry, x = medians$length, line = median_label),
  data.frame(
    entry = factor(labels, levels = labels),
    x = body_length / 1000,
    line = body_label
  )
)
lines_df$line <- factor(lines_df$line, levels = c(median_label, body_label))

p <- ggplot(copies, aes(x = length)) +
  geom_histogram(
    bins = 50,
    fill = "#cccccc",
    colour = "black",
    linewidth = 0.2
  ) +
  geom_vline(
    data = lines_df,
    aes(xintercept = x, colour = line, linetype = line)
  ) +
  scale_colour_manual(
    values = setNames(c("#dd3b3b", "black"), c(median_label, body_label))
  ) +
  scale_linetype_manual(
    values = setNames(c("dashed", "dotted"), c(median_label, body_label))
  ) +
  facet_wrap(~entry, scales = "free", ncol = 3) +
  labs(x = "Copy length (kb)", y = "Copies", colour = NULL, linetype = NULL) +
  theme_cowplot(12) +
  theme(strip.background = element_blank(), legend.position = "bottom")

n_rows <- ceiling(length(labels) / 3)
ggsave(
  snakemake@output[[1]],
  p,
  width = 10,
  height = 1 + 2.5 * n_rows,
  bg = "white"
)
