# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Repeat subfamilies in the differential H3K9me3 regions of a contrast
# against their differential accessibility (TE module).
#
# - subfamilies.tsv, correlation.tsv, scatter.pdf: per subfamily with at
#   least min_copies copies and tested by DESeq2, the percentage of copies in
#   regions that lose (and gain) H3K9me3 against the shrunken log2 fold change
#   of the insertions in all copies; Spearman correlation per H3K9me3 change.
#   The dashed line is the percentage of the genome in the regions (expected
#   for copies placed at random).
# - full_length.tsv, full_length.pdf: full-length copies of the te: profile:
#   subfamilies entries by H3K9me3 change and accessibility (TE peaks), with
#   Fisher's exact test of the copies in the H3K9me3 regions against the
#   copies with a more accessible TE peak

library(ggplot2)
library(ggrepel)
library(cowplot)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
min_copies <- snakemake@params[["min_copies"]]
changes <- intersect(c("lost", "gained"), snakemake@params[["changes"]])
n_labels <- 15

message_plot <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message) +
    theme_void(base_size = 12)
}
stars <- function(p) {
  ifelse(
    is.na(p) | p >= 0.05,
    "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
  )
}

# Subfamilies: copies in the H3K9me3 regions and accessibility
# -----------------------------------------------------
copies <- read.delim(
  snakemake@input[["copies"]],
  colClasses = c(subfamily = "character")
)
atac <- read.delim(
  snakemake@input[["atac"]],
  colClasses = c(subfamily = "character")
)
atac <- atac[!is.na(atac$padj), ]
subfamilies <- merge(
  atac[, c(
    "subfamily",
    "family",
    "class",
    "baseMean",
    "log2FoldChange",
    "log2FoldChange_shrunken",
    "padj",
    "change"
  )],
  copies,
  by = "subfamily"
)
subfamilies <- subfamilies[subfamilies$n_copies >= min_copies, ]
for (change in c("lost", "gained")) {
  subfamilies[[paste0("percent_copies_", change)]] <- 100 *
    subfamilies[[paste0("n_copies_", change)]] /
    subfamilies$n_copies
}
cat(sprintf(
  "%d subfamilies tested by DESeq2 with at least %d copies\n",
  nrow(subfamilies),
  min_copies
))
subfamilies <- subfamilies[
  order(
    -subfamilies[[paste0(
      "percent_copies_",
      changes[1]
    )]]
  ),
]
write.table(
  subfamilies,
  snakemake@output[["subfamilies"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

correlation <- do.call(
  rbind,
  lapply(changes, function(change) {
    x <- subfamilies[[paste0("percent_copies_", change)]]
    y <- subfamilies$log2FoldChange_shrunken
    test <- if (length(x) >= 3) {
      suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    } else {
      list(estimate = NA, p.value = NA)
    }
    data.frame(
      h3k9me3_change = change,
      n_subfamilies = length(x),
      spearman_rho = unname(test$estimate),
      p = test$p.value
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

change_labels <- c(
  "lost" = paste("H3K9me3 lost in", treatment),
  "gained" = paste("H3K9me3 gained in", treatment)
)
if (nrow(subfamilies) == 0) {
  scatter <- message_plot(paste(
    "No subfamilies tested by DESeq2 with at least",
    min_copies,
    "copies"
  ))
} else {
  plot_data <- do.call(
    rbind,
    lapply(changes, function(change) {
      percent <- subfamilies[[paste0("percent_copies_", change)]]
      # Label the differential subfamilies with most copies in the regions
      # and the largest change (sqrt(percent) x |log2 fold change|)
      score <- sqrt(percent) * abs(subfamilies$log2FoldChange_shrunken)
      ranked <- order(subfamilies$change == "unchanged", -score)
      labelled <- subfamilies$subfamily[head(ranked, n_labels)]
      data.frame(
        h3k9me3_change = change_labels[[change]],
        subfamily = subfamilies$subfamily,
        class = subfamilies$class,
        percent = percent,
        log2FoldChange_shrunken = subfamilies$log2FoldChange_shrunken,
        atac_change = subfamilies$change,
        label = ifelse(
          subfamilies$subfamily %in% labelled,
          subfamilies$subfamily,
          ""
        )
      )
    })
  )
  plot_data$h3k9me3_change <- factor(
    plot_data$h3k9me3_change,
    levels = change_labels[changes]
  )
  plot_data$significant <- factor(
    ifelse(plot_data$atac_change == "unchanged", "no", "yes"),
    levels = c("no", "yes")
  )
  plot_data <- plot_data[order(plot_data$significant), ]
  expected <- data.frame(
    h3k9me3_change = factor(
      change_labels[changes],
      levels = levels(
        plot_data$h3k9me3_change
      )
    ),
    percent = sapply(changes, function(change) {
      100 * copies[[paste0("genome_fraction_", change)]][1]
    })
  )
  caption <- paste(
    sprintf(
      "%s: Spearman rho = %.2f (p = %.2g, n = %d)",
      change_labels[correlation$h3k9me3_change],
      correlation$spearman_rho,
      correlation$p,
      correlation$n_subfamilies
    ),
    collapse = "\n"
  )
  scatter <- ggplot(
    plot_data,
    aes(percent, log2FoldChange_shrunken, label = label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(
      data = expected,
      aes(xintercept = percent),
      linetype = "dashed",
      colour = "grey50"
    ) +
    geom_point(
      aes(fill = class, alpha = significant),
      shape = 21,
      size = 2.5,
      colour = "black",
      stroke = 0.3
    ) +
    geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
    facet_wrap(~h3k9me3_change, scales = "free_x") +
    scale_x_continuous(trans = "sqrt") +
    scale_fill_brewer(palette = "Set1") +
    scale_alpha_manual(values = c("no" = 0.35, "yes" = 1), drop = FALSE) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      alpha = guide_legend(order = 2, override.aes = list(fill = "grey50"))
    ) +
    labs(
      x = "Copies in the H3K9me3 regions (%, square-root scale)",
      y = paste(
        "ATAC-seq shrunken log2 fold change (all copies),",
        treatment,
        "vs",
        control
      ),
      fill = "Class",
      alpha = "Differentially\naccessible",
      caption = paste0(
        "Repeat subfamilies with >= ",
        min_copies,
        " copies; copies with >= 50% in a region; ",
        "dashed line: % of the genome in the regions\n",
        caption
      )
    ) +
    theme_cowplot(12) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
}
ggsave(
  snakemake@output[["scatter"]],
  scatter,
  width = 3 + 5 * max(1, length(changes)),
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# Full-length copies: H3K9me3 change and accessibility
# -----------------------------------------------------
full_length <- read.delim(
  snakemake@input[["full_length"]],
  colClasses = c(chr = "character", group = "character")
)
groups <- unique(full_length$group)
more <- full_length$accessibility == "more accessible"
summary <- do.call(
  rbind,
  lapply(groups, function(group) {
    do.call(
      rbind,
      lapply(changes, function(change) {
        keep <- full_length$group == group
        in_regions <- full_length$h3k9me3_change[keep] == change
        open <- more[keep]
        table_2x2 <- table(
          factor(in_regions, levels = c(TRUE, FALSE)),
          factor(open, levels = c(TRUE, FALSE))
        )
        test <- if (
          all(rowSums(table_2x2) > 0) && all(colSums(table_2x2) > 0)
        ) {
          fisher.test(table_2x2)
        } else {
          list(estimate = NA, conf.int = c(NA, NA), p.value = NA)
        }
        data.frame(
          group = group,
          h3k9me3_change = change,
          n_copies = sum(keep),
          n_in_regions = sum(in_regions),
          n_more_accessible = sum(open),
          n_both = sum(in_regions & open),
          odds_ratio = unname(test$estimate),
          ci_low = test$conf.int[1],
          ci_high = test$conf.int[2],
          p = test$p.value
        )
      })
    )
  })
)
if (!is.null(summary)) {
  summary$padj <- p.adjust(summary$p, method = "BH")
  print(summary)
}
write.table(
  if (is.null(summary)) {
    data.frame(
      group = character(),
      h3k9me3_change = character(),
      n_copies = integer()
    )
  } else {
    summary
  },
  snakemake@output[["full_length_tests"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

if (nrow(full_length) == 0) {
  full_length_plot <- message_plot(
    "No full-length copies (te: profile: subfamilies)"
  )
} else {
  # Main change of H3K9me3 (lost if present)
  change <- changes[1]
  in_regions <- full_length$h3k9me3_change == change
  category_levels <- c(
    paste0("H3K9me3 ", change, ", more accessible"),
    paste0("H3K9me3 ", change, " only"),
    "more accessible only",
    "neither"
  )
  full_length$category <- factor(
    ifelse(
      in_regions & more,
      category_levels[1],
      ifelse(
        in_regions,
        category_levels[2],
        ifelse(more, category_levels[3], category_levels[4])
      )
    ),
    levels = rev(category_levels)
  )
  category_colours <- setNames(
    c("#984EA3", "#377EB8", "#dd3b3b", "#cccccc"),
    category_levels
  )
  full_length$group <- factor(full_length$group, levels = groups)
  tests <- summary[summary$h3k9me3_change == change, ]
  annotations <- data.frame(
    group = factor(tests$group, levels = groups),
    label = paste0("n = ", tests$n_copies, "\n", stars(tests$padj))
  )
  full_length_plot <- ggplot(full_length, aes(group, fill = category)) +
    geom_bar(
      position = "fill",
      colour = "black",
      linewidth = 0.3,
      width = 0.7
    ) +
    geom_text(
      data = annotations,
      aes(group, 1, label = label),
      inherit.aes = FALSE,
      vjust = -0.2,
      size = 3.5,
      lineheight = 0.9
    ) +
    # All categories in the legend, also those without copies
    scale_fill_manual(values = category_colours, drop = FALSE) +
    scale_y_continuous(
      labels = function(x) paste0(100 * x, "%"),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = NULL,
      y = "Full-length copies",
      fill = paste(treatment, "vs", control),
      caption = paste0(
        "More accessible: overlapped by a more accessible TE peak\n",
        "Fisher's exact test, H3K9me3 ",
        change,
        " vs more accessible; BH-adjusted p: * < 0.05, ** < 0.01, *** < 0.001"
      )
    ) +
    theme_cowplot(12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0, size = 10)
    )
}
ggsave(
  snakemake@output[["full_length"]],
  full_length_plot,
  width = 4 + 0.8 * max(1, length(groups)),
  height = 6,
  dpi = 300,
  bg = "white"
)
