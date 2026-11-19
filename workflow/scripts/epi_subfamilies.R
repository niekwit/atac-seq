# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Repeat subfamilies in the differential regions of an epigenomic mark
# (H3K9me3, DNA methylation) of a contrast against their differential
# accessibility (TE module).
#
# - subfamilies.tsv, correlation.tsv, scatter.pdf: per subfamily with at
#   least min_copies copies and tested by DESeq2, the percentage of copies in
#   regions where the mark decreases (and increases) against the shrunken log2
#   fold change of the insertions in all copies; Spearman correlation per
#   change of the mark.
#   The dashed line is the percentage of the genome in the regions, as a
#   reference (not the expected percentage of copies under random placement,
#   which also depends on the lengths of the copies and regions).
# - full_length.tsv, full_length.pdf: full-length copies of the te: profile:
#   subfamilies entries by change of the mark and accessibility (TE peaks),
#   with Fisher's exact test of the copies in the regions against the copies
#   with a more accessible TE peak

library(ggplot2)
library(ggrepel)
library(cowplot)

treatment <- snakemake@params[["treatment"]]
control <- snakemake@params[["control"]]
min_copies <- snakemake@params[["min_copies"]]
# Mark (h3k9me3, methylation): prefix of its columns; change codes of a
# decrease and an increase in the treatment; plot labels
prefix <- snakemake@params[["prefix"]]
down <- snakemake@params[["down"]]
up <- snakemake@params[["up"]]
labels <- snakemake@params[["labels"]]
changes <- intersect(c(down, up), snakemake@params[["changes"]])
n_labels <- 15

# Writes a table with the column names of the mark (e.g. h3k9me3_change)
write_tsv <- function(data, path) {
  names(data) <- sub("^mark_", paste0(prefix, "_"), names(data))
  write.table(data, path, quote = FALSE, sep = "\t", row.names = FALSE)
}

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

# Subfamilies: copies in the regions and accessibility
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
for (change in c(down, up)) {
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
write_tsv(
  subfamilies,
  snakemake@output[["subfamilies"]]
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
      mark_change = change,
      n_subfamilies = length(x),
      spearman_rho = unname(test$estimate),
      p = test$p.value
    )
  })
)
print(correlation)
write_tsv(
  correlation,
  snakemake@output[["correlation"]]
)

change_labels <- setNames(c(labels[["down"]], labels[["up"]]), c(down, up))
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
        mark_change = change_labels[[change]],
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
  plot_data$mark_change <- factor(
    plot_data$mark_change,
    levels = change_labels[changes]
  )
  plot_data$significant <- factor(
    ifelse(plot_data$atac_change == "unchanged", "no", "yes"),
    levels = c("no", "yes")
  )
  plot_data <- plot_data[order(plot_data$significant), ]
  expected <- data.frame(
    mark_change = factor(
      change_labels[changes],
      levels = levels(
        plot_data$mark_change
      )
    ),
    percent = sapply(changes, function(change) {
      100 * copies[[paste0("genome_fraction_", change)]][1]
    })
  )
  caption <- paste(
    sprintf(
      "%s: Spearman rho = %.2f (p = %.2g, n = %d)",
      change_labels[correlation$mark_change],
      correlation$spearman_rho,
      correlation$p,
      correlation$n_subfamilies
    ),
    collapse = "\n"
  )
  classes <- sort(unique(plot_data$class))
  class_colours <- setNames(
    if (length(classes) <= 9) {
      RColorBrewer::brewer.pal(9, "Set1")[seq_along(classes)]
    } else {
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(
        length(classes)
      )
    },
    classes
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
    facet_wrap(~mark_change, scales = "free_x") +
    scale_x_continuous(trans = "sqrt") +
    # Set1, extended for more than 9 classes
    scale_fill_manual(values = class_colours) +
    scale_alpha_manual(values = c("no" = 0.35, "yes" = 1), drop = FALSE) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      alpha = guide_legend(order = 2, override.aes = list(fill = "grey50"))
    ) +
    labs(
      x = paste0(
        "Copies in the ",
        labels[["regions"]],
        " (%, square-root scale)"
      ),
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

# Full-length copies: change of the mark and accessibility
# -----------------------------------------------------
full_length <- read.delim(
  snakemake@input[["full_length"]],
  colClasses = c(chr = "character", group = "character")
)
names(full_length) <- sub(paste0("^", prefix, "_"), "mark_", names(full_length))
groups <- unique(full_length$group)
more <- full_length$accessibility == "more accessible"
summary <- do.call(
  rbind,
  lapply(groups, function(group) {
    do.call(
      rbind,
      lapply(changes, function(change) {
        keep <- full_length$group == group
        in_regions <- full_length$mark_change[keep] == change
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
          mark_change = change,
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
write_tsv(
  if (is.null(summary)) {
    data.frame(
      group = character(),
      mark_change = character(),
      n_copies = integer()
    )
  } else {
    summary
  },
  snakemake@output[["full_length_tests"]]
)

if (nrow(full_length) == 0) {
  full_length_plot <- message_plot(
    "No full-length copies (te: profile: subfamilies)"
  )
} else {
  # Main change of the mark (the decrease if present)
  change <- changes[1]
  change_short <- labels[[if (change == down) "down_short" else "up_short"]]
  in_regions <- full_length$mark_change == change
  category_levels <- c(
    paste0(change_short, ", more accessible"),
    paste0(change_short, " only"),
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
  tests <- summary[summary$mark_change == change, ]
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
      # Fixed breaks: the space above the bars for the labels would otherwise
      # give breaks at 0, 30, 60 and 90%
      breaks = seq(0, 1, 0.25),
      labels = function(x) paste0(100 * x, "%"),
      expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = NULL,
      y = "Full-length copies",
      fill = paste(treatment, "vs", control),
      caption = paste0(
        "More accessible: overlapped by a more accessible TE peak\n",
        "Fisher's exact test, ",
        tolower(change_short),
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
  # Many te: profile: subfamilies entries make a wide plot
  limitsize = FALSE,
  height = 6,
  dpi = 300,
  bg = "white"
)
