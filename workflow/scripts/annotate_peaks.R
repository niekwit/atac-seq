# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

# Load Snakemake parameters
peak_file <- snakemake@input[["peaks"]]
edb_file <- snakemake@input[["edb"]]
txdb_file <- snakemake@input[["txdb"]]
txt <- snakemake@output[["txt"]]

# Load narrowPeak file (0-based BED coordinates)
peaks <- read.delim(
  peak_file,
  header = FALSE,
  col.names = c(
    "chr",
    "start",
    "end",
    "peak_id",
    "score",
    "strand",
    "signal_value",
    "p_value",
    "q_value",
    "summit"
  ),
  colClasses = c("character", rep(NA, 9))
) %>%
  mutate(start = start + 1, strand = "*")
peaks <- makeGRangesFromDataFrame(
  peaks,
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = FALSE
)

# Annotate peaks
txdb <- AnnotationDbi::loadDb(txdb_file)
seqlevels(peaks, pruning.mode = "coarse") <- intersect(
  seqlevels(peaks),
  seqlevels(txdb)
)
peakAnno <- annotatePeak(peaks, tssRegion = c(-3000, 3000), TxDb = txdb)

# Add gene names and gene biotype to annotation
load(edb_file)
df <- as.data.frame(peakAnno) %>%
  left_join(edb, by = "geneId")

# Write annotation to file
write.table(df, file = txt, quote = FALSE, sep = "\t", row.names = FALSE)
