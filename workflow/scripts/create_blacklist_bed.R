# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


library(rtracklayer)
library(AnnotationHub)
# https://doi.org/10.1093/bioinformatics/btad198

# Get genome
genome_assembly <- snakemake@params[["genome"]]

print(paste("Fetching blacklist regions for:", genome_assembly))

# Get blacklist regions from AnnotationHub
ah <- AnnotationHub()
query <- query(ah, c(genome_assembly, "blacklist"))
blacklist <- query[[1]]

# Convert to Ensemble style chromosome names
seqlevelsStyle(blacklist) <- "Ensembl"

# Fix all 0 coordinates
# export function annoyingly 1-based 0 to -1 0-based (BED format)
# Convert this to 1 to avoid negative coordinates in the BED file
start(blacklist) <- pmax(start(blacklist), 1)

# Export to BED file
output_file <- snakemake@output[["bed"]]
export(blacklist, output_file, format = "BED")
