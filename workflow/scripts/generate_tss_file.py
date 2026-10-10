"""
Writes one transcription start site (TSS) per protein coding gene (its 5'
end) to a BED file, for the TSS enrichment (ENCODE and ataqv) of genomes
without an ENCODE TSS file. This is how the ENCODE TSS files (e.g.
ENCFF493CCB for GRCh38) were made from GENCODE.

TSS enrichment scores are lower with one TSS per transcript or with all
genes, as many non-coding genes have little accessibility at their TSS.
Mitochondrial genes are left out.

The GTF file is read line by line, so that memory use stays low.
"""

import logging
import re
import sys

if "snakemake" in globals():
    gtf = snakemake.input["gtf"]
    output_bed = snakemake.output["bed"]
    log = snakemake.log[0]
else:
    # For testing: python generate_tss_file.py <annotation.gtf> <output.bed>
    if len(sys.argv) != 3:
        sys.exit("Usage: python generate_tss_file.py <annotation.gtf> <output.bed>")
    gtf, output_bed = sys.argv[1:]
    log = f"{output_bed}.log"

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

MITO = {"MT", "chrM", "M"}
ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def five_end(chrom, start, end, strand):
    """0-based BED interval of the 5' end of a 1-based GTF feature"""
    pos = start - 1 if strand == "+" else end - 1
    return chrom, pos, pos + 1, strand


genes = {}  # gene_id: TSS
logging.info(f"Reading {gtf}")
with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[2] != "gene" or fields[0] in MITO:
            continue
        values = dict(ATTRIBUTE.findall(fields[8]))
        # Ensembl: gene_biotype, GENCODE: gene_type
        if values.get("gene_biotype", values.get("gene_type")) != "protein_coding":
            continue
        genes[values["gene_id"]] = five_end(
            fields[0], int(fields[3]), int(fields[4]), fields[6]
        )
logging.info(f"{len(genes)} protein coding genes")
if not genes:
    sys.exit(f"No protein coding genes found in {gtf}")

with open(output_bed, "w") as fh:
    for gene_id, (chrom, start, end, strand) in sorted(
        genes.items(), key=lambda x: (x[1], x[0])
    ):
        fh.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")
logging.info("Done!")
