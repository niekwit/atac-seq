import sys
import logging
import time
import pyranges1 as pr

"""
Note to self about 1-based GTF vs 0-based BED:

Some genomic coordinate formats (GFF, GTF) use a 1-based, start-and-end-included format. Pyranges takes care of converting between these conventions when loading and writing files in these formats.

Source:
https://pyranges1.readthedocs.io/en/latest/tutorial.html

"""

if "snakemake" not in globals():
    # For testing pruposes
    # Get command line arguments
    args = sys.argv[1:]
    if len(args) != 2:
        print("Usage: python generate_tss_file.py <annotations.gtf> <output.bed>")
        sys.exit(1)
    gtf = args[0]
    output_bed = args[1]
    date_time = time.strftime("%Y-%m-%d_%H:%M:%S")
    log = f"generate_tss_file_{date_time}.log"
else:
    gtf = snakemake.input["gtf"]
    output_bed = snakemake.output["bed"]
    log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# Load the GTF
logging.info(f"Loading GTF file {gtf}")
gr = pr.read_gtf(gtf)

# Filter for protein coding genes: one TSS per protein coding gene, as in the
# ENCODE TSS files (e.g. ENCFF493CCB for GRCh38). TSS enrichment scores are
# lower with one TSS per transcript or with all genes, as many non-coding
# genes have little accessibility at their TSS.
logging.info("Filtering for protein coding genes...")
genes = gr[(gr.Feature == "gene") & (gr.gene_biotype == "protein_coding")]

# Remove mitochondrial genes
logging.info("Removing mitochondrial genes...")
genes = genes[genes.Chromosome != "MT"]

# Calculate TSS
# .five_end() handles the strand-specific logic
logging.info("Calculating TSS...")
tss = genes.five_end()

logging.info("Using Gene ID for the Name column...")
tss["Name"] = tss.gene_id

# Add a dummy score
tss["Score"] = 0

# Select only the columns needed for BED
tss = tss[["Chromosome", "Start", "End", "Name", "Score", "Strand"]]

# Export to BED
logging.info(f"Exporting TSS to BED file: {output_bed}")
output_bed = output_bed
tss.to_bed(output_bed)
logging.info("Done!")
