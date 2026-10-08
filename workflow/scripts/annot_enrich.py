"""
Fraction of reads in annotated regions, as in the ENCODE ATAC-seq pipeline
(encode_task_annot_enrich.py): the fraction of Tn5 shifted reads (tagAlign)
that overlap universal DNase I hypersensitive sites, blacklisted regions,
promoters and enhancers. The ENCODE region files are only available for
GRCh38 and mm10; for other genomes only the blacklist fraction is reported.
"""

import gzip
import logging
import subprocess

# Load Snakemake variables
tagalign = snakemake.input["ta"]
regions = {
    "universal_DHS": snakemake.input.get("dnase"),
    "blacklist": snakemake.input["blacklist"],
    "promoter": snakemake.input.get("prom"),
    "enhancer": snakemake.input.get("enh"),
}
qc_out = snakemake.output[0]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def reads_in_regions(reads, regions_bed):
    """Number of reads overlapping the (sorted and merged) regions"""
    command = (
        f"set -o pipefail; bedtools sort -i {regions_bed} | "
        "bedtools merge -i stdin | "
        f"bedtools intersect -u -nonamecheck -a {reads} -b stdin | wc -l"
    )
    logging.info(f"Running: {command}")
    with open(log, "a") as log_fh:
        result = subprocess.run(
            command,
            shell=True,
            executable="/bin/bash",
            check=True,
            stdout=subprocess.PIPE,
            stderr=log_fh,
            text=True,
        )
    return int(result.stdout.strip())


with gzip.open(tagalign, "rt") as fh:
    total = sum(1 for _ in fh)
logging.info(f"{total} reads in {tagalign}")

with open(qc_out, "w") as fh:
    fh.write("metric\treads\tfraction\n")
    for name, regions_bed in regions.items():
        if not regions_bed:
            continue
        n = reads_in_regions(tagalign, regions_bed)
        fh.write(f"fraction_of_reads_in_{name}_regions\t{n}\t{n / total}\n")
        logging.info(f"{name}: {n} reads ({n / total:.4f})")
logging.info("Done!")
