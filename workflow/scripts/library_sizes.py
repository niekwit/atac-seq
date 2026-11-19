"""
Read pairs (read 1 of each pair, samtools view -c -f 64) of the filtered
BAM files of all samples, for subsample_reads: min (subsampling of the
tagAlign files to the smallest library).
"""

import logging
import subprocess

samples = snakemake.params["samples"]
bams = snakemake.input["bams"]
threads = str(snakemake.threads)

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(snakemake.log[0])],
    force=True,
)

rows = ["sample\tread_pairs"]
for sample, bam in zip(samples, bams):
    command = ["samtools", "view", "-c", "-f", "64", "-@", threads, bam]
    pairs = int(
        subprocess.run(command, check=True, capture_output=True, text=True).stdout
    )
    logging.info(f"{sample}: {pairs} read pairs ({bam})")
    rows.append(f"{sample}\t{pairs}")
empty = [row.split("\t")[0] for row in rows[1:] if row.endswith("\t0")]
if empty:
    # Stops before any tagAlign file is made: a minimum of 0 cannot be used
    raise ValueError(
        f"subsample_reads: min: no read pairs in {', '.join(empty)}; remove these "
        "samples or set subsample_reads to a number"
    )
with open(snakemake.output[0], "w") as fh:
    fh.write("\n".join(rows) + "\n")
