"""
Fingerprint plot and Jensen-Shannon distance (JSD) of the replicates of a
condition with deepTools plotFingerprint, as in the ENCODE ATAC-seq pipeline
(encode_task_jsd.py):

1. Reads overlapping blacklisted regions are removed from the filtered,
   deduplicated BAM files
2. plotFingerprint on all replicates together (reads with MAPQ >= 30)

ATAC-seq has no control, so the metrics compare each replicate with a
synthetic, uniformly covered sample (e.g. "Synthetic JS Distance").
"""

import logging
import os
import subprocess

# Load Snakemake variables
bams = list(snakemake.input["bams"])
blacklist = snakemake.input["blacklist"]
samples = snakemake.params["samples"]
mapq = snakemake.params["mapq"]
plot_png = snakemake.output["plot"]
metrics = snakemake.output["metrics"]
threads = str(snakemake.threads)
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def run(command, **kwargs):
    logging.info(f"Running: {' '.join(command)}")
    with open(log, "a") as log_fh:
        subprocess.run(command, check=True, stderr=log_fh, **kwargs)


# Remove reads in blacklisted regions
# -----------------------------------------------------
tmp_dir = os.path.dirname(metrics)
filtered_bams = []
try:
    for sample, bam in zip(samples, bams):
        filtered = os.path.join(tmp_dir, f"{sample}.jsd_tmp.bfilt.bam")
        with open(filtered, "wb") as fh:
            run(
                ["bedtools", "intersect", "-nonamecheck", "-v", "-abam", bam]
                + ["-b", blacklist],
                stdout=fh,
            )
        run(["samtools", "index", "-@", threads, filtered])
        filtered_bams.append(filtered)

    # Fingerprint and JSD
    # -----------------------------------------------------
    run(
        ["plotFingerprint", "-b"]
        + filtered_bams
        + ["--labels"]
        + samples
        + [
            "--outQualityMetrics", metrics,
            "--minMappingQuality", str(mapq),
            "-T", "Fingerprints of different samples",
            "--numberOfProcessors", threads,
            "--plotFile", plot_png,
        ]  # fmt: skip
    )
finally:
    for bam in filtered_bams:
        for f in [bam, f"{bam}.bai"]:
            if os.path.exists(f):
                os.remove(f)
logging.info("Done!")
