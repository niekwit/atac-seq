"""
Generates fold enrichment and p-value signal tracks (bigWig) from the MACS2
pileup and lambda bedGraphs, as in the ENCODE ATAC-seq pipeline
(encode_task_macs2_signal_track_atac.py):

- Fold enrichment: pileup / local lambda (macs2 bdgcmp -m FE)
- P-value: -log10 Poisson p-value of the pileup given the local lambda
  (macs2 bdgcmp -m ppois)

The bedGraphs are genome-wide (tens of millions of lines), so sorting and
format conversion are left to the UCSC tools and GNU sort.
"""

import gzip
import logging
import os
import subprocess

# Load Snakemake variables
treat_bdg = snakemake.input["treat"]
control_bdg = snakemake.input["control"]
tagalign = snakemake.input["ta"]
chrom_sizes = snakemake.input["chrom_sizes"]
fc_bigwig = snakemake.output["fc"]
pval_bigwig = snakemake.output["pval"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

prefix = fc_bigwig.removesuffix(".fc.signal.bigwig")


def run(command):
    """Runs a shell command (pipes allowed), with stderr to the log"""
    logging.info(f"Running: {command}")
    with open(log, "a") as log_fh:
        subprocess.run(
            f"set -o pipefail; {command}",
            shell=True,
            executable="/bin/bash",
            check=True,
            stderr=log_fh,
        )


def remove_overlaps(sorted_bedgraph, out_bedgraph):
    """
    Writes a sorted bedGraph without overlapping intervals:
    an interval is dropped if it starts before the end of the previous
    interval on the same chromosome (bedGraphToBigWig rejects overlaps).
    """
    prev_chrom, prev_end = None, 0
    kept = dropped = 0
    with open(sorted_bedgraph) as fh_in, open(out_bedgraph, "w") as fh_out:
        for line in fh_in:
            chrom, start, end = line.split("\t", 3)[:3]
            if chrom != prev_chrom or int(start) >= prev_end:
                fh_out.write(line)
                kept += 1
            else:
                dropped += 1
            prev_chrom, prev_end = chrom, int(end)
    logging.info(f"Kept {kept} intervals, removed {dropped} overlapping intervals")


def bedgraph_to_bigwig(bedgraph, bigwig):
    """Clips a MACS2 bedGraph to the chromosome sizes and converts it to bigWig"""
    clipped = f"{bedgraph}.clipped"
    sorted_bdg = f"{bedgraph}.sorted"
    no_overlaps = f"{bedgraph}.no_overlaps"

    # Clip intervals to 0-chromosome size (bedtools slop -b 0 clamps
    # coordinates to the chromosome; bedClip drops what is still outside)
    run(
        f"bedtools slop -i {bedgraph} -g {chrom_sizes} -b 0 | "
        f"bedClip stdin {chrom_sizes} {clipped}"
    )

    # Sort by chromosome and start, as bedGraphToBigWig requires
    run(f"LC_COLLATE=C sort -k1,1 -k2,2n {clipped} > {sorted_bdg}")
    os.remove(clipped)

    remove_overlaps(sorted_bdg, no_overlaps)
    os.remove(sorted_bdg)

    run(f"bedGraphToBigWig {no_overlaps} {chrom_sizes} {bigwig}")
    os.remove(no_overlaps)
    os.remove(bedgraph)


# Fold enrichment over local background
# -----------------------------------------------------
run(f"macs2 bdgcmp -t {treat_bdg} -c {control_bdg} --o-prefix {prefix} -m FE")
bedgraph_to_bigwig(f"{prefix}_FE.bdg", fc_bigwig)

# -log10(p-value)
# -----------------------------------------------------
# The pileup is per million reads (--SPMR), so -S scales it back to read
# counts: sval is the number of reads in millions
with gzip.open(tagalign, "rt") as fh:
    sval = sum(1 for _ in fh) / 1e6
logging.info(f"Reads (millions) for p-value scaling: {sval}")

run(
    f"macs2 bdgcmp -t {treat_bdg} -c {control_bdg} --o-prefix {prefix} "
    f"-m ppois -S {sval}"
)
bedgraph_to_bigwig(f"{prefix}_ppois.bdg", pval_bigwig)
logging.info("Done!")
