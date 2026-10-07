"""
Splits a paired-end tagAlign file into two pseudoreplicates by randomly
assigning read pairs, as in the ENCODE ATAC-seq pipeline (encode_task_spr.py).

Why: IDR and naive overlap need two peak sets to compare. Next to the
comparisons of the true biological replicates, ENCODE compares two
pseudoreplicates of each replicate (self-consistency) and of the pooled
replicates. Each pseudoreplicate holds a random half of the reads, so the
peaks found in both halves are those that do not depend on sampling noise
alone. The number of reproducible peaks between pseudoreplicates is compared
with that between true replicates (rescue and self-consistency ratios) and
gives the optimal and conservative peak sets.

What:
1. Count the read pairs. Both mates of a pair are on consecutive lines (as
   written by bam2ta), and are kept together so that each fragment ends up
   in one pseudoreplicate only
2. Randomly assign half of the pairs (rounded up) to pseudoreplicate 1 and
   the others to pseudoreplicate 2
3. Write each pair to its pseudoreplicate, keeping the original order

The random seed is the size of the uncompressed tagAlign file when
seed is 0, as in ENCODE, so that reruns give the same pseudoreplicates.
ENCODE shuffles with `shuf` and an openssl keystream as random source; this
script uses Python's random number generator instead, so the split differs
from that of the ENCODE pipeline but is equally random and reproducible.
Only the assignments (one byte per pair) are kept in memory, instead of the
whole file as `shuf` does.
"""

import gzip
import logging
import random

# Load Snakemake variables
tagalign = snakemake.input[0]
pr1 = snakemake.output["pr1"]
pr2 = snakemake.output["pr2"]
seed = int(snakemake.params["seed"])
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def open_tagalign(path):
    """Opens a gzipped or plain tagAlign file for reading in binary mode (as zcat -f)"""
    with open(path, "rb") as fh:
        is_gzipped = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if is_gzipped else open(path, "rb")


def gzip_writer(raw_fh):
    """
    Writes gzip to an open file without file name and time stamp in the
    header (as gzip -n), so that reruns give identical files
    """
    return gzip.GzipFile(
        fileobj=raw_fh, mode="wb", filename="", mtime=0, compresslevel=6
    )


try:
    # Count lines and uncompressed size (the default seed)
    # -----------------------------------------------------
    n_lines = 0
    n_bytes = 0
    with open_tagalign(tagalign) as fh:
        for line in fh:
            n_lines += 1
            n_bytes += len(line)
    if n_lines % 2 != 0:
        raise ValueError(
            f"{tagalign} has an odd number of lines ({n_lines}), "
            "so mates are not all paired"
        )
    n_pairs = n_lines // 2
    n_pr1 = (n_pairs + 1) // 2
    logging.info(f"{n_pairs} read pairs: {n_pr1} to pr1, {n_pairs - n_pr1} to pr2")

    if seed == 0:
        seed = n_bytes
    logging.info(f"Random seed for pseudoreplication: {seed}")

    # Randomly assign pairs to a pseudoreplicate
    # -----------------------------------------------------
    # 1: pr1, 0: pr2
    assignment = bytearray([1]) * n_pr1 + bytearray(n_pairs - n_pr1)
    random.Random(seed).shuffle(assignment)

    # Write the pairs to their pseudoreplicate
    # -----------------------------------------------------
    # GzipFile does not close the file it writes to, so the raw files are
    # opened (and closed) here
    with (
        open_tagalign(tagalign) as fh_in,
        open(pr1, "wb") as raw_pr1,
        open(pr2, "wb") as raw_pr2,
        gzip_writer(raw_pr1) as fh_pr1,
        gzip_writer(raw_pr2) as fh_pr2,
    ):
        for to_pr1 in assignment:
            mates = next(fh_in) + next(fh_in)
            (fh_pr1 if to_pr1 else fh_pr2).write(mates)
except Exception:
    logging.exception("Pseudoreplication failed")
    raise

logging.info("Done!")
