"""
Filters the TE module alignment (bowtie2 default mode, one alignment per read
pair, read name order):

1. Keep properly paired alignments
2. Remove PCR duplicates: read pairs with the same chromosome, fragment start,
   fragment length and strand as an earlier read pair
3. Remove alignments on filter_chrs (e.g. MT)
4. Mark multimapping read pairs with NH:i:2, so that featureCounts without
   -M counts uniquely mapping read pairs only (other read pairs get NH:i:1).
   A read pair is multimapping when a mate has an equally good alternative
   alignment (XS >= AS) that its mate does not resolve (MAPQ < 2, which
   bowtie2 also gives to some uniquely aligned reads with low scores);
   bowtie2 then picked one of the best alignments at random.

Duplicates are found in a first pass over the BAM file, which stores a 64-bit
hash of each read pair's fragment (8 bytes per read pair), and removed in a
second pass.
"""

import itertools
import logging
from array import array

import numpy as np
import pysam

# Load Snakemake variables
bam = snakemake.input[0]
bam_out = snakemake.output["bam"]
stats_out = snakemake.output["stats"]
filter_chrs = set(snakemake.params["filter_chrs"])
threads = snakemake.threads
log = snakemake.log[0]

# bowtie2 gives MAPQ 0 or 1 to reads with an equally good alternative alignment
MIN_UNIQUE_MAPQ = 2

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def read_pairs(fh):
    """
    Yields the properly paired alignment (read 1, read 2) of each read pair
    that has one, in file order
    """
    for _, records in itertools.groupby(
        fh.fetch(until_eof=True), key=lambda r: r.query_name
    ):
        r1 = r2 = None
        for r in records:
            if r.is_secondary or r.is_supplementary or not r.is_proper_pair:
                continue
            if r.is_unmapped or r.mate_is_unmapped:
                continue
            if r.is_read1:
                r1 = r
            elif r.is_read2:
                r2 = r
        if r1 is not None and r2 is not None:
            yield r1, r2


def multimapping(r):
    """
    True if the mate has an equally good alternative alignment that the
    pairing does not resolve
    """
    return (
        r.mapping_quality < MIN_UNIQUE_MAPQ
        and r.has_tag("XS")
        and r.get_tag("XS") >= r.get_tag("AS")
    )


def fragment(r1):
    """Chromosome, start, length and strand of read 1 of an alignment pair"""
    return (
        r1.reference_id,
        min(r1.reference_start, r1.next_reference_start),
        abs(r1.template_length),
        r1.is_reverse,
    )


# Pass 1: hash of the fragment of each read pair
# -----------------------------------------------------
logging.info(f"Pass 1: finding duplicates in {bam}")
hashes = array("Q")
mask = (1 << 64) - 1
with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=threads) as fh:
    for r1, _ in read_pairs(fh):
        hashes.append(hash(fragment(r1)) & mask)

hashes = np.frombuffer(hashes, dtype=np.uint64)
n_pairs = len(hashes)
# Keep the first read pair of each fragment
_, first = np.unique(hashes, return_index=True)
keep = np.zeros(n_pairs, dtype=bool)
keep[first] = True
del hashes
logging.info(f"{n_pairs} properly paired read pairs, {n_pairs - len(first)} duplicates")

# Pass 2: write the filtered alignments with NH tags
# -----------------------------------------------------
logging.info(f"Pass 2: writing {bam_out}")
stats = {
    "proper_pairs": n_pairs,
    "duplicate_pairs": n_pairs - len(first),
    "filter_chrs_pairs": 0,
    "kept_pairs": 0,
    "unique_pairs": 0,
    "multimapping_pairs": 0,
}
with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=threads) as fh_in:
    with pysam.AlignmentFile(bam_out, "wb", template=fh_in, threads=threads) as fh_out:
        for i, (r1, r2) in enumerate(read_pairs(fh_in)):
            if not keep[i]:
                continue
            if r1.reference_name in filter_chrs:
                stats["filter_chrs_pairs"] += 1
                continue
            stats["kept_pairs"] += 1
            unique = not (multimapping(r1) or multimapping(r2))
            stats["unique_pairs" if unique else "multimapping_pairs"] += 1
            for r in (r1, r2):
                r.set_tag("NH", 1 if unique else 2, "i")
                fh_out.write(r)

# Each read of a kept pair is one Tn5 insertion
stats["kept_reads"] = 2 * stats["kept_pairs"]
with open(stats_out, "w") as fh:
    fh.write("\t".join(stats) + "\n")
    fh.write("\t".join(str(v) for v in stats.values()) + "\n")
for key, value in stats.items():
    logging.info(f"{key}: {value}")
logging.info("Done!")
