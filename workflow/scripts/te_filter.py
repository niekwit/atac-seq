"""
Filters the multimapper alignments (bowtie2 -k, read name order) of the TE
module:

1. Keep the best properly paired alignments of each read pair: in -k mode
   bowtie2 reports any valid alignment, not only the best ones, so only the
   alignments with the highest alignment score (AS of both mates) are kept
2. Remove PCR duplicates: read pairs whose set of alignments (chromosome,
   fragment start, fragment length and strand of each alignment) is identical
   to that of an earlier read pair. Using all alignments, instead of only the
   primary one (as Picard does), makes this work for multimappers, whose
   primary alignment bowtie2 picks at random.
3. Remove alignments on filter_chrs (e.g. MT)
4. Add the NH tag (number of best alignments of the read pair), which
   featureCounts uses to weight multimappers

Duplicates are found in a first pass over the BAM file, which stores a 64-bit
hash of each read pair's alignments (8 bytes per read pair), and removed in a
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
max_alignments = int(snakemake.params["max_alignments"])
threads = snakemake.threads
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def mate_pairs(records):
    """
    Returns the properly paired alignments of a read pair as (read 1, read 2)
    tuples, matching the mates on their positions
    """
    read2 = {}
    for r in records:
        if r.is_read2:
            key = (r.reference_id, r.next_reference_start, r.reference_start)
            read2.setdefault(key, []).append(r)
    pairs = []
    for r in records:
        if r.is_read1:
            key = (r.reference_id, r.reference_start, r.next_reference_start)
            if read2.get(key):
                pairs.append((r, read2[key].pop()))
    return pairs


def read_pairs(fh):
    """
    Yields, for each read pair with at least one properly paired alignment
    (in file order), its alignment pairs with the highest alignment score and
    the number of alignment pairs reported by bowtie2
    """
    for _, records in itertools.groupby(
        fh.fetch(until_eof=True), key=lambda r: r.query_name
    ):
        proper = [
            r
            for r in records
            if r.is_proper_pair and not r.is_unmapped and not r.mate_is_unmapped
        ]
        pairs = mate_pairs(proper)
        if pairs:
            scores = [r1.get_tag("AS") + r2.get_tag("AS") for r1, r2 in pairs]
            best = max(scores)
            yield [p for p, score in zip(pairs, scores) if score == best], len(pairs)


def fragments(pairs):
    """(chromosome, start, length, strand of read 1) of each alignment pair"""
    return tuple(
        sorted(
            (
                r1.reference_id,
                min(r1.reference_start, r1.next_reference_start),
                abs(r1.template_length),
                r1.is_reverse,
            )
            for r1, _ in pairs
        )
    )


# Pass 1: hash of the alignments of each read pair
# -----------------------------------------------------
logging.info(f"Pass 1: finding duplicates in {bam}")
hashes = array("Q")
mask = (1 << 64) - 1
with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=threads) as fh:
    for pairs, _ in read_pairs(fh):
        hashes.append(hash(fragments(pairs)) & mask)

hashes = np.frombuffer(hashes, dtype=np.uint64)
n_pairs = len(hashes)
# Keep the first read pair of each set of identical alignments
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
    "filter_chrs_only_pairs": 0,
    "kept_pairs": 0,
    "unique_pairs": 0,
    "multimapping_pairs": 0,
    "pairs_at_max_alignments": 0,
}
with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=threads) as fh_in:
    with pysam.AlignmentFile(bam_out, "wb", template=fh_in, threads=threads) as fh_out:
        for i, (pairs, n_reported) in enumerate(read_pairs(fh_in)):
            if not keep[i]:
                continue
            if n_reported >= max_alignments:
                stats["pairs_at_max_alignments"] += 1
            pairs = [p for p in pairs if p[0].reference_name not in filter_chrs]
            n = len(pairs)
            if n == 0:
                stats["filter_chrs_only_pairs"] += 1
                continue
            stats["kept_pairs"] += 1
            stats["unique_pairs" if n == 1 else "multimapping_pairs"] += 1
            for pair in pairs:
                for r in pair:
                    r.set_tag("NH", n, "i")
                    fh_out.write(r)

# Each read of a kept pair is one Tn5 insertion
stats["kept_reads"] = 2 * stats["kept_pairs"]
with open(stats_out, "w") as fh:
    fh.write("\t".join(stats) + "\n")
    fh.write("\t".join(str(v) for v in stats.values()) + "\n")
for key, value in stats.items():
    logging.info(f"{key}: {value}")
logging.info("Done!")
