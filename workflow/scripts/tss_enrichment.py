"""
Calculates the TSS enrichment score as in the ENCODE ATAC-seq pipeline
(encode_task_tss_enrich.py, which uses metaseq), so that it can be graded
against the ENCODE standards:

1. Each read of the filtered, deduplicated BAM file is shifted upstream by
   half the read length, so that it is centred on its 5' end (the Tn5
   cutting site)
2. Read coverage is sampled at 400 evenly spaced points in a window of
   +/- 2 kb around each TSS (minus strand TSSs are flipped)
3. The coverage is averaged over all TSSs and divided by the background:
   the mean of the outer 100 bp at both ends of the window
4. The TSS enrichment score is the maximum of this normalised profile

metaseq selects the reads of a window before shifting them, from the window
"padded" by the shift. As the shift is negative, this padded window is
smaller than the window itself, so reads that are shifted into the outer
half read length of the window are not counted. This lowers the background
and raises the score (by about a third for 100 bp reads). It is reproduced
here, as the ENCODE thresholds are based on scores calculated this way.
"""

import gzip
import logging

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam

# Load Snakemake variables
bam_file = snakemake.input["bam"]
tss_file = snakemake.input["tss"]
fastq = snakemake.input["fastq"]
score_out = snakemake.output["score"]
profile_out = snakemake.output["profile"]
plot_out = snakemake.output["plot"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

# ENCODE settings
BP_EDGE = 2000  # window of +/- 2 kb around the TSS
BINS = 400  # coverage sampled at 400 points (every ~10 bp)
EDGE_BINS = int(100 / (2 * BP_EDGE / BINS))  # 10 bins = outer 100 bp


def read_length(fastq, n_reads=1_000_000):
    """Longest read in the first million reads, as ENCODE's get_read_length"""
    max_length = 0
    with gzip.open(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                max_length = max(max_length, len(line.strip()))
            if i >= 4 * n_reads:
                break
    return max_length


def load_tss(tss_file):
    """Returns {chrom: [(position, strand), ...]} of 0-based TSS positions"""
    tss = {}
    with open(tss_file) as fh:
        for line in fh:
            chrom, start, _, _, _, strand = line.rstrip("\n").split("\t")[:6]
            tss.setdefault(chrom, []).append((int(start), strand))
    return tss


def load_reads(bam, chrom):
    """Start, end and strand of the reads of a chromosome, sorted by start"""
    starts, ends, reverse = [], [], []
    for read in bam.fetch(chrom):
        if read.is_unmapped:
            continue
        starts.append(read.reference_start)
        ends.append(read.reference_end)
        reverse.append(read.is_reverse)
    return (
        np.array(starts, dtype=np.int64),
        np.array(ends, dtype=np.int64),
        np.array(reverse, dtype=bool),
    )


def window_coverage(starts, ends, reverse, max_len, start, stop, shift):
    """
    Read coverage of the window [start, stop), as metaseq _local_coverage
    with shift_width=shift
    """
    size = stop - start

    # Reads overlapping the padded window (BAM is sorted by start, so reads
    # overlapping it start between padded_start - max_len and padded_stop)
    padded_start = max(start - shift, 0)
    padded_stop = stop + shift
    i0 = np.searchsorted(starts, padded_start - max_len, side="left")
    i1 = np.searchsorted(starts, padded_stop, side="left")
    s, e, r = starts[i0:i1], ends[i0:i1], reverse[i0:i1]
    keep = e > padded_start
    s, e, r = s[keep], e[keep], r[keep]

    # Shift reads towards their 5' end (upstream for a negative shift)
    offset = np.where(r, -shift, shift)
    s, e = s + offset - start, e + offset - start

    # Clip to the window and skip reads shifted outside it
    s, e = np.maximum(s, 0), np.minimum(e, size)
    inside = (s < size) & (e > s)
    s, e = s[inside], e[inside]

    # Coverage from a difference array: +1 at each read start, -1 at each end
    diff = np.bincount(s, minlength=size + 1) - np.bincount(e, minlength=size + 1)
    return np.cumsum(diff[:size])


# Read length and shift
# -----------------------------------------------------
read_len = read_length(fastq)
shift = int(-read_len / 2)
logging.info(f"Read length: {read_len}; reads shifted by {shift} bp")

# Coverage profile around the TSSs
# -----------------------------------------------------
tss = load_tss(tss_file)

# Window of 2 * BP_EDGE + 1 bp (bedtools slop -b BP_EDGE of a 1 bp TSS),
# sampled at BINS evenly spaced points by linear interpolation (metaseq rebin)
window = 2 * BP_EDGE + 1
sample_at = np.linspace(0, window - 1, BINS)
lower = np.floor(sample_at).astype(int)
upper = np.minimum(lower + 1, window - 1)
weight = sample_at - lower

profile_sum = np.zeros(BINS)
n_tss = n_skipped = n_reads = 0
with pysam.AlignmentFile(bam_file) as bam:
    chrom_lengths = dict(zip(bam.references, bam.lengths))
    for chrom, sites in tss.items():
        if chrom not in chrom_lengths:
            n_skipped += len(sites)
            continue
        length = chrom_lengths[chrom]
        starts, ends, reverse = load_reads(bam, chrom)
        max_len = int((ends - starts).max()) if len(starts) else 0
        n_reads += len(starts)

        for position, strand in sites:
            start = position - BP_EDGE
            if start < 0 or start + window > length:
                # Window runs past the chromosome end
                n_skipped += 1
                continue
            values = window_coverage(
                starts, ends, reverse, max_len, start, start + window, shift
            )
            row = values[lower] * (1 - weight) + values[upper] * weight
            if strand == "-":
                row = row[::-1]
            profile_sum += row
            n_tss += 1
        logging.info(f"{chrom}: {len(starts)} reads, {len(sites)} TSSs")

logging.info(f"Used {n_tss} TSSs ({n_skipped} skipped) and {n_reads} reads")

# Normalise to the background at the window edges
# -----------------------------------------------------
profile = profile_sum / n_tss
background = (profile[:EDGE_BINS].sum() + profile[-EDGE_BINS:].sum()) / (2 * EDGE_BINS)
profile = profile / background
tss_enrichment = profile.max()
logging.info(f"TSS enrichment: {tss_enrichment}")

with open(score_out, "w") as fh:
    fh.write(f"{tss_enrichment}\n")

positions = np.linspace(-BP_EDGE, BP_EDGE, BINS)
with open(profile_out, "w") as fh:
    fh.write("distance_to_tss\tenrichment\n")
    for x, y in zip(positions, profile):
        fh.write(f"{x:.1f}\t{y:.4f}\n")

# Plot
# -----------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(positions, profile, color="firebrick")
ax.axvline(0, linestyle=":", color="black")
ax.set_xlabel("Distance from TSS (bp)")
ax.set_ylabel("TSS enrichment")
ax.set_title(f"TSS enrichment: {tss_enrichment:.2f}")
fig.tight_layout()
fig.savefig(plot_out)
logging.info("Done!")
