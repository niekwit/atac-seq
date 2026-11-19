"""
Signal of H3K9me3 and DNA methylation (BigWig files of the control and the
treatment of a contrast) over all copies of each repeat subfamily
(resources/te/repeats.saf), for the integration of all data sets of the
contrast (TE module).

Per subfamily and BigWig file, the signal is pooled over the bases of all
copies:

- coverage signal (missing_as_zero, e.g. H3K9me3): the summed signal of all
  bases of the copies divided by their total length (bases without signal
  are 0)
- methylation (not missing_as_zero, per CpG): the summed methylation of the
  bases with a value divided by the number of those bases (the mean over the
  covered CpGs); NA without a covered base

The sums are exact: per chromosome, the BigWig intervals are read in windows
and integrated over each copy with cumulative sums. Intervals with a NaN
value count as missing. If the chromosomes of a BigWig file are UCSC names
(chr1), the Ensembl names of the copies (1, MT) are converted.

Output columns: subfamily, n_copies, bp, and per BigWig file (names from
params: names) {name} (the mean signal) and {name}_covered_bp.
"""

import logging

import numpy as np
import pandas as pd
import pyBigWig

saf_file = snakemake.input["saf"]
bigwigs = list(snakemake.input["bigwigs"])
names = list(snakemake.params["names"])
missing_as_zero = list(snakemake.params["missing_as_zero"])
out_file = snakemake.output[0]
log = snakemake.log[0]

# Copies whose starts fall in a window are integrated together (bp)
WINDOW = 5_000_000

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# SAF: 1-based starts
saf = pd.read_csv(saf_file, sep="\t", dtype={"Chr": str, "GeneID": str})
copies = pd.DataFrame(
    {
        "chr": saf["Chr"],
        "start": saf["Start"].astype(np.int64) - 1,
        "end": saf["End"].astype(np.int64),
        "subfamily": saf["GeneID"],
    }
).sort_values(["chr", "start"], ignore_index=True)
logging.info(f"{len(copies)} copies of {copies['subfamily'].nunique()} subfamilies")


def integrate(bw, chrom, starts, ends):
    """
    Summed signal and number of bases with a value of each region [start,
    end) of a chromosome of a BigWig file (starts sorted)
    """
    sums = np.zeros(len(starts))
    covered = np.zeros(len(starts))
    first = 0
    while first < len(starts):
        last = np.searchsorted(starts, starts[first] + WINDOW, side="left")
        last = max(last, first + 1)
        lo, hi = int(starts[first]), int(ends[first:last].max())
        intervals = bw.intervals(chrom, lo, hi)
        if intervals:
            iv = np.array(intervals, dtype=float)
            iv = iv[np.isfinite(iv[:, 2])]
        if intervals and len(iv):
            s, e, v = iv[:, 0], iv[:, 1], iv[:, 2]
            # Cumulative signal and bases up to the end of each interval
            cum_signal = np.concatenate([[0.0], np.cumsum(v * (e - s))])
            cum_bases = np.concatenate([[0.0], np.cumsum(e - s)])

            def up_to(x, weights, cum):
                """Integral from the chromosome start to x"""
                # Intervals that end at or before x are complete
                i = np.searchsorted(e, x, side="right")
                total = cum[i]
                # The next interval is partly before x if it starts before x
                partial = i < len(s)
                j = np.where(partial, i, 0)
                extra = np.where(partial & (s[j] < x), weights[j] * (x - s[j]), 0.0)
                return total + extra

            r_start, r_end = starts[first:last], ends[first:last]
            sums[first:last] = up_to(r_end, v, cum_signal) - up_to(
                r_start, v, cum_signal
            )
            ones = np.ones_like(v)
            covered[first:last] = up_to(r_end, ones, cum_bases) - up_to(
                r_start, ones, cum_bases
            )
        first = last
    return sums, covered


result = (
    copies.assign(bp=copies["end"] - copies["start"])
    .groupby("subfamily")
    .agg(n_copies=("bp", "size"), bp=("bp", "sum"))
)
for path, name, as_zero in zip(bigwigs, names, missing_as_zero):
    bw = pyBigWig.open(path)
    bw_chroms = bw.chroms()
    # Ensembl names of the copies to the names in the BigWig file
    if set(copies["chr"]) & set(bw_chroms):
        rename = {c: c for c in copies["chr"].unique()}
    else:
        rename = {c: "chrM" if c == "MT" else f"chr{c}" for c in copies["chr"].unique()}
        logging.info(f"{name}: UCSC chromosome names")
    sums = np.zeros(len(copies))
    covered = np.zeros(len(copies))
    for chrom, rows in copies.groupby("chr", sort=False).indices.items():
        bw_chrom = rename[chrom]
        if bw_chrom not in bw_chroms:
            continue
        # Copies sorted by start within the chromosome; ends within its length
        length = bw_chroms[bw_chrom]
        starts = copies["start"].to_numpy()[rows].clip(0, length)
        ends = copies["end"].to_numpy()[rows].clip(0, length)
        sums[rows], covered[rows] = integrate(bw, bw_chrom, starts, ends)
    bw.close()
    by_subfamily = (
        pd.DataFrame(
            {"sum": sums, "covered": covered, "subfamily": copies["subfamily"]}
        )
        .groupby("subfamily")[["sum", "covered"]]
        .sum()
    )
    if as_zero:
        result[name] = by_subfamily["sum"] / result["bp"]
    else:
        result[name] = by_subfamily["sum"] / by_subfamily["covered"].replace(0, np.nan)
    result[f"{name}_covered_bp"] = by_subfamily["covered"].astype(np.int64)
    logging.info(
        f"{name}: {(by_subfamily['covered'] > 0).sum()} subfamilies with "
        f"signal; {covered.sum() / result['bp'].sum():.2%} of the copy bases "
        "have a value"
    )
result.reset_index().to_csv(out_file, sep="\t", index=False)
logging.info("Done!")
