"""
Annotates the TE consensus peaks with:

1. frac_unique: the fraction of the Tn5 insertions in the peak (all samples)
   that come from uniquely mapping read pairs. Peaks with a low fraction get
   their signal from multimappers, which sit at a random copy among their
   best alignments, so the copy (locus) of these peaks is uncertain, although
   the subfamily is usually right.
2. The RepeatMasker repeat with the largest overlap with the peak: subfamily,
   family, class and overlap (bp); empty for peaks outside repeats.

Writes a table with one row per peak (peak_id first) and the peaks as BED6
with frac_unique as score (0-1000).
"""

import logging
import os
import subprocess

import pandas as pd

# Load Snakemake variables
peaks_bed = snakemake.input["peaks"]
counts_multi = snakemake.input["multi"]
counts_unique = snakemake.input["unique"]
saf = snakemake.input["saf"]
subfamilies_file = snakemake.input["subfamilies"]
table_out = snakemake.output["table"]
bed_out = snakemake.output["bed"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def total_insertions(path):
    """Insertions per peak, summed over the samples (featureCounts table)"""
    counts = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    return counts.iloc[:, 5:].sum(axis=1)


peaks = pd.read_csv(
    peaks_bed,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "peak_id"],
    dtype={"chr": str},
)
multi = total_insertions(counts_multi)
unique = total_insertions(counts_unique)
peaks["insertions"] = peaks["peak_id"].map(multi).fillna(0).astype(int)
peaks["frac_unique"] = (
    (peaks["peak_id"].map(unique) / peaks["peak_id"].map(multi)).round(4)
).where(peaks["insertions"] > 0)

# Repeat with the largest overlap with each peak (SAF is 1-based)
repeats = pd.read_csv(saf, sep="\t", dtype={"Chr": str})
repeats_bed = f"{table_out}.repeats.tmp.bed"
pd.DataFrame(
    {
        "chr": repeats["Chr"],
        "start": repeats["Start"] - 1,
        "end": repeats["End"],
        "subfamily": repeats["GeneID"],
    }
).to_csv(repeats_bed, sep="\t", header=False, index=False)
command = ["bedtools", "intersect", "-wo", "-a", peaks_bed, "-b", repeats_bed]
overlaps = subprocess.run(command, check=True, capture_output=True, text=True).stdout
os.remove(repeats_bed)

columns = ["peak_id", "repeat_subfamily", "repeat_overlap_bp"]
rows = [line.split("\t") for line in overlaps.splitlines()]
overlaps = pd.DataFrame(
    [(r[3], r[7], int(r[8])) for r in rows], columns=columns
).sort_values(["peak_id", "repeat_overlap_bp"], ascending=[True, False])
best = overlaps.drop_duplicates("peak_id")

subfamilies = pd.read_csv(subfamilies_file, sep="\t")
best = best.merge(
    subfamilies[["subfamily", "family", "class"]].rename(
        columns={
            "subfamily": "repeat_subfamily",
            "family": "repeat_family",
            "class": "repeat_class",
        }
    ),
    on="repeat_subfamily",
    how="left",
)
peaks = peaks.merge(best, on="peak_id", how="left")
peaks["repeat_overlap_bp"] = peaks["repeat_overlap_bp"].fillna(0).astype(int)
logging.info(
    f"{len(peaks)} peaks, {(peaks['repeat_overlap_bp'] > 0).sum()} overlap a repeat; "
    f"median fraction of insertions from uniquely mapping reads: "
    f"{peaks['frac_unique'].median():.2f}"
)

columns = [
    "peak_id",
    "insertions",
    "frac_unique",
    "repeat_subfamily",
    "repeat_family",
    "repeat_class",
    "repeat_overlap_bp",
]
peaks[columns].to_csv(table_out, sep="\t", index=False)

# BED6 with frac_unique as score
bed = peaks[["chr", "start", "end", "peak_id"]].copy()
bed["score"] = (peaks["frac_unique"].fillna(0) * 1000).round().astype(int)
bed["strand"] = "."
bed.to_csv(bed_out, sep="\t", header=False, index=False)
logging.info("Done!")
