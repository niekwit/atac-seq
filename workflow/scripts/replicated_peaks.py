"""
Peaks of a condition that are reproducible in at least k of its n replicates
(not part of the ENCODE pipeline, whose optimal and conservative peak sets
come from one pair of replicates). k is peaks: min_replicate_fraction * n,
rounded up, and at least 2.

Peaks are peaks of the pooled sample (blacklist-filtered), identified by
chromosome, start, end and summit:
- overlap: the pooled peak overlaps a peak of at least k replicates by >= 50%
  (of either peak), the criterion of the naive overlap of ENCODE, counted
  with bedtools intersect -f 0.5 -F 0.5 -e. The pooled peak lines are written.
- IDR: the IDR peaks of the pairs of true replicates are subsets of the
  pooled peaks. A peak is kept if there are k replicates of which all pairs
  passed the IDR threshold for the peak. The lines of the first pair with
  the peak are written, with as score (column 5) the lowest score (scaled
  IDR) of these pairs.

With k = n (default), a peak is in the peak sets of all pairs of replicates;
with two replicates, this is the conservative peak set.
"""

import logging
import math
import re
import subprocess
from itertools import combinations

# Load Snakemake variables
method = snakemake.wildcards["method"]
replicate_peaks = list(snakemake.input["replicates"])
pair_peaks = list(snakemake.input["pairs"])
pooled = snakemake.input["pooled"]
pairs = snakemake.params["pairs"]
fraction = snakemake.params["min_fraction"]
output = snakemake.output[0]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

n = len(replicate_peaks)
k = min(n, max(2, math.ceil(fraction * n)))
logging.info(f"Peaks reproducible in at least {k} of {n} replicates ({method})")


def peak_key(fields):
    """Chromosome, start, end and summit of a narrowPeak line"""
    return (fields[0], fields[1], fields[2], fields[9])


def read_peaks(text):
    """Returns {key: narrowPeak fields}"""
    fields = [line.split("\t") for line in text.splitlines() if line.strip()]
    return {peak_key(p): p for p in fields}


with open(pooled) as fh:
    pooled_peaks = read_peaks(fh.read())
logging.info(f"{pooled}: {len(pooled_peaks)} peaks")

kept = []
if method == "overlap":
    # Number of replicates with an overlapping peak, for each pooled peak
    counts = dict.fromkeys(pooled_peaks, 0)
    for peaks in replicate_peaks:
        command = [
            "bedtools", "intersect", "-nonamecheck", "-u",
            "-f", "0.5", "-F", "0.5", "-e",
            "-a", pooled, "-b", peaks,
        ]  # fmt: skip
        overlapping = read_peaks(
            subprocess.run(command, check=True, capture_output=True, text=True).stdout
        )
        logging.info(f"{peaks}: overlaps {len(overlapping)} pooled peaks")
        for key in overlapping:
            counts[key] += 1
    kept = [fields for key, fields in pooled_peaks.items() if counts[key] >= k]
else:
    # IDR peaks of each pair of replicates (rep{i}_vs_rep{j})
    pair_reps = [tuple(map(int, re.findall(r"\d+", pair))) for pair in pairs]
    idr_peaks = {}
    for f, reps in zip(pair_peaks, pair_reps):
        with open(f) as fh:
            idr_peaks[reps] = read_peaks(fh.read())
        logging.info(f"{f}: {len(idr_peaks[reps])} peaks")

    # Keys in the order of the first pair, then of the other pairs
    keys = list(dict.fromkeys(key for peaks in idr_peaks.values() for key in peaks))
    for key in keys:
        passed = {reps for reps, peaks in idr_peaks.items() if key in peaks}
        for group in combinations(range(1, n + 1), k):
            group_pairs = list(combinations(group, 2))
            if all(pair in passed for pair in group_pairs):
                fields = list(next(idr_peaks[p][key] for p in pair_reps if p in passed))
                fields[4] = str(min(int(idr_peaks[p][key][4]) for p in group_pairs))
                kept.append(fields)
                break

with open(output, "w") as fh:
    fh.writelines("\t".join(fields) + "\n" for fields in kept)
logging.info(f"{len(kept)} peaks in at least {k} replicates written to {output}")
