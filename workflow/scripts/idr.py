"""
Runs IDR (Irreproducible Discovery Rate) on two peak sets and keeps the
reproducible peaks, as in the ENCODE ATAC-seq pipeline (encode_task_idr.py).

Why
---
Peak callers report many peaks at a lenient threshold (here MACS2 with
p < 0.01, top 300,000 peaks), most of which are noise. Rather than picking
an arbitrary significance cutoff for each sample, IDR (Li et al., Ann. Appl.
Stat. 2011) uses the agreement between two replicates: real peaks are found
in both and have a similar rank (by p-value) in each, whereas noise peaks
are ranked inconsistently. IDR fits a copula mixture model with a
reproducible and an irreproducible component to the ranks of the matched
peaks, and gives each peak the probability that it belongs to the
irreproducible component. The global IDR of a peak is the expected fraction
of irreproducible peaks among all peaks at least as reproducible, so a
threshold of 0.05 is like a 5% false discovery rate for reproducibility.

The workflow runs this script on three kinds of pairs (see the idr rule):
- true replicates (repX_vs_repY)
- the two pseudoreplicates of each replicate (repX-pr1_vs_repX-pr2):
  self-consistency
- the pooled pseudoreplicates of a condition (pooled-pr1_vs_pooled-pr2)
The number of IDR peaks of these pairs gives the optimal and conservative
peak sets, and the rescue and self-consistency ratios
(reproducibility rule).

What
----
1. IDR on the two peak sets, with the peaks of the pooled sample (pooled
   replicates or pooled pseudoreplicates) as the list of peaks to test:
   - --input-file-type narrowPeak, --rank p.value: peaks are ranked by
     their MACS2 -log10(p-value) (column 8)
   - --peak-list: the pooled peaks are the "oracle" peaks; each is matched
     with the overlapping peak of each of the two samples, and peaks of the
     pooled sample without a match in both samples are left out
   - --soft-idr-threshold: peaks with an IDR below this value are
     summarised in the IDR log (all peaks are written to the output)
   - --use-best-multisummit-IDR: MACS2 --call-summits reports a peak with
     several summits as several lines with the same coordinates; these get
     the best IDR of the group
   - --plot: diagnostic plots of the ranks and IDR values (png)
2. Peaks are clipped to the chromosome sizes, as `bedClip -truncate`:
   negative starts become 0, ends past the chromosome end are truncated,
   and empty peaks (start >= end) are removed. Unlike bedClip, peaks that
   start past the chromosome end are removed instead of written with
   start > end. A chromosome that is not in the chromosome sizes file is an
   error, as in bedClip.
3. All clipped peaks are written to the gzipped unthresholded file
   (all IDR output columns, see below).
4. Peaks with a global IDR <= threshold are kept: column 12 (-log10 global
   IDR) >= -log10(threshold). Duplicate lines are removed, peaks are sorted
   by -log10(p-value) (column 8), highest first, and the first 10 columns
   (narrowPeak) are written.

IDR output columns (narrowPeak input):
 1-10  merged peak in narrowPeak format: chrom, start, end, name, score
       (scaled IDR: min(int(-125 * log2(IDR)), 1000)), strand,
       signalValue, -log10(p-value), -log10(q-value), summit offset
 11    -log10(local IDR)
 12    -log10(global IDR)
 13-16 sample 1 start, end, signalValue, summit offset
 17-20 sample 2 start, end, signalValue, summit offset

Differences from ENCODE:
- ENCODE sorts with `sort | uniq | sort -grk8,8` in the system locale; peaks
  with the same p-value are here ordered by the whole line in reverse byte
  order (as `LC_ALL=C sort -grk8,8`), so that the order does not depend on
  the locale.
- Peaks starting past the chromosome end are removed (see 2.); they do not
  occur in practice, as the MACS2 peaks are already clipped.
"""

import gzip
import logging
import math
import os
import subprocess

# Load Snakemake variables
peak1 = snakemake.input["peak1"]
peak2 = snakemake.input["peak2"]
pooled = snakemake.input["pooled"]
chrom_sizes_file = snakemake.input["chrom_sizes"]
narrowpeak = snakemake.output["narrowpeak"]
unthresholded = snakemake.output["unthresholded"]
plot = snakemake.output["plot"]
idr_log = snakemake.output["idr_log"]
threshold = float(snakemake.params["threshold"])
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

idr_out = f"{narrowpeak}.unthresholded-peaks.txt"


def clip_peak(peak, chrom_sizes):
    """
    Clips a peak (list of columns) to its chromosome, as bedClip -truncate.
    Returns None for peaks that are removed.
    """
    chrom = peak[0]
    if chrom not in chrom_sizes:
        raise ValueError(f"Chromosome {chrom} is not in {chrom_sizes_file}")
    size = chrom_sizes[chrom]
    start, end = int(peak[1]), int(peak[2])
    if end < 0:
        raise ValueError(f"Negative peak end: {' '.join(peak[:3])}")
    start = max(start, 0)
    if start >= end or start >= size:
        logging.warning(f"Removed peak outside chromosome: {' '.join(peak[:3])}")
        return None
    peak[1], peak[2] = str(start), str(min(end, size))
    return peak


try:
    # Run IDR
    # -----------------------------------------------------
    command = [
        "idr",
        "--samples", peak1, peak2,
        "--peak-list", pooled,
        "--input-file-type", "narrowPeak",
        "--output-file", idr_out,
        "--rank", "p.value",
        "--soft-idr-threshold", str(threshold),
        "--plot",
        "--use-best-multisummit-IDR",
        "--log-output-file", idr_log,
    ]  # fmt: skip
    logging.info(f"Running: {' '.join(command)}")
    with open(log, "a") as log_fh:
        subprocess.run(command, check=True, stderr=log_fh)
    os.replace(f"{idr_out}.png", plot)

    # Clip peaks to the chromosome sizes
    # -----------------------------------------------------
    with open(chrom_sizes_file) as fh:
        chrom_sizes = {
            chrom: int(size) for chrom, size in (line.split()[:2] for line in fh)
        }
    with open(idr_out) as fh:
        peaks = [line.rstrip("\n").split("\t") for line in fh]
    logging.info(f"IDR reported {len(peaks)} peaks")
    peaks = [p for p in (clip_peak(peak, chrom_sizes) for peak in peaks) if p]

    # Write all peaks (gzip -n: no file name and time stamp in the header)
    with open(unthresholded, "wb") as raw_fh, gzip.GzipFile(
        fileobj=raw_fh, mode="wb", filename="", mtime=0, compresslevel=6
    ) as fh:
        fh.writelines(("\t".join(peak) + "\n").encode() for peak in peaks)

    # Keep peaks with global IDR <= threshold
    # -----------------------------------------------------
    neg_log10_threshold = -math.log10(threshold)
    passed = {
        "\t".join(peak[:12]) for peak in peaks if float(peak[11]) >= neg_log10_threshold
    }
    # Sort by -log10(p-value), highest first; ties by the whole line (reverse)
    passed = sorted(
        passed, key=lambda line: (float(line.split("\t")[7]), line), reverse=True
    )
    with open(narrowpeak, "w") as fh:
        fh.writelines("\t".join(line.split("\t")[:10]) + "\n" for line in passed)
    logging.info(f"{len(passed)} peaks with IDR <= {threshold} written to {narrowpeak}")
    if not passed:
        logging.warning(
            "No IDR peaks found. The IDR threshold might be too stringent "
            "or replicates have very poor concordance."
        )
except Exception:
    logging.exception("IDR failed")
    raise
finally:
    if os.path.exists(idr_out):
        os.remove(idr_out)

logging.info("Done!")
