"""
Calls peaks on a Tn5 shifted tagAlign file with MACS2, as in the ENCODE
ATAC-seq pipeline (encode_task_macs2_atac.py):

1. MACS2 callpeak on the reads, extended around the Tn5 cutting sites
2. Peaks sorted by p-value (strongest first) and renamed Peak_1, Peak_2, ...
3. Only the top cap_num_peak peaks are kept (ENCODE: 300,000), so that
   IDR has a large, lenient peak set to work with
4. Peaks clipped to the chromosome sizes

The pileup and lambda bedGraphs written by MACS2 are kept for the
signal tracks (macs2_signal_track rule).
"""

import logging
import os
import subprocess

# Load Snakemake variables
tagalign = snakemake.input["ta"]
chrom_sizes_file = snakemake.input["chrom_sizes"]
narrowpeak = snakemake.output["narrowpeak"]
treat_bdg = snakemake.output["treat"]
control_bdg = snakemake.output["control"]
gsize = snakemake.params["gsize"]
pval_thresh = snakemake.params["pval_thresh"]
smooth_win = int(snakemake.params["smooth_win"])
cap_num_peak = int(snakemake.params["cap_num_peak"])
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# MACS2 writes {outdir}/{name}_peaks.narrowPeak etc.
outdir = os.path.dirname(narrowpeak)
name = os.path.basename(narrowpeak).removesuffix(".narrowPeak")
prefix = os.path.join(outdir, name)

# Call peaks
# -----------------------------------------------------
# Reads are shifted by -smooth_win/2 and extended to smooth_win bp, so that
# each read becomes a smooth_win bp fragment centred on the Tn5 cutting site.
# -B --SPMR: also write pileup bedGraphs, normalised per million reads
# --keep-dup all: duplicates were already removed
shift = -int(round(smooth_win / 2))
command = [
    "macs2", "callpeak",
    "-t", tagalign,
    "-f", "BED",
    "-n", name,
    "--outdir", outdir,
    "-g", str(gsize),
    "-p", str(pval_thresh),
    "--shift", str(shift),
    "--extsize", str(smooth_win),
    "--nomodel",
    "-B", "--SPMR",
    "--keep-dup", "all",
    "--call-summits",
]  # fmt: skip
logging.info(f"Running: {' '.join(command)}")
with open(log, "a") as log_fh:
    subprocess.run(command, check=True, stderr=log_fh)

# Sort, rename and cap peaks
# -----------------------------------------------------
# narrowPeak columns (0-based): 0 chrom, 1 start, 2 end, 3 name, 4 score,
# 5 strand, 6 signal value, 7 -log10(p-value), 8 -log10(q-value), 9 summit
with open(f"{prefix}_peaks.narrowPeak") as fh:
    peaks = [line.rstrip("\n").split("\t") for line in fh]
logging.info(f"MACS2 called {len(peaks)} peaks")

# Sort by -log10(p-value), highest first. Ties are ordered by the whole line,
# as ENCODE's `LC_COLLATE=C sort -k 8gr,8gr`: sort by line first, then do a
# stable sort on the p-value.
peaks.sort(key=lambda p: "\t".join(p))
peaks.sort(key=lambda p: float(p[7]), reverse=True)

peaks = peaks[:cap_num_peak]
logging.info(f"Kept the top {len(peaks)} peaks (cap: {cap_num_peak})")

for i, peak in enumerate(peaks, start=1):
    peak[3] = f"Peak_{i}"
    start = max(int(peak[1]), 0)
    end = max(int(peak[2]), 0)
    peak[1], peak[2] = str(start), str(end)
    # Summit offset -1 means no summit: use the peak centre
    if peak[9] == "-1":
        peak[9] = str(start + int((end - start + 1) / 2.0))

# Clip peaks to chromosome sizes
# -----------------------------------------------------
# As `bedClip -truncate`: peak ends past the chromosome end are truncated;
# peaks on unknown chromosomes or starting past the chromosome end are dropped.
with open(chrom_sizes_file) as fh:
    chrom_sizes = {
        chrom: int(size) for chrom, size in (line.split()[:2] for line in fh)
    }

clipped = []
for peak in peaks:
    size = chrom_sizes.get(peak[0])
    if size is None or int(peak[1]) >= size:
        logging.warning(f"Dropped peak outside chromosome: {' '.join(peak[:3])}")
        continue
    peak[2] = str(min(int(peak[2]), size))
    clipped.append(peak)

with open(narrowpeak, "w") as fh:
    fh.writelines("\t".join(peak) + "\n" for peak in clipped)
logging.info(f"Wrote {len(clipped)} peaks to {narrowpeak}")

# Clean up
# -----------------------------------------------------
# Keep the bedGraphs for the signal tracks, remove the other MACS2 output
if f"{prefix}_treat_pileup.bdg" != treat_bdg:
    os.replace(f"{prefix}_treat_pileup.bdg", treat_bdg)
    os.replace(f"{prefix}_control_lambda.bdg", control_bdg)
for suffix in ["_peaks.xls", "_summits.bed", "_peaks.narrowPeak"]:
    os.remove(f"{prefix}{suffix}")
logging.info("Done!")
