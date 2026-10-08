"""
Fragment length distribution and nucleosomal QC, as in the ENCODE ATAC-seq
pipeline (encode_task_fraglen_stat_pe.py):

1. Picard CollectInsertSizeMetrics on the first 5 million reads of the
   filtered, deduplicated BAM file (fragments up to 1 kb)
2. QC checks on the fragment length histogram (ENCODE standards: a
   nucleosome free region (NFR) and a mononucleosome peak should be present):
   - Fraction of fragments in the NFR (< 150 bp) >= 0.4
   - NFR / mononucleosome (150-300 bp) fragments >= 2.5
   - Presence of a NFR peak (20-90 bp), a mononucleosome peak (120-250 bp)
     and a dinucleosome peak (300-500 bp)
"""

import logging
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks_cwt

# Load Snakemake variables
bam = snakemake.input["bam"]
metrics = snakemake.output["metrics"]
histogram_pdf = snakemake.output["histogram_pdf"]
qc_out = snakemake.output["qc"]
plot_png = snakemake.output["plot"]
java_heap = snakemake.params["java_heap"]
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

NFR_UPPER_LIMIT = 150
MONO_NUC_LOWER_LIMIT = 150
MONO_NUC_UPPER_LIMIT = 300


def read_picard_histogram(metrics_file):
    """Returns the insert size histogram (size, count) of Picard's metrics"""
    with open(metrics_file) as fh:
        for line in fh:
            if line.startswith("## HISTOGRAM"):
                break
        return np.loadtxt(fh, skiprows=1, ndmin=2)


# Insert size histogram
# -----------------------------------------------------
command = [
    "picard", f"-Xmx{java_heap}", "CollectInsertSizeMetrics",
    "-I", bam,
    "-O", metrics,
    "-H", histogram_pdf,
    "-W", "1000",
    "--STOP_AFTER", "5000000",
    "--VERBOSITY", "ERROR",
    "--QUIET", "true",
]  # fmt: skip
logging.info(f"Running: {' '.join(command)}")
with open(log, "a") as log_fh:
    subprocess.run(command, check=True, stderr=log_fh)

data = read_picard_histogram(metrics)
sizes, counts = data[:, 0], data[:, 1]

# QC checks
# -----------------------------------------------------
results = []  # (metric, pass, value)

nfr = counts[sizes < NFR_UPPER_LIMIT].sum()
frac_nfr = nfr / counts.sum()
results.append(("Fraction of reads in NFR", frac_nfr >= 0.4, frac_nfr))

mono_nuc = counts[
    (sizes > MONO_NUC_LOWER_LIMIT) & (sizes <= MONO_NUC_UPPER_LIMIT)
].sum()
nfr_vs_mono_nuc = nfr / mono_nuc
results.append(("NFR / mono-nuc reads", nfr_vs_mono_nuc >= 2.5, nfr_vs_mono_nuc))

# Peaks in the histogram. find_peaks_cwt returns indices into the histogram,
# so (as in ENCODE) the ranges are shifted by the smallest fragment size
peaks = find_peaks_cwt(counts, np.array([25]))
start = sizes[0]
for metric, lower, upper in [
    ("Presence of NFR peak", 20, 90),
    ("Presence of Mono-Nuc peak", 120, 250),
    ("Presence of Di-Nuc peak", 300, 500),
]:
    found = any(lower - start <= peak <= upper - start for peak in peaks)
    results.append((metric, found, "OK" if found else f"none in [{lower}, {upper}]"))

with open(qc_out, "w") as fh:
    fh.write("metric\tpass\tvalue\n")
    for metric, qc_pass, value in results:
        fh.write(f"{metric}\t{qc_pass}\t{value}\n")
for result in results:
    logging.info("\t".join(str(x) for x in result))

# Plot
# -----------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(sizes, counts, width=1, color="steelblue")
ax.set_xlim(0, 1000)
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Fragments")
fig.tight_layout()
fig.savefig(plot_png)
logging.info("Done!")
