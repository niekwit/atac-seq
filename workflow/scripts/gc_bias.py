"""
GC bias of the filtered, deduplicated reads with Picard CollectGcBiasMetrics,
as in the ENCODE ATAC-seq pipeline (encode_task_gc_bias.py). The Picard
output is replotted as in ENCODE: normalised coverage, mean base quality and
the fraction of genomic windows at each GC percentage.
"""

import logging
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Load Snakemake variables
bam = snakemake.input["bam"]
fasta = snakemake.input["fasta"]
metrics = snakemake.output["metrics"]
summary = snakemake.output["summary"]
chart_pdf = snakemake.output["chart_pdf"]
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

# GC bias metrics
# -----------------------------------------------------
command = [
    "picard", f"-Xmx{java_heap}", "CollectGcBiasMetrics",
    "-R", fasta,
    "-I", bam,
    "-O", metrics,
    "-CHART", chart_pdf,
    "-S", summary,
    "--ASSUME_SORTED", "false",
    "--VERBOSITY", "ERROR",
    "--QUIET", "true",
]  # fmt: skip
logging.info(f"Running: {' '.join(command)}")
with open(log, "a") as log_fh:
    subprocess.run(command, check=True, stderr=log_fh)

# Plot
# -----------------------------------------------------
data = pd.read_table(metrics, comment="#")

fig, ax = plt.subplots()
ax.set_xlim(0, 100)
lines = ax.plot(
    data["GC"], data["NORMALIZED_COVERAGE"], label="Normalized coverage", color="r"
)
ax.set_xlabel("GC%")
ax.set_ylabel("Normalized coverage")

ax2 = ax.twinx()
lines += ax2.plot(
    data["GC"],
    data["MEAN_BASE_QUALITY"],
    label="Mean base quality at GC%",
    color="b",
)
ax2.set_ylabel("Mean base quality at GC%")

ax3 = ax.twinx()
lines += ax3.plot(
    data["GC"],
    data["WINDOWS"] / data["WINDOWS"].sum(),
    label="Windows at GC%",
    color="g",
)
ax3.get_yaxis().set_visible(False)

ax.legend(lines, [line.get_label() for line in lines], loc="best")
fig.tight_layout()
fig.savefig(plot_png)
logging.info("Done!")
