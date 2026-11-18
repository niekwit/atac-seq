"""
Profile and heatmap of the signal of the control and the treatment of a
contrast over groups of regions (deepTools computeMatrix, plotProfile and
plotHeatmap; one panel per group, one line per condition):

- ATAC-seq Tn5 insertions (CPM) over the H3K9me3 regions that are lost,
  gained or unchanged (scale-regions: regions scaled to body_length bp)
- H3K9me3 signal around the centre of the more, less and unchanged
  accessible peaks (reference-point)

Regions come from a table with chr, start and end columns, grouped by
group_column. If the chromosomes of the BigWig files are UCSC names (chr1),
the Ensembl names of the regions (1, MT) are converted. Writes a PDF with a
message if no group has regions.
"""

import gzip
import logging
import os
import subprocess
import tempfile

import pandas as pd
import pyBigWig

table_file = snakemake.input["table"]
bigwigs = list(snakemake.input["bigwigs"])
group_column = snakemake.params["group_column"]
groups = dict(snakemake.params["groups"])
mode = snakemake.params["mode"]
flank = int(snakemake.params["flank"])
body_length = int(snakemake.params.get("body_length", 0) or 0)
labels = list(snakemake.params["labels"])
colours = list(snakemake.params["colours"])
y_label = snakemake.params["y_label"]
region_name = snakemake.params["region_name"]
extra_matrix = snakemake.params["extra_computematrix"]
extra_plot = snakemake.params["extra_plotprofile"]
out_matrix = snakemake.output["matrix"]
out_profile = snakemake.output["profile"]
out_heatmap = snakemake.output["heatmap"]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

table = pd.read_csv(table_file, sep="\t", dtype={"chr": str})

# Chromosome names of the BigWig files
bigwig_chroms = set()
for path in bigwigs:
    bw = pyBigWig.open(path)
    bigwig_chroms |= set(bw.chroms())
    bw.close()
if not set(table["chr"]) & bigwig_chroms:
    ucsc = ("chr" + table["chr"]).replace({"chrMT": "chrM"})
    if set(ucsc) & bigwig_chroms:
        logging.info("BigWig files have UCSC chromosome names: regions converted")
        table["chr"] = ucsc
    else:
        logging.warning("No chromosome of the regions is in the BigWig files")

tmp = tempfile.mkdtemp(dir=os.path.dirname(out_matrix))
# Heatmap labels without the numbers, which are in the profile
beds, region_labels, heatmap_labels = [], [], []
for value, label in groups.items():
    selected = table[table[group_column] == value]
    logging.info(f"{label}: {len(selected)} regions")
    if selected.empty:
        continue
    path = os.path.join(tmp, f"{value.replace(' ', '_')}.bed")
    selected[["chr", "start", "end"]].sort_values(["chr", "start"]).to_csv(
        path, sep="\t", header=False, index=False
    )
    beds.append(path)
    region_labels.append(f"{label} (n = {len(selected)})")
    heatmap_labels.append(label)


def run(command):
    logging.info(" ".join(command))
    with open(log, "a") as fh:
        subprocess.run(command, check=True, stdout=fh, stderr=fh)


if not beds:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logging.warning(f"No {region_name}: no profile")
    for path in (out_profile, out_heatmap):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, f"No {region_name}", ha="center")
        fig.savefig(path)
        plt.close(fig)
    with gzip.open(out_matrix, "wt"):
        pass
else:
    if mode == "scale-regions":
        layout = [
            "scale-regions",
            "--regionBodyLength",
            str(body_length),
            "--unscaled5prime",
            "0",
            "--unscaled3prime",
            "0",
        ]
        point_labels = ["--startLabel", "start", "--endLabel", "end"]
    else:
        layout = ["reference-point", "--referencePoint", "center"]
        point_labels = ["--refPointLabel", "centre"]
    run(
        [
            "computeMatrix",
            *layout,
            "-R",
            *beds,
            "-S",
            *bigwigs,
            "--samplesLabel",
            *labels,
            "-b",
            str(flank),
            "-a",
            str(flank),
            "--binSize",
            "50",
            "--missingDataAsZero",
            "-p",
            str(snakemake.threads),
            "-o",
            out_matrix,
            *extra_matrix.split(),
        ]
    )
    run(
        [
            "plotProfile",
            "-m",
            out_matrix,
            "-o",
            out_profile,
            "--perGroup",
            "--regionsLabel",
            *region_labels,
            "--colors",
            *colours,
            *point_labels,
            "--numPlotsPerRow",
            str(len(beds)),
            "--plotHeight",
            "7",
            "--plotWidth",
            "8",
            "--yAxisLabel",
            y_label,
            *extra_plot.split(),
        ]
    )
    run(
        [
            "plotHeatmap",
            "-m",
            out_matrix,
            "-o",
            out_heatmap,
            "--regionsLabel",
            *heatmap_labels,
            *point_labels,
            "--colorMap",
            "Reds",
            "--sortRegions",
            "descend",
            # The profile is in its own plot
            "--whatToShow",
            "heatmap and colorbar",
            "--heatmapHeight",
            "14",
            "--heatmapWidth",
            "4",
            "--xAxisLabel",
            "",
            "--yAxisLabel",
            y_label,
        ]
    )
subprocess.run(["rm", "-r", tmp], check=True)
