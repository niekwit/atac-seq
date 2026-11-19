"""
Profile and heatmap of the signal of the control and the treatment of a
contrast over groups of regions (deepTools computeMatrix, plotProfile and
plotHeatmap; one panel per group, one line per condition):

- ATAC-seq Tn5 insertions (TE signal tracks, CPM or DESeq2-scaled) over the
  regions where an epigenomic mark (H3K9me3, DNA methylation) decreases,
  increases or does not change
  (scale-regions: regions scaled to body_length bp)
- signal of the mark (e.g. H3K9me3 coverage, CpG methylation) around the
  centre of the more, less and unchanged accessible peaks (reference-point)

Bins without signal are 0 for coverage signal (missing_as_zero), and missing
(left out of the mean) for methylation, which is only defined at CpGs.

Regions come from a table with chr, start and end columns, grouped by
group_column. If the chromosomes of the BigWig files are UCSC names (chr1),
the Ensembl names of the regions (1, MT) are converted. Writes a PDF with a
message if no group has regions.
"""

import gzip
import json
import logging
import os
import subprocess
import tempfile

import numpy as np
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
missing_as_zero = snakemake.params.get("missing_as_zero", True)
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
# Peaks not tested by DESeq2 (padj NA) are left out, as in the other analyses;
# otherwise they would count as unchanged
if "padj" in table.columns:
    untested = table["padj"].isna()
    logging.info(f"{untested.sum()} regions not tested by DESeq2 left out")
    table = table[~untested]

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

# Removed at exit, also after an error
_tmp = tempfile.TemporaryDirectory(dir=os.path.dirname(out_matrix))
tmp = _tmp.name
# BED files of the groups, and the label of each (by file name)
beds, names = [], {}
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
    names[os.path.basename(path)] = label


def run(command):
    logging.info(" ".join(command))
    with open(log, "a") as fh:
        subprocess.run(command, check=True, stdout=fh, stderr=fh)


def no_data(message):
    """PDFs with a message instead of the profile and the heatmap"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logging.warning(f"{message}: no profile")
    for path in (out_profile, out_heatmap):
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.axis("off")
        ax.text(0.5, 0.5, message, ha="center")
        fig.savefig(path)
        plt.close(fig)


def has_signal(matrix):
    """
    Whether the computeMatrix matrix has a finite value (without
    --missingDataAsZero, e.g. methylation without CpGs in any window, all
    bins are nan, which plotHeatmap does not plot)
    """
    try:
        values = pd.read_csv(matrix, sep="\t", header=None, skiprows=1)
    except pd.errors.EmptyDataError:
        # No regions left (e.g. --skipZeros in extra: computematrix)
        return False
    values = values.iloc[:, 6:]
    values = values.apply(pd.to_numeric, errors="coerce")
    return bool(np.isfinite(values.to_numpy(dtype=float)).any())


if not beds:
    no_data(f"No {region_name}")
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
            *(["--missingDataAsZero"] if missing_as_zero else []),
            "-p",
            str(snakemake.threads),
            "-o",
            out_matrix,
            *extra_matrix.split(),
        ]
    )
    if not has_signal(out_matrix):
        no_data(f"No signal over the {region_name}")
    else:
        # Regions in the matrix per group: fewer than in the BED files with, e.g.,
        # --skipZeros in extra: computematrix; groups without regions are dropped
        with gzip.open(out_matrix, "rt") as fh:
            header = json.loads(fh.readline()[1:])
        boundaries = header["group_boundaries"]
        region_labels, heatmap_labels = [], []
        for i, group in enumerate(header["group_labels"]):
            # By file name; by position if deepTools renamed the groups (a single
            # group is "genes") and no group was dropped
            if group in names:
                label = names[group]
            elif len(header["group_labels"]) == len(beds):
                label = list(names.values())[i]
            else:
                label = group
            n = boundaries[i + 1] - boundaries[i]
            logging.info(f"{label}: {n} regions in the matrix")
            region_labels.append(f"{label} (n = {n})")
            # Heatmap labels without the numbers, which are in the profile
            heatmap_labels.append(label)
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
                str(len(region_labels)),
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
