"""
Signal of all data sets of a contrast over the full-length copies of the
te: profile: subfamilies entries (resources/te/full_length_copies.bed, one
group per entry), for the integration of all data sets (TE module):
ATAC-seq Tn5 insertions (TE signal tracks, mean of the replicates), H3K9me3
and CpG methylation of the control and the treatment.

- matrix.gz: deepTools computeMatrix scale-regions matrix of the six BigWig
  files (copies scaled to body_length bp, 5' to 3', with flank bp on each
  side). Coverage signal (ATAC-seq, H3K9me3) is computed with
  --missingDataAsZero, methylation without (bins without a CpG are missing),
  so in two runs joined with computeMatrixOperations cbind. Within each group,
  the copies are sorted by the change in accessibility (below).
- copies.tsv: per copy, the mean of the body bins of each BigWig file, the
  log2 ratio of the ATAC-seq signal of the treatment and the control (with a
  pseudocount of 1% of the median signal of the copies), and, from the
  integration of each mark (full_length_copies.tsv), the change of the mark
  in the region the copy overlaps most and its accessibility (TE peaks)
- heatmap.pdf: plotHeatmap of the matrix (one colour scale per data set,
  shared by the control and the treatment, from 0 to the 98th percentile of
  the bins; methylation 0 to 100; copies in the order of the matrix, most
  gained accessibility first, 5' to 3' from left to right)
- profile.tsv: mean signal per bin of each group and BigWig file (bins
  without a value left out), for the profile plot

Copies on chromosomes missing from a BigWig file are left out, so that both
computeMatrix runs have the same regions.
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

bed_file = snakemake.input["bed"]
bigwigs = list(snakemake.input["bigwigs"])
mark_tables = {
    "h3k9me3": snakemake.input["h3k9me3"],
    "methylation": snakemake.input["methylation"],
}
# Per BigWig file (in the order of bigwigs): data set (label and column name:
# atac, h3k9me3, methylation), role (control, treatment) and condition,
# whether missing bins are 0, and the colour map of the heatmap
datasets = list(snakemake.params["datasets"])
columns = list(snakemake.params["columns"])
roles = list(snakemake.params["roles"])
conditions = list(snakemake.params["conditions"])
missing_as_zero = list(snakemake.params["missing_as_zero"])
colour_maps = list(snakemake.params["colour_maps"])
flank = int(snakemake.params["flank"])
body_length = int(snakemake.params["body_length"])
bin_size = int(snakemake.params["bin_size"])
out_matrix = snakemake.output["matrix"]
out_copies = snakemake.output["copies"]
out_heatmap = snakemake.output["heatmap"]
out_profile = snakemake.output["profile"]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# Removed at exit, also after an error
_tmp = tempfile.TemporaryDirectory(dir=os.path.dirname(out_matrix))
tmp = _tmp.name


def run(command):
    logging.info(" ".join(command))
    with open(log, "a") as fh:
        subprocess.run(command, check=True, stdout=fh, stderr=fh)


# Full-length copies: the copies of a group are followed by "#<entry> (n=...)"
# -----------------------------------------------------
rows, pending = [], []
for line in open(bed_file):
    if line.startswith("#"):
        group = line[1:].split("\t")[0].rsplit(" (n=", 1)[0].strip()
        rows += [(*row, group) for row in pending]
        pending = []
    elif line.strip():
        f = line.rstrip("\n").split("\t")
        pending.append((f[0], int(f[1]), int(f[2]), f[3], f[5]))
rows += [(*row, row[3]) for row in pending]
copies = pd.DataFrame(
    rows, columns=["chr", "start", "end", "subfamily", "strand", "group"]
).astype({"chr": str})
groups = list(dict.fromkeys(copies["group"]))

# Copies on chromosomes of all BigWig files (Ensembl or UCSC names)
for path in bigwigs:
    bw = pyBigWig.open(path)
    chroms = set(bw.chroms())
    bw.close()
    ucsc = ("chr" + copies["chr"]).replace({"chrMT": "chrM"})
    present = copies["chr"].isin(chroms) | ucsc.isin(chroms)
    if not present.all():
        logging.warning(
            f"{path}: {(~present).sum()} copies on chromosomes not in the file "
            "left out"
        )
    copies = copies[present]
copies = copies.reset_index(drop=True)
copies["copy_id"] = [f"copy_{i}" for i in range(len(copies))]
for group in groups:
    logging.info(f"{group}: {(copies['group'] == group).sum()} copies")
groups = [g for g in groups if (copies["group"] == g).any()]
if not groups:
    raise ValueError("No full-length copies on the chromosomes of the BigWigs")

# BED file per group (deepTools group label: the file name)
beds = []
for i, group in enumerate(groups):
    path = os.path.join(tmp, f"group{i}.bed")
    selected = copies[copies["group"] == group]
    selected.assign(score=0)[
        ["chr", "start", "end", "copy_id", "score", "strand"]
    ].to_csv(path, sep="\t", header=False, index=False)
    beds.append(path)

# Matrix: coverage and methylation signal, joined
# -----------------------------------------------------
labels = [f"{d} {c}" for d, c in zip(datasets, conditions)]
parts = []
for as_zero in (True, False):
    selected = [i for i, z in enumerate(missing_as_zero) if z == as_zero]
    if not selected:
        continue
    path = os.path.join(tmp, f"matrix_{int(as_zero)}.gz")
    run(
        [
            "computeMatrix",
            "scale-regions",
            "-R",
            *beds,
            "-S",
            *[bigwigs[i] for i in selected],
            "--samplesLabel",
            *[labels[i] for i in selected],
            "--regionBodyLength",
            str(body_length),
            "-b",
            str(flank),
            "-a",
            str(flank),
            "--binSize",
            str(bin_size),
            "--sortRegions",
            "keep",
            *(["--missingDataAsZero"] if as_zero else []),
            "-p",
            str(snakemake.threads),
            "-o",
            path,
        ]
    )
    parts.append((selected, path))
combined = os.path.join(tmp, "combined.gz")
if len(parts) > 1:
    run(
        [
            "computeMatrixOperations",
            "cbind",
            "-m",
            *[p for _, p in parts],
            "-o",
            combined,
        ]
    )
else:
    combined = parts[0][1]
order = [i for selected, _ in parts for i in selected]

with gzip.open(combined, "rt") as fh:
    header = json.loads(fh.readline()[1:])
    matrix = pd.read_csv(fh, sep="\t", header=None, dtype={0: str})
regions = matrix.iloc[:, :6].copy()
regions.columns = ["chr", "start", "end", "copy_id", "score", "strand"]
values = matrix.iloc[:, 6:].apply(pd.to_numeric, errors="coerce").to_numpy(float)
if len(regions) != len(copies):
    raise ValueError(
        f"{len(regions)} regions in the matrix, {len(copies)} copies: "
        "computeMatrix left out regions"
    )
bounds = header["sample_boundaries"]
n_flank = flank // bin_size
n_body = body_length // bin_size

# Mean of the body bins per copy and BigWig file
means = {}
for k, i in enumerate(order):
    block = values[:, bounds[k] : bounds[k + 1]]
    body = block[:, n_flank : n_flank + n_body]
    with np.errstate(all="ignore"):
        means[i] = np.nanmean(body, axis=1)
table = regions[["copy_id"]].copy()
for i in range(len(bigwigs)):
    table[f"{columns[i]}_{roles[i]}"] = means[i]
atac = table[["atac_control", "atac_treatment"]].to_numpy()
positive = atac[np.isfinite(atac) & (atac > 0)]
pseudocount = 0.01 * np.median(positive) if positive.size else 1.0
logging.info(f"ATAC-seq pseudocount: {pseudocount:g}")
table["atac_log2ratio"] = np.log2(
    (table["atac_treatment"] + pseudocount) / (table["atac_control"] + pseudocount)
)
copies = copies.merge(table, on="copy_id")

# Change of each mark and accessibility, from the integration with the mark
for mark, path in mark_tables.items():
    marked = pd.read_csv(path, sep="\t", dtype={"chr": str})
    keep = ["chr", "start", "end", "subfamily", "group", f"{mark}_change"]
    if mark == "h3k9me3":
        keep.append("accessibility")
    copies = copies.merge(
        marked[keep].drop_duplicates(["chr", "start", "end", "group"]),
        on=["chr", "start", "end", "subfamily", "group"],
        how="left",
    )

# Copies of each group sorted by the change in accessibility
# -----------------------------------------------------
rank = copies.sort_values(
    ["atac_log2ratio", "atac_treatment"], ascending=False
).reset_index(drop=True)
rank["position"] = range(len(rank))
position = dict(zip(rank["copy_id"], rank["position"]))
new_rows, new_bounds = [], [0]
group_bounds = header["group_boundaries"]
for g in range(len(group_bounds) - 1):
    idx = np.arange(group_bounds[g], group_bounds[g + 1])
    idx = idx[np.argsort([position[c] for c in regions["copy_id"].iloc[idx]])]
    new_rows.extend(idx)
    new_bounds.append(new_bounds[-1] + len(idx))
matrix = matrix.iloc[new_rows]
with gzip.open(out_matrix, "wt") as fh:
    fh.write("@" + json.dumps(header) + "\n")
    matrix.to_csv(fh, sep="\t", header=False, index=False, na_rep="nan")

copies["group"] = pd.Categorical(copies["group"], categories=groups, ordered=True)
copies["position"] = copies["copy_id"].map(position)
copies = copies.sort_values(["group", "position"]).drop(columns=["copy_id", "position"])
copies.to_csv(out_copies, sep="\t", index=False, na_rep="NA")
for mark in mark_tables:
    logging.info(
        f"Full-length copies by {mark} change:\n"
        + pd.crosstab(copies["group"], copies[f"{mark}_change"]).to_string()
    )

# Heatmap: colour scale per data set, from 0 to the 98th percentile of the
# body bins (methylation: 0 to 100)
# -----------------------------------------------------
z_min, z_max = [], []
for i in range(len(bigwigs)):
    pair = [j for j in range(len(bigwigs)) if datasets[j] == datasets[i]]
    if not missing_as_zero[i]:
        high = np.nanmax([np.nanmax(means[j]) for j in pair] + [0])
        z_min.append(0)
        z_max.append(100 if high > 1 else 1)
        continue
    k = [order.index(j) for j in pair]
    pooled = np.concatenate([values[:, bounds[m] : bounds[m + 1]].ravel() for m in k])
    pooled = pooled[np.isfinite(pooled)]
    high = np.percentile(pooled, 98) if pooled.size else 1
    z_min.append(0)
    z_max.append(float(high) if high > 0 else 1)
group_labels = [f"{g} (n={(copies['group'] == g).sum()})" for g in groups]
run(
    [
        "plotHeatmap",
        "-m",
        out_matrix,
        "-o",
        out_heatmap,
        "--sortRegions",
        "keep",
        "--colorMap",
        *[colour_maps[i] for i in order],
        "--zMin",
        *[f"{z_min[i]:g}" for i in order],
        "--zMax",
        *[f"{z_max[i]:g}" for i in order],
        "--missingDataColor",
        "#eeeeee",
        "--regionsLabel",
        *group_labels,
        # Copies 5' to 3' from left to right; the labels would overlap the
        # flank ticks
        "--startLabel",
        "",
        "--endLabel",
        "",
        "--whatToShow",
        "heatmap and colorbar",
        "--heatmapHeight",
        "22",
        "--heatmapWidth",
        "4",
        "--xAxisLabel",
        "",
    ]
)

# Mean profile per group and BigWig file
# -----------------------------------------------------
n_bins = 2 * n_flank + n_body
position = (np.arange(n_bins) + 0.5) * bin_size - flank
profile = []
for g, group in enumerate(groups):
    rows_g = values[group_bounds[g] : group_bounds[g + 1]]
    for k, i in enumerate(order):
        block = rows_g[:, bounds[k] : bounds[k + 1]]
        with np.errstate(all="ignore"):
            mean = np.nanmean(block, axis=0)
        profile.append(
            pd.DataFrame(
                {
                    "group": group,
                    "n_copies": len(rows_g),
                    "dataset": datasets[i],
                    "condition": conditions[i],
                    "position": position,
                    "signal": mean,
                }
            )
        )
pd.concat(profile).to_csv(out_profile, sep="\t", index=False, na_rep="NA")
logging.info("Done!")
