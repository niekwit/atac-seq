"""
Profile of the Tn5 insertions (CPM, mean of the replicates of each
condition) around the TSS of the up-, down- and unchanged genes of the
RNA-seq data of a contrast (deepTools computeMatrix reference-point and
plotProfile, one panel per group of genes, one line per condition).

The TSS of a gene is the 5' end of its GTF gene line (start on +, end on -).
The groups are the RNA-seq calls of rnaseq_prepare (rnaseq: fdr and
lfc_threshold); genes not tested by DESeq2 are left out. Writes a PDF with
a message if no group has genes.
"""

import gzip
import logging
import subprocess

import pandas as pd

rnaseq_file = snakemake.input["rnaseq"]
gtf_file = snakemake.input["gtf"]
bigwigs = list(snakemake.input["bigwigs"])
labels = list(snakemake.params["labels"])
colours = list(snakemake.params["colours"])
flank = int(snakemake.params["flank"])
extra_matrix = snakemake.params["extra_computematrix"]
extra_plot = snakemake.params["extra_plotprofile"]
out_beds = {
    "up": snakemake.output["up"],
    "down": snakemake.output["down"],
    "unchanged": snakemake.output["unchanged"],
}
out_matrix = snakemake.output["matrix"]
out_pdf = snakemake.output["pdf"]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# TSS of the genes of the GTF
# -----------------------------------------------------
gtf = pd.read_csv(
    gtf_file,
    sep="\t",
    comment="#",
    header=None,
    usecols=[0, 2, 3, 4, 6, 8],
    names=["chr", "feature", "start", "end", "strand", "attributes"],
    dtype={"chr": str},
)
gtf = gtf[gtf["feature"] == "gene"].copy()
gtf["gene_id"] = gtf["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
gtf = gtf.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
# BED: 0-based start of the 1 bp at the 5' end
gtf["tss"] = gtf["start"].where(gtf["strand"] == "+", gtf["end"]) - 1
tss = gtf.set_index("gene_id")[["chr", "tss", "strand"]]

# Genes of each RNA-seq group
# -----------------------------------------------------
rna = pd.read_csv(rnaseq_file, sep="\t", dtype={"gene_id": str})
rna = rna[(rna["type"] == "gene") & rna["gene_id"].isin(tss.index)]
counts = {}
for group, path in out_beds.items():
    ids = rna.loc[rna["rna_change"] == group, "gene_id"]
    bed = tss.loc[ids].reset_index()
    bed = pd.DataFrame(
        {
            "chr": bed["chr"],
            "start": bed["tss"],
            "end": bed["tss"] + 1,
            "name": bed["gene_id"],
            "score": 0,
            "strand": bed["strand"],
        }
    ).sort_values(["chr", "start"])
    bed.to_csv(path, sep="\t", header=False, index=False)
    counts[group] = len(bed)
    logging.info(f"{group}: {len(bed)} genes")

groups = [g for g in out_beds if counts[g] > 0]
if not groups:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    logging.warning("No up-, down- or unchanged genes: no profile")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.axis("off")
    ax.text(0.5, 0.5, "No genes tested in the RNA-seq data", ha="center")
    fig.savefig(out_pdf)
    with gzip.open(out_matrix, "wt"):
        pass
else:
    region_labels = [f"{g} (n = {counts[g]})" for g in groups]

    def run(command):
        logging.info(" ".join(command))
        with open(log, "a") as fh:
            subprocess.run(command, check=True, stdout=fh, stderr=fh)

    run(
        [
            "computeMatrix",
            "reference-point",
            "--referencePoint",
            "TSS",
            "-R",
            *[out_beds[g] for g in groups],
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
            out_pdf,
            "--perGroup",
            "--regionsLabel",
            *region_labels,
            "--colors",
            *colours,
            "--refPointLabel",
            "TSS",
            "--numPlotsPerRow",
            str(len(groups)),
            "--plotHeight",
            "7",
            "--plotWidth",
            "8",
            "--yAxisLabel",
            "Tn5 insertions per million",
            *extra_plot.split(),
        ]
    )
