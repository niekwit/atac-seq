"""
Prepares the differential H3K9me3 regions of a contrast (h3k9me3: contrasts)
for the integration with the differential accessibility.

- DiffBind table (diffbind: true): every site with its log2 fold change and
  FDR; sites with FDR below h3k9me3: fdr are lost (log2 fold change < 0) or
  gained (> 0) in the treatment, the others unchanged. DiffBind reports
  1-based starts (GRanges), which are converted to BED (0-based) starts.
- BED file (diffbind: false): regions that change in one direction
  (direction: lost or gained), e.g. MACS2 peaks of the control called with
  the treatment as control sample. Overlapping regions (e.g. consensus peaks
  of replicates) are merged; there is no log2 fold change or FDR.

Chromosome names are converted to Ensembl names (chr1 -> 1, chrM -> MT), and
regions on chromosomes that are not in the genome are removed.

Writes regions.tsv (region_id, chr, start, end, width, log2FoldChange, FDR,
change) and BED files of the lost and gained regions (empty if none).
"""

import logging
import os
import subprocess
import tempfile

import pandas as pd

regions_file = snakemake.input["regions"]
chrom_sizes = snakemake.input["chrom_sizes"]
diffbind = snakemake.params["diffbind"]
direction = snakemake.params["direction"]
columns = snakemake.params["columns"]
flip_fold = snakemake.params["flip_fold"]
fdr = float(snakemake.params["fdr"])
out_table = snakemake.output["table"]
out_beds = {"lost": snakemake.output["lost"], "gained": snakemake.output["gained"]}
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

chroms = set(pd.read_csv(chrom_sizes, sep="\t", header=None, dtype={0: str})[0])


def ensembl_names(chrs):
    """Ensembl chromosome names (chr1 -> 1, chrM -> MT)"""
    return chrs.astype(str).str.removeprefix("chr").replace({"M": "MT"})


if diffbind:
    sep = "," if regions_file.lower().endswith(".csv") else "\t"
    table = pd.read_csv(regions_file, sep=sep)
    missing = [c for c in columns.values() if c not in table.columns]
    if missing:
        raise ValueError(
            f"Columns {missing} (h3k9me3: columns) not in {regions_file} "
            f"(columns: {', '.join(table.columns)})"
        )
    regions = pd.DataFrame(
        {
            "chr": ensembl_names(table[columns["chr"]]),
            # GRanges (1-based) to BED (0-based) start
            "start": table[columns["start"]].astype(int) - 1,
            "end": table[columns["end"]].astype(int),
            "log2FoldChange": table[columns["log2fc"]].astype(float),
            "FDR": table[columns["fdr"]].astype(float),
        }
    )
    if flip_fold:
        regions["log2FoldChange"] = -regions["log2FoldChange"]
    significant = regions["FDR"] < fdr
    regions["change"] = "unchanged"
    regions.loc[significant & (regions["log2FoldChange"] < 0), "change"] = "lost"
    regions.loc[significant & (regions["log2FoldChange"] > 0), "change"] = "gained"
    logging.info(f"{len(regions)} DiffBind sites in {regions_file}")
else:
    bed = pd.read_csv(
        regions_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chr", "start", "end"],
        dtype={"chr": str},
        comment="#",
    )
    bed = bed[~bed["chr"].str.startswith(("track", "browser"))]
    bed["chr"] = ensembl_names(bed["chr"])
    bed = bed.astype({"start": int, "end": int})
    logging.info(f"{len(bed)} regions in {regions_file}")
    # Merge overlapping regions (>= 1 bp; adjacent regions stay separate)
    with tempfile.NamedTemporaryFile(
        "w", suffix=".bed", dir=os.path.dirname(out_table), delete=False
    ) as fh:
        bed.sort_values(["chr", "start", "end"]).to_csv(
            fh, sep="\t", header=False, index=False
        )
        tmp = fh.name
    merged = subprocess.run(
        ["bedtools", "merge", "-d", "-1", "-i", tmp],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    os.remove(tmp)
    regions = pd.DataFrame(
        [line.split("\t")[:3] for line in merged.splitlines()],
        columns=["chr", "start", "end"],
    ).astype({"start": int, "end": int})
    logging.info(f"{len(regions)} regions after merging overlapping regions")
    regions["log2FoldChange"] = float("nan")
    regions["FDR"] = float("nan")
    regions["change"] = direction

off_genome = ~regions["chr"].isin(chroms)
if off_genome.any():
    logging.warning(
        f"{off_genome.sum()} regions on chromosomes not in the genome removed: "
        f"{', '.join(sorted(regions.loc[off_genome, 'chr'].unique()))}"
    )
regions = regions[~off_genome].sort_values(["chr", "start", "end"])
if regions.empty:
    raise ValueError(f"No H3K9me3 regions of {regions_file} on the genome")
regions.insert(0, "region_id", [f"h3k9me3_{i}" for i in range(1, len(regions) + 1)])
regions.insert(4, "width", regions["end"] - regions["start"])
regions.to_csv(out_table, sep="\t", index=False)

for change, path in out_beds.items():
    selected = regions[regions["change"] == change]
    selected[["chr", "start", "end", "region_id"]].to_csv(
        path, sep="\t", header=False, index=False
    )
    logging.info(
        f"{change}: {len(selected)} regions, {selected['width'].sum()} bp "
        f"(median width {selected['width'].median():.0f} bp)"
        if len(selected)
        else f"{change}: no regions"
    )
logging.info(f"unchanged: {(regions['change'] == 'unchanged').sum()} regions")
