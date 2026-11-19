"""
Prepares the differential regions of an epigenomic mark of a contrast
({mark}: contrasts: H3K9me3, DNA methylation) for the integration with the
differential accessibility.

- Results table (h3k9me3: diffbind: true, methylation: methylkit: true): every
  tested region with its effect size (log2 fold change, or methylation
  difference in percentage points) and FDR (q-value); regions with FDR below
  {mark}: fdr and |effect| above min_effect decrease (lost, hypo) or increase
  (gained, hyper) in the treatment, the others are unchanged; rows without
  a q-value or effect are removed. DiffBind and
  methylKit report 1-based starts, which are converted to BED (0-based).
- BED file: regions that change in one direction (direction), e.g. MACS2
  peaks of the control called with the treatment as control sample, or the
  DMRs of one direction. Overlapping regions (e.g. consensus peaks of
  replicates) are merged; there is no effect size or FDR.

Chromosome names are converted to Ensembl names (chr1 -> 1, chrM -> MT), and
regions on chromosomes that are not in the genome, or outside their
chromosome, are removed.

Writes regions.tsv (region_id, chr, start, end, width, the effect, FDR,
change) and BED files of the regions where the mark decreases and increases
(empty if none).
"""

import logging
import os
import subprocess
import tempfile

import pandas as pd

regions_file = snakemake.input["regions"]
chrom_sizes = snakemake.input["chrom_sizes"]
name = snakemake.params["name"]
# Mark (h3k9me3, methylation), the prefix of the region IDs
prefix = snakemake.params["prefix"]
is_table = snakemake.params["table"]
direction = snakemake.params["direction"]
down = snakemake.params["down"]
up = snakemake.params["up"]
effect_name = snakemake.params["effect_name"]
columns = snakemake.params["columns"]
flip = snakemake.params["flip"]
fdr = float(snakemake.params["fdr"])
min_effect = float(snakemake.params["min_effect"])
out_table = snakemake.output["table"]
out_beds = {down: snakemake.output["decreased"], up: snakemake.output["increased"]}
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

lengths = (
    pd.read_csv(chrom_sizes, sep="\t", header=None, usecols=[0, 1], dtype={0: str})
    .set_index(0)[1]
    .astype(int)
)
chroms = set(lengths.index)


def ensembl_names(chrs):
    """Ensembl chromosome names (chr1 -> 1, chrM -> MT)"""
    return chrs.astype(str).str.removeprefix("chr").replace({"M": "MT"})


def on_genome(regions):
    """
    Removes regions on chromosomes not in the genome, outside their
    chromosome or with start >= end (0 <= start < end <= length); logs them
    """
    off_genome = ~regions["chr"].isin(chroms)
    if off_genome.any():
        logging.warning(
            f"{off_genome.sum()} regions on chromosomes not in the genome removed: "
            f"{', '.join(sorted(regions.loc[off_genome, 'chr'].unique()))}"
        )
    regions = regions[~off_genome]
    length = regions["chr"].map(lengths)
    invalid = (regions["start"] < 0) | (regions["start"] >= regions["end"])
    invalid |= regions["end"] > length
    if invalid.any():
        logging.warning(
            f"{invalid.sum()} regions outside their chromosome or with start >= "
            "end removed:\n"
            + regions.loc[invalid, ["chr", "start", "end"]]
            .head(20)
            .to_string(index=False)
        )
    return regions[~invalid]


if is_table:
    sep = "," if regions_file.lower().endswith(".csv") else "\t"
    table = pd.read_csv(regions_file, sep=sep)
    missing = [c for c in columns.values() if c not in table.columns]
    if missing:
        raise ValueError(
            f"Columns {missing} (columns of the config) not in {regions_file} "
            f"(columns: {', '.join(table.columns)})"
        )
    regions = pd.DataFrame(
        {
            "chr": ensembl_names(table[columns["chr"]]),
            # GRanges (1-based) to BED (0-based) start
            "start": table[columns["start"]].astype(int) - 1,
            "end": table[columns["end"]].astype(int),
            effect_name: table[columns["effect"]].astype(float),
            "FDR": table[columns["fdr"]].astype(float),
        }
    )
    # Rows without a test result (missing q-value or effect) cannot be called
    # changed or unchanged, so they are left out
    untested = regions[[effect_name, "FDR"]].isna().any(axis=1)
    if untested.any():
        logging.warning(
            f"{untested.sum()} regions without {columns['fdr']} or "
            f"{columns['effect']} removed"
        )
        regions = regions[~untested].copy()
    if regions.empty:
        raise ValueError(f"No regions with a test result in {regions_file}")
    if flip:
        regions[effect_name] = -regions[effect_name]
    significant = (regions["FDR"] < fdr) & (regions[effect_name].abs() > min_effect)
    regions["change"] = "unchanged"
    regions.loc[significant & (regions[effect_name] < 0), "change"] = down
    regions.loc[significant & (regions[effect_name] > 0), "change"] = up
    logging.info(f"{len(regions)} tested {name} regions in {regions_file}")
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
    # Before merging, so that an invalid record does not remove valid ones
    bed = on_genome(bed)
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
    regions[effect_name] = float("nan")
    regions["FDR"] = float("nan")
    regions["change"] = direction

regions = on_genome(regions).sort_values(["chr", "start", "end"])
if regions.empty:
    raise ValueError(f"No {name} regions of {regions_file} on the genome")
regions.insert(0, "region_id", [f"{prefix}_{i}" for i in range(1, len(regions) + 1)])
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
