"""
Annotates the differential accessibility of the peaks of a contrast with the
differential H3K9me3 regions of the same contrast (h3k9me3_prepare).

Each peak tested by DESeq2 (padj not NA) gets the H3K9me3 region it overlaps
most (bedtools intersect -wo, >= 1 bp), with the change of H3K9me3 in that
region (lost, gained or unchanged), or "outside" if it overlaps no region.
Peaks not tested by DESeq2 are left out, as they would count as unchanged.
"""

import io
import logging
import os
import subprocess
import tempfile

import pandas as pd

peaks_file = snakemake.input["atac"]
regions_file = snakemake.input["regions"]
output = snakemake.output["peaks"]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

atac_columns = [
    "peak_id",
    "chr",
    "start",
    "end",
    "baseMean",
    "log2FoldChange",
    "log2FoldChange_shrunken",
    "padj",
    "change",
]
peaks = pd.read_csv(peaks_file, sep="\t", dtype={"chr": str})[atac_columns]
untested = peaks["padj"].isna()
logging.info(
    f"{len(peaks)} peaks, {untested.sum()} not tested by DESeq2 (padj NA) left out"
)
peaks = peaks[~untested]
regions = pd.read_csv(regions_file, sep="\t", dtype={"chr": str})

# Removed at exit, also after an error
_tmp = tempfile.TemporaryDirectory(dir=os.path.dirname(output))
tmp = _tmp.name
peaks_bed = os.path.join(tmp, "peaks.bed")
regions_bed = os.path.join(tmp, "regions.bed")
peaks[["chr", "start", "end", "peak_id"]].to_csv(
    peaks_bed, sep="\t", header=False, index=False
)
regions[["chr", "start", "end", "region_id"]].to_csv(
    regions_bed, sep="\t", header=False, index=False
)
result = subprocess.run(
    ["bedtools", "intersect", "-wo", "-a", peaks_bed, "-b", regions_bed],
    check=True,
    capture_output=True,
    text=True,
).stdout

overlaps = pd.read_csv(
    io.StringIO(result),
    sep="\t",
    header=None,
    usecols=[3, 7, 8],
    names=["peak_id", "region_id", "overlap_bp"],
)
# Region with the largest overlap of each peak
overlaps = overlaps.sort_values("overlap_bp", ascending=False).drop_duplicates(
    "peak_id"
)
overlaps = overlaps.merge(
    regions[["region_id", "width", "log2FoldChange", "FDR", "change"]].rename(
        columns={
            "region_id": "h3k9me3_region",
            "width": "h3k9me3_width",
            "log2FoldChange": "h3k9me3_log2FoldChange",
            "FDR": "h3k9me3_FDR",
            "change": "h3k9me3_change",
        }
    ),
    left_on="region_id",
    right_on="h3k9me3_region",
).drop(columns="region_id")
overlaps = overlaps.rename(columns={"overlap_bp": "h3k9me3_overlap_bp"})

annotated = peaks.merge(overlaps, on="peak_id", how="left")
annotated["h3k9me3_change"] = annotated["h3k9me3_change"].fillna("outside")
annotated["h3k9me3_overlap_bp"] = annotated["h3k9me3_overlap_bp"].fillna(0).astype(int)
annotated.to_csv(output, sep="\t", index=False, na_rep="NA")
logging.info(
    "Peaks by accessibility and H3K9me3 change:\n"
    + pd.crosstab(annotated["change"], annotated["h3k9me3_change"]).to_string()
)
