"""
Repeat copies in the differential H3K9me3 regions of a contrast (TE module).

- copies.tsv: per subfamily (resources/te/repeats.saf, all RepeatMasker
  copies of the counted classes), the number of copies and the copies with
  at least half their length in regions that lose and gain H3K9me3
  (bedtools intersect -f 0.5), and the fraction of the genome covered by
  those regions (the expected fraction of copies, if copies were placed at
  random).
- full_length.tsv: each full-length copy of the te: profile: subfamilies
  entries (resources/te/full_length_copies.bed) with the H3K9me3 change of
  the region it overlaps most (lost, gained, unchanged or outside) and its
  accessibility: more or less accessible if a differential TE peak overlaps
  it (more accessible if both), unchanged if only unchanged TE peaks do, no
  peak otherwise. TE peaks not tested by DESeq2 are left out.
"""

import io
import logging
import os
import subprocess
import tempfile

import pandas as pd

saf_file = snakemake.input["saf"]
regions_file = snakemake.input["regions"]
full_length_file = snakemake.input["full_length"]
peaks_file = snakemake.input["peaks"]
chrom_sizes = snakemake.input["chrom_sizes"]
out_copies = snakemake.output["copies"]
out_full_length = snakemake.output["full_length"]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log)],
    force=True,
)

tmp = tempfile.mkdtemp(dir=os.path.dirname(out_copies))


def write_bed(frame, name):
    path = os.path.join(tmp, f"{name}.bed")
    frame.sort_values(["chr", "start"]).to_csv(
        path, sep="\t", header=False, index=False
    )
    return path


def intersect(a, b, *options):
    """bedtools intersect -wo of two BED files, as a data frame (no header)"""
    result = subprocess.run(
        ["bedtools", "intersect", "-wo", *options, "-a", a, "-b", b],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if not result:
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(result), sep="\t", header=None)


regions = pd.read_csv(regions_file, sep="\t", dtype={"chr": str})
genome_bp = (
    pd.read_csv(chrom_sizes, sep="\t", header=None, usecols=[1])[1].astype(int).sum()
)

# Copies per subfamily in regions that lose or gain H3K9me3
# -----------------------------------------------------
# SAF: 1-based starts
saf = pd.read_csv(saf_file, sep="\t", dtype={"Chr": str})
copies = pd.DataFrame(
    {
        "chr": saf["Chr"],
        "start": saf["Start"] - 1,
        "end": saf["End"],
        "copy": range(len(saf)),
        "subfamily": saf["GeneID"],
    }
)
copies_bed = write_bed(copies, "copies")
counts = copies.groupby("subfamily").size().rename("n_copies").to_frame()
for change in ("lost", "gained"):
    selected = regions[regions["change"] == change]
    counts[f"n_copies_{change}"] = 0
    counts[f"genome_fraction_{change}"] = selected["width"].sum() / genome_bp
    if selected.empty:
        continue
    bed = write_bed(selected[["chr", "start", "end", "region_id"]], change)
    hits = intersect(copies_bed, bed, "-f", "0.5")
    if hits.empty:
        continue
    # Each copy once (a copy can overlap two adjacent regions)
    hits = hits.drop_duplicates(3)
    counts[f"n_copies_{change}"] = (
        hits.groupby(4).size().reindex(counts.index, fill_value=0).astype(int)
    )
    logging.info(
        f"{len(hits)} of {len(copies)} copies with >= 50% in regions with "
        f"H3K9me3 {change} ({selected['width'].sum() / genome_bp:.2%} of the genome)"
    )
counts.reset_index().to_csv(out_copies, sep="\t", index=False)

# Full-length copies: H3K9me3 change and accessibility
# -----------------------------------------------------
# deepTools groups: the copies of an entry are followed by "#<entry> (n=...)"
rows, pending = [], []
with open(full_length_file) as fh:
    for line in fh:
        if line.startswith("#"):
            group = line[1:].split("\t")[0].rsplit(" (n=", 1)[0].strip()
            rows += [(*row, group) for row in pending]
            pending = []
        elif line.strip():
            fields = line.rstrip("\n").split("\t")
            pending.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
# Copies after the last group line (none in files of te_full_length_copies)
rows += [(*row, row[3]) for row in pending]
full_length = pd.DataFrame(
    rows, columns=["chr", "start", "end", "subfamily", "group"]
).astype({"chr": str})
full_length["copy"] = range(len(full_length))
full_length["h3k9me3_change"] = "outside"
full_length["accessibility"] = "no peak"
if not full_length.empty:
    fl_bed = write_bed(full_length[["chr", "start", "end", "copy"]], "full_length")
    regions_bed = write_bed(
        regions[["chr", "start", "end", "region_id", "change"]], "regions"
    )
    hits = intersect(fl_bed, regions_bed)
    if not hits.empty:
        # Region with the largest overlap of each copy
        hits = hits.sort_values(9, ascending=False).drop_duplicates(3)
        full_length.loc[hits[3].values, "h3k9me3_change"] = hits[8].values

    peaks = pd.read_csv(peaks_file, sep="\t", dtype={"chr": str})
    peaks = peaks[peaks["padj"].notna()]
    peaks_bed = write_bed(peaks[["chr", "start", "end", "change"]], "peaks")
    hits = intersect(fl_bed, peaks_bed)
    if not hits.empty:
        for change in ("unchanged", "less accessible", "more accessible"):
            # Later changes take precedence
            copies_with = hits.loc[hits[7] == change, 3].unique()
            full_length.loc[copies_with, "accessibility"] = change
full_length = full_length.drop(columns="copy")
full_length.to_csv(out_full_length, sep="\t", index=False)
if not full_length.empty:
    logging.info(
        "Full-length copies by H3K9me3 change and accessibility:\n"
        + pd.crosstab(
            [full_length["group"], full_length["h3k9me3_change"]],
            full_length["accessibility"],
        ).to_string()
    )
subprocess.run(["rm", "-r", tmp], check=True)
