"""
Tests whether the differential TE peaks of a contrast overlap regions of
other data sets (te: correlation: datasets, e.g. EM-seq DMRs) more than
expected, with GAT (Heger et al., 2013, doi:10.1093/bioinformatics/btt343).

- Segments: the peaks more accessible and less accessible in the treatment,
  as two segment tracks.
- Annotations: the regions of each data set, named after the data set.
  Chromosome names are converted to Ensembl names (chr1 -> 1, chrM -> MT),
  and regions on chromosomes that are not in the genome are removed.
- Workspace: all TE consensus peaks, so the differential peaks are compared
  with accessible regions in general. GAT places the differential peaks at
  random within the consensus peaks and compares the observed overlap (bp)
  with the overlap of the random placements. q-values are over all tests
  of the contrast.

Also counts the peaks that overlap at least one region of a data set and the
regions that overlap at least one peak, and, for the Venn diagrams, merged
regions as in ChIPpeakAnno (Zhu et al., 2010, doi:10.1186/1471-2105-11-237):
peaks and regions that overlap (>= 1 bp) are merged, and the merged regions
with only peaks, only data set regions, or both are counted. One measure for
both sets, so a peak overlapping several regions counts once.
"""

import logging
import os
import subprocess
import tempfile

import pandas as pd

# Load Snakemake variables
more_bed = snakemake.input["more"]
less_bed = snakemake.input["less"]
workspace = snakemake.input["workspace"]
chrom_sizes = snakemake.input["chrom_sizes"]
datasets = snakemake.params["datasets"]
n_samples = snakemake.params["samples"]
output = snakemake.output["table"]
threads = snakemake.threads
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

TRACKS = {"more accessible": more_bed, "less accessible": less_bed}
chroms = set(pd.read_csv(chrom_sizes, sep="\t", header=None, dtype={0: str})[0])
tmp = tempfile.mkdtemp(dir=os.path.dirname(output))


def read_bed(path):
    """First three columns of a BED file, with Ensembl chromosome names"""
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chr", "start", "end"],
        dtype={"chr": str},
        comment="#",
    )
    bed = bed[~bed["chr"].str.startswith(("track", "browser"))]
    bed["chr"] = bed["chr"].str.removeprefix("chr").replace({"M": "MT"})
    return bed[bed["chr"].isin(chroms)].astype({"start": int, "end": int})


def count_overlapping(a, b):
    """Number of intervals in BED file a that overlap an interval in b"""
    command = ["bedtools", "intersect", "-u", "-a", a, "-b", b]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return len(result.stdout.splitlines())


def count_merged(peaks, regions):
    """
    Merged regions (overlap >= 1 bp) with only peaks, only regions, or both
    """
    bed = pd.concat([peaks.assign(set="peak"), regions.assign(set="region")])
    bed = bed.sort_values(["chr", "start", "end"])
    bed_file = os.path.join(tmp, "merge.bed")
    bed.to_csv(bed_file, sep="\t", header=False, index=False)
    command = ["bedtools", "merge", "-d", "-1", "-c", "4", "-o", "distinct"]
    result = subprocess.run(
        command + ["-i", bed_file], check=True, capture_output=True, text=True
    )
    sets = pd.Series([line.split("\t")[3] for line in result.stdout.splitlines()])
    return {
        "n_merged_peaks_only": int((sets == "peak").sum()),
        "n_merged_regions_only": int((sets == "region").sum()),
        "n_merged_both": int((sets == "peak,region").sum()),
    }


# Workspace: the consensus peaks as one track (BED3)
workspace_file = os.path.join(tmp, "workspace.bed")
read_bed(workspace).to_csv(workspace_file, sep="\t", header=False, index=False)

# Segments: one track per direction
segments = []
for track, path in TRACKS.items():
    bed = read_bed(path)
    logging.info(f"{track}: {len(bed)} peaks")
    bed["name"] = track
    segments.append(bed)
segments = pd.concat(segments)
segments_file = os.path.join(tmp, "segments.bed")
segments.to_csv(segments_file, sep="\t", header=False, index=False)

# Annotations: one per data set
annotations = []
counts = []
for name, path in datasets.items():
    bed = read_bed(path)
    logging.info(f"{name} ({path}): {len(bed)} regions on the genome's chromosomes")
    if bed.empty:
        raise ValueError(
            f"No regions of {name} ({path}) on the chromosomes of the genome"
        )
    region_file = os.path.join(tmp, f"dataset_{len(annotations)}.bed")
    bed.to_csv(region_file, sep="\t", header=False, index=False)
    bed["name"] = name
    annotations.append(bed)
    for track, peaks in TRACKS.items():
        merged = count_merged(read_bed(peaks), bed[["chr", "start", "end"]])
        counts.append(
            {
                "track": track,
                "annotation": name,
                "n_peaks": sum(1 for line in open(peaks) if line.strip()),
                "n_peaks_overlap": count_overlapping(peaks, region_file),
                "n_regions": len(bed),
                "n_regions_overlap": count_overlapping(region_file, peaks),
                **merged,
            }
        )
annotations_file = os.path.join(tmp, "annotations.bed")
pd.concat(annotations).to_csv(annotations_file, sep="\t", header=False, index=False)
counts = pd.DataFrame(counts)

# GAT
gat_columns = [
    "observed",
    "expected",
    "CI95low",
    "CI95high",
    "fold",
    "l2fold",
    "pvalue",
    "qvalue",
]
if segments.empty:
    logging.warning("No differential peaks: GAT is not run")
    results = counts.assign(**{c: float("nan") for c in gat_columns})
else:
    gat_file = os.path.join(tmp, "gat.tsv")
    command = [
        "gat-run.py",
        "--with-segment-tracks",
        f"--segments={segments_file}",
        f"--annotations={annotations_file}",
        f"--workspace={workspace_file}",
        f"--num-samples={n_samples}",
        "--random-seed=42",
        f"--num-threads={threads}",
        f"--stdout={gat_file}",
        f"--log={log}.gat",
    ]
    logging.info(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True)
    gat = pd.read_csv(gat_file, sep="\t")
    results = counts.merge(
        gat[["track", "annotation"] + gat_columns],
        on=["track", "annotation"],
        how="left",
    )

results = results.rename(
    columns={
        "track": "peaks",
        "annotation": "dataset",
        "observed": "observed_overlap_bp",
        "expected": "expected_overlap_bp",
    }
)
results.to_csv(output, sep="\t", index=False)
subprocess.run(["rm", "-r", tmp], check=True)
logging.info("Done!")
