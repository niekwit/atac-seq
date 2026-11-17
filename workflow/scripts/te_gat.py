"""
Tests whether the differential TE peaks of a contrast overlap regions of
other data sets (te: correlation: datasets, e.g. EM-seq DMRs) more than
expected, with GAT (Heger et al., 2013, doi:10.1093/bioinformatics/btt343).

- Segments: the peaks more accessible and less accessible in the treatment,
  as two segment tracks.
- Annotations: the regions of each data set, named after the data set.
  Chromosome names are converted to Ensembl names (chr1 -> 1, chrM -> MT),
  and regions on chromosomes that are not in the genome are removed.
- Workspace: the whole genome (all chromosomes of the genome), without
  filter_chrs (e.g. MT) and the blacklist, where no peaks can be called. GAT
  places the differential peaks at random within the workspace and compares
  the observed overlap (bp) with the overlap of the random placements.
  q-values are over all tests of the contrast.

For the Venn diagrams, the peaks and the regions of each data set are first
merged within each set (only intervals that overlap by >= 1 bp, so adjacent
tiles, e.g. 1 kb DMR tiles, stay separate). Counted are the merged peaks, the
merged regions, the intersections (overlapping peak-region pairs, as
bedtools intersect reports them), the peaks that overlap at least one region
and the regions that overlap at least one peak.
"""

import logging
import os
import subprocess
import tempfile

import pandas as pd

# Load Snakemake variables
more_bed = snakemake.input["more"]
less_bed = snakemake.input["less"]
chrom_sizes = snakemake.input["chrom_sizes"]
blacklist = snakemake.input["blacklist"]
filter_chrs = set(snakemake.params["filter_chrs"])
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


def merge_within(bed, name):
    """Merges the intervals of one set that overlap (>= 1 bp); writes a BED3"""
    path = os.path.join(tmp, f"{name}.bed")
    bed.sort_values(["chr", "start", "end"]).to_csv(
        path, sep="\t", header=False, index=False
    )
    command = ["bedtools", "merge", "-d", "-1", "-i", path]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    merged_path = os.path.join(tmp, f"{name}.merged.bed")
    with open(merged_path, "w") as fh:
        fh.write(result.stdout)
    return merged_path, len(result.stdout.splitlines())


def count_lines(command):
    """Number of output lines of a bedtools command"""
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return len(result.stdout.splitlines())


def count_overlaps(peaks, regions):
    """
    Venn counts of two merged BED files: intersections (overlapping pairs),
    peaks that overlap a region and regions that overlap a peak
    """
    intersect = ["bedtools", "intersect"]
    return {
        "n_intersections": count_lines(intersect + ["-a", peaks, "-b", regions]),
        "n_peaks_overlap": count_lines(intersect + ["-u", "-a", peaks, "-b", regions]),
        "n_regions_overlap": count_lines(
            intersect + ["-u", "-a", regions, "-b", peaks]
        ),
    }


# Workspace: the whole genome without filter_chrs and the blacklist (BED3)
genome = pd.read_csv(
    chrom_sizes, sep="\t", header=None, usecols=[0, 1], names=["chr", "end"]
).astype({"chr": str})
genome = genome[~genome["chr"].isin(filter_chrs)].assign(start=0)
genome_file = os.path.join(tmp, "genome.bed")
genome.sort_values("chr")[["chr", "start", "end"]].to_csv(
    genome_file, sep="\t", header=False, index=False
)
workspace_file = os.path.join(tmp, "workspace.bed")
if os.path.getsize(blacklist) > 0:
    blacklist_file, _ = merge_within(read_bed(blacklist), "blacklist")
    with open(workspace_file, "w") as fh:
        subprocess.run(
            ["bedtools", "subtract", "-a", genome_file, "-b", blacklist_file],
            check=True,
            stdout=fh,
        )
else:
    os.replace(genome_file, workspace_file)
logging.info(
    f"Workspace: {len(genome)} chromosomes without filter_chrs "
    f"({', '.join(sorted(filter_chrs))}) and the blacklist"
)

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
    region_file, n_regions = merge_within(
        bed[["chr", "start", "end"]], f"dataset_{len(annotations)}"
    )
    bed["name"] = name
    annotations.append(bed)
    for track, peaks in TRACKS.items():
        peak_file, n_peaks = merge_within(
            read_bed(peaks), track.replace(" ", "_") + "_peaks"
        )
        counts.append(
            {
                "track": track,
                "annotation": name,
                "n_peaks": n_peaks,
                "n_regions": n_regions,
                **count_overlaps(peak_file, region_file),
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
