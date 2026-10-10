"""
Converts the UCSC RepeatMasker annotation (rmsk.txt.gz) to:

1. A SAF file for featureCounts with one feature per repeat copy, named
   after its subfamily (repName), so that featureCounts sums the copies of
   each subfamily
2. A table of the subfamilies with their family, class, number of copies and
   total length (bp)

Chromosome names are converted to Ensembl names (chr1 -> 1, chrM -> MT).
Repeats on chromosomes that are not in the genome (e.g. UCSC alt/random
contigs) or on filter_chrs (e.g. MT) are removed. If classes is not empty,
only repeats of these RepeatMasker classes are kept. Uncertain class
assignments (e.g. "LTR?") count as their class; families keep RepeatMasker's
names.
"""

import logging

import pandas as pd

# Load Snakemake variables
rmsk = snakemake.input["rmsk"]
chrom_sizes = snakemake.input["chrom_sizes"]
saf_out = snakemake.output["saf"]
subfamilies_out = snakemake.output["subfamilies"]
classes = snakemake.params["classes"]
filter_chrs = snakemake.params["filter_chrs"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# Columns of the UCSC rmsk table (https://genome.ucsc.edu/cgi-bin/hgTables)
RMSK_COLUMNS = {
    5: "chrom",
    6: "start",  # 0-based
    7: "end",
    9: "strand",
    10: "subfamily",  # repName
    11: "class",  # repClass
    12: "family",  # repFamily
}

df = pd.read_csv(
    rmsk,
    sep="\t",
    header=None,
    usecols=list(RMSK_COLUMNS),
    dtype={5: str, 10: str, 11: str, 12: str},
).rename(columns=RMSK_COLUMNS)
logging.info(f"{len(df)} repeats in {rmsk}")

# Ensembl chromosome names
df["chrom"] = df["chrom"].str.removeprefix("chr").replace({"M": "MT"})
chroms = pd.read_csv(chrom_sizes, sep="\t", header=None, dtype={0: str})[0]
keep_chroms = set(chroms) - set(filter_chrs)
df = df[df["chrom"].isin(keep_chroms)]
logging.info(f"{len(df)} repeats on the chromosomes of the genome")

# Uncertain classifications (e.g. "LTR?") count as their class
df["class"] = df["class"].str.removesuffix("?")
if classes:
    df = df[df["class"].isin(classes)]
    logging.info(f"{len(df)} repeats of classes {', '.join(classes)}")

if df.empty:
    raise ValueError("No repeats left; check genome and te:classes in the config")

# SAF (1-based, inclusive)
saf = pd.DataFrame(
    {
        "GeneID": df["subfamily"],
        "Chr": df["chrom"],
        "Start": df["start"] + 1,
        "End": df["end"],
        "Strand": df["strand"],
    }
)
saf.to_csv(saf_out, sep="\t", index=False)

# Subfamilies; a few subfamily names occur in more than one family or class,
# these get their most frequent one
df["length"] = df["end"] - df["start"]
subfamilies = (
    df.groupby("subfamily")
    .agg(
        family=("family", lambda x: x.mode().iloc[0]),
        repeat_class=("class", lambda x: x.mode().iloc[0]),
        n_copies=("length", "size"),
        bp=("length", "sum"),
    )
    .rename(columns={"repeat_class": "class"})
    .reset_index()
)
subfamilies.to_csv(subfamilies_out, sep="\t", index=False)
logging.info(f"{len(subfamilies)} subfamilies")
logging.info("Done!")
