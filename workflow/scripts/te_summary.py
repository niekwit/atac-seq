"""
Summarises the Tn5 insertion counts in repeats per subfamily, family and
class. For each sample:

- count: insertions, multimappers weighted 1/NH (featureCounts -M --fraction)
- count_unique: insertions of uniquely mapping read pairs only
- frac_unique: count_unique / count; low values mean that the signal of a
  repeat comes mostly from multimappers, and that it would be largely missed
  by the ENCODE (unique/low multimapping) alignments
- CPM: count per million insertions (reads of kept read pairs)
- log2_enrichment: log2 of the fraction of insertions in the repeat divided
  by the fraction of the genome it covers (pseudocount 1); > 0 means more
  accessible than the genome average

Counts are fractional; round them for count based tests (e.g. DESeq2).
"""

import logging

import numpy as np
import pandas as pd

# Load Snakemake variables
samples = snakemake.params["samples"]
filter_chrs = set(snakemake.params["filter_chrs"])
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def read_counts(f):
    """featureCounts table as subfamily x sample counts"""
    header = pd.read_csv(f, sep="\t", comment="#", nrows=0).columns
    bams = header[6:]
    df = pd.read_csv(
        f, sep="\t", comment="#", usecols=["Geneid", *bams], index_col="Geneid"
    )
    # results/te/{sample}.filtered.bam -> sample
    df.columns = [b.split("/")[-1].removesuffix(".filtered.bam") for b in df.columns]
    return df[samples]


subfamilies = pd.read_csv(snakemake.input["subfamilies"], sep="\t")
multi = read_counts(snakemake.input["multi"])
unique = read_counts(snakemake.input["unique"])

chroms = pd.read_csv(
    snakemake.input["chrom_sizes"], sep="\t", header=None, dtype={0: str}
)
genome_bp = chroms.loc[~chroms[0].isin(filter_chrs), 1].sum()

stats = pd.concat(
    [pd.read_csv(f, sep="\t") for f in snakemake.input["stats"]], ignore_index=True
)
stats.insert(0, "sample", samples)
stats = stats.set_index("sample")
insertions = stats["kept_reads"]

# Subfamily table with counts
sub = subfamilies.set_index("subfamily")
counts = {
    "count": multi.reindex(sub.index, fill_value=0),
    "count_unique": unique.reindex(sub.index, fill_value=0),
}

LEVELS = {
    "subfamily": ["subfamily", "family", "class"],
    "family": ["family", "class"],
    "class": ["class"],
}
for level, columns in LEVELS.items():
    annot = sub.reset_index()
    group = annot.groupby(columns, sort=True)
    out = group.agg(n_copies=("n_copies", "sum"), bp=("bp", "sum"))
    if level != "subfamily":
        out.insert(0, "n_subfamilies", group.size())
    out["frac_genome"] = out["bp"] / genome_bp
    keys = pd.MultiIndex.from_frame(annot[columns])
    for sample in samples:
        c = counts["count"][sample].set_axis(keys).groupby(level=columns).sum()
        u = counts["count_unique"][sample].set_axis(keys).groupby(level=columns).sum()
        out[f"{sample}_count"] = c
        out[f"{sample}_count_unique"] = u
        out[f"{sample}_frac_unique"] = (u / c).where(c > 0)
        out[f"{sample}_CPM"] = c / insertions[sample] * 1e6
        out[f"{sample}_log2_enrichment"] = np.log2(
            (c + 1) / insertions[sample] / out["frac_genome"]
        )
    out = out.reset_index()
    out.to_csv(
        snakemake.output["repeat_class" if level == "class" else level],
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    logging.info(f"{level}: {len(out)} rows")

# Library statistics
for sample in samples:
    stats.loc[sample, "insertions_in_repeats"] = multi[sample].sum()
    stats.loc[sample, "insertions_in_repeats_unique"] = unique[sample].sum()
stats["frac_insertions_in_repeats"] = stats["insertions_in_repeats"] / insertions
stats["frac_genome_in_repeats"] = subfamilies["bp"].sum() / genome_bp
stats["frac_duplicate_pairs"] = stats["duplicate_pairs"] / stats["proper_pairs"]
stats["frac_multimapping_pairs"] = stats["multimapping_pairs"] / stats["kept_pairs"]
stats.to_csv(snakemake.output["stats"], sep="\t", float_format="%.6g")
logging.info("Done!")
