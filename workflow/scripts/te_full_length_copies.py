"""
Writes the full-length copies of selected repeat subfamilies (e.g. young
LINE-1s) from the UCSC RepeatMasker annotation as a BED file for deepTools
computeMatrix, with one group of regions per entry of te: profile:
subfamilies. An entry is a subfamily (repName, e.g. L1MdTf_I) or, if no
subfamily has that name, the prefix of a group of subfamilies (e.g. L1MdTf for
L1MdTf_I, L1MdTf_II and L1MdTf_III: all repNames that start with the entry
followed by "_"), whose copies are then plotted together. "family:<name>" and
"class:<name>" select all subfamilies of a RepeatMasker family (repFamily,
e.g. family:L1) or class (repClass, e.g. class:LINE).

A copy is full-length if its alignment covers at least min_consensus_fraction
of the subfamily's consensus sequence (repStart, repEnd and repLeft of the
rmsk table). The consensus length is not a column of rmsk: each row gives it
as its last aligned consensus position plus the consensus bases after it, and
rows of one subfamily do not always agree (e.g. 6,579 bp for most L1MdTf_I
rows, but 3,914 to 8,042 bp for some). The consensus length of a subfamily is
therefore the most common value of its rows (the largest if tied), so that a
row reporting a shorter consensus does not pass as full-length. RepeatMasker
can split one element into several rows; these are not joined, so an element
counts only if one row covers enough of the consensus. Chromosome names are converted to Ensembl names and copies on
filter_chrs or on chromosomes not in the genome are removed, as in
prepare_te_annotation.py.

In the BED file, each group of regions is followed by a line "#<label>",
the group label (te_profile_matrix.py splits the groups into one BED file per
label for computeMatrix). The TSV file has all copies of each entry
(full-length or not) with their length and covered fraction, for the length QC
plot (te_profile_lengths.R).
"""

import logging

import pandas as pd

# Load Snakemake variables
rmsk = snakemake.input["rmsk"]
chrom_sizes = snakemake.input["chrom_sizes"]
bed_out = snakemake.output["bed"]
copies_out = snakemake.output["copies"]
subfamilies = snakemake.params["subfamilies"]
min_fraction = snakemake.params["min_consensus_fraction"]
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
    13: "rep_start",
    14: "rep_end",
    15: "rep_left",
}

df = pd.read_csv(
    rmsk,
    sep="\t",
    header=None,
    usecols=list(RMSK_COLUMNS),
    dtype={5: str, 10: str, 11: str, 12: str},
).rename(columns=RMSK_COLUMNS)


def members(entry):
    """
    Subfamilies of an entry: family:<name> or class:<name>, the subfamily
    itself, or the subfamilies entry_* if there is no subfamily entry
    """
    for level in ("family", "class"):
        if entry.startswith(f"{level}:"):
            name = entry.removeprefix(f"{level}:")
            return sorted(set(df.loc[df[level] == name, "subfamily"]))
    names = set(df["subfamily"])
    if entry in names:
        return [entry]
    return sorted(name for name in names if name.startswith(f"{entry}_"))


groups = {entry: members(entry) for entry in subfamilies}
for entry, group in groups.items():
    if not group:
        logging.warning(f"{entry}: no matching subfamilies")
    elif group != [entry]:
        logging.info(f"{entry}: {len(group)} subfamilies: {', '.join(group)}")
df = df[df["subfamily"].isin({name for group in groups.values() for name in group})]

# Consensus positions covered by each copy. On the + strand the alignment
# runs from repStart to repEnd and repLeft is minus the consensus bases after
# it; on the - strand repStart is minus the bases after repEnd and repLeft is
# the start of the alignment
plus = df["strand"] == "+"
first = df["rep_start"].where(plus, df["rep_left"])
after = -df["rep_left"].where(plus, df["rep_start"])
row_consensus_length = df["rep_end"] + after


def most_common(lengths):
    """
    Most common consensus length of the rows of a subfamily (the largest if
    tied)
    """
    counts = lengths.value_counts()
    return counts[counts == counts.max()].index.max()


# Consensus length per subfamily, from all its rows (before the chromosome
# filter)
consensus_lengths = row_consensus_length.groupby(df["subfamily"]).agg(most_common)
for name, length in consensus_lengths.items():
    rows = df["subfamily"] == name
    other = (row_consensus_length[rows] != length).sum()
    logging.info(
        f"{name}: consensus length {length} bp ({other} of {rows.sum()} rows "
        "report another length)"
    )
df["consensus_length"] = df["subfamily"].map(consensus_lengths)
df["fraction"] = (df["rep_end"] - first + 1) / df["consensus_length"]

# Ensembl chromosome names
df["chrom"] = df["chrom"].str.removeprefix("chr").replace({"M": "MT"})
chroms = pd.read_csv(chrom_sizes, sep="\t", header=None, dtype={0: str})[0]
df = df[df["chrom"].isin(set(chroms) - set(filter_chrs))]

df["full_length"] = df["fraction"] >= min_fraction
full_length = df[df["full_length"]]

# All copies of each entry, for the length QC plot
copies = pd.concat(
    [
        df[df["subfamily"].isin(group)].assign(entry=entry)
        for entry, group in groups.items()
    ]
)
copies["length"] = copies["end"] - copies["start"]
columns = ["entry", "subfamily", "chrom", "start", "end", "strand", "length"]
columns += ["consensus_length", "fraction", "full_length"]
copies[columns].to_csv(copies_out, sep="\t", index=False, float_format="%.3f")

with open(bed_out, "w") as bed:
    for entry, group in groups.items():
        n_all = df["subfamily"].isin(group).sum()
        copies = full_length[full_length["subfamily"].isin(group)]
        logging.info(
            f"{entry}: {len(copies)} full-length copies out of {n_all} "
            f"(at least {min_fraction:.0%} of the consensus)"
        )
        if copies.empty:
            logging.warning(f"{entry}: no full-length copies, left out")
            continue
        copies = copies.sort_values(["chrom", "start"])
        for row in copies.itertuples():
            bed.write(
                f"{row.chrom}\t{row.start}\t{row.end}\t{row.subfamily}\t"
                f"{row.fraction:.3f}\t{row.strand}\n"
            )
        bed.write(f"#{entry} (n={len(copies)})\n")

if full_length.empty:
    raise ValueError(
        "No full-length copies of the subfamilies in te: profile: subfamilies"
    )
logging.info("Done!")
