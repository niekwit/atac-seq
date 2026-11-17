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
rmsk table). RepeatMasker can split one element into several rows; these are
not joined, so an element counts only if one row covers enough of the
consensus. Chromosome names are converted to Ensembl names and copies on
filter_chrs or on chromosomes not in the genome are removed, as in
prepare_te_annotation.py.

In the BED file, each group of regions is followed by a line "#<label>",
which deepTools uses as the group label.
"""

import logging

import pandas as pd

# Load Snakemake variables
rmsk = snakemake.input["rmsk"]
chrom_sizes = snakemake.input["chrom_sizes"]
bed_out = snakemake.output["bed"]
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

# Ensembl chromosome names
df["chrom"] = df["chrom"].str.removeprefix("chr").replace({"M": "MT"})
chroms = pd.read_csv(chrom_sizes, sep="\t", header=None, dtype={0: str})[0]
df = df[df["chrom"].isin(set(chroms) - set(filter_chrs))]

# Consensus positions covered by each copy. On the + strand the alignment
# runs from repStart to repEnd and repLeft is minus the consensus bases after
# it; on the - strand repStart is minus the bases after repEnd and repLeft is
# the start of the alignment
plus = df["strand"] == "+"
first = df["rep_start"].where(plus, df["rep_left"])
after = -df["rep_left"].where(plus, df["rep_start"])
consensus_length = df["rep_end"] + after
df["fraction"] = (df["rep_end"] - first + 1) / consensus_length
full_length = df[df["fraction"] >= min_fraction]

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
