"""
Reads the RNA-seq DESeq2 results tables of a contrast (rnaseq: tables: one
table, or several, e.g. genes and TEs in separate files, which are combined)
and matches their rows
to the genes of the GTF annotation and, with the TE module, to the repeat
subfamilies:

- Genes: Ensembl gene IDs (a version suffix, e.g. ENSMUSG00000000001.5, is
  removed) or gene names (rnaseq: id_type: ensembl, symbol, or auto: whichever
  matches more rows).
- TEs (TE module only): rows whose ID, up to the first te_separator
  (default ":"), is a RepeatMasker subfamily, e.g. "L1MdTf_I:L1:LINE"
  (TEtranscripts, TEcount) or "L1MdTf_I". Only the first part is used, as
  family names can be subfamily names too (e.g. MIRb:MIR:SINE).

Rows that match neither are reported in the log and left out. Duplicated
genes are reduced to the row with the lowest padj.

The output has one row per gene or subfamily: id (as in the table), type
(gene or te), gene_id, gene_name, subfamily, baseMean, log2FoldChange, padj
and rna_change: up or down (padj < fdr and |log2FoldChange| >
lfc_threshold), unchanged, or not tested (padj NA, e.g. independent
filtering by DESeq2).
"""

import logging
import re

import numpy as np
import pandas as pd

# Load Snakemake variables
# One table (a string) or several (a list)
table_files = snakemake.input["rnaseq"]
table_files = [table_files] if isinstance(table_files, str) else list(table_files)
table_file = ", ".join(table_files)
gtf_file = snakemake.input["gtf"]
# A list (empty without the TE module)
subfamilies_file = next(iter(snakemake.input.get("subfamilies", [])), None)
columns = snakemake.params["columns"]
id_type = snakemake.params["id_type"]
fdr = snakemake.params["fdr"]
lfc_threshold = snakemake.params["lfc_threshold"]
separator = snakemake.params["te_separator"]
output = snakemake.output["table"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def read_table(path):
    """DESeq2 results as CSV (.csv, .csv.gz) or tab-separated (otherwise)"""
    name = path.lower().removesuffix(".gz")
    sep = "," if name.endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep)


def read_genes(path):
    """Gene IDs and names of the gene lines of a GTF file"""
    gtf = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[2, 8],
        names=["feature", "attributes"],
        dtype=str,
    )
    attributes = gtf.loc[gtf["feature"] == "gene", "attributes"]
    genes = pd.DataFrame(
        {
            "gene_id": attributes.str.extract(r'gene_id "([^"]+)"', expand=False),
            "gene_name": attributes.str.extract(r'gene_name "([^"]+)"', expand=False),
        }
    )
    return genes.dropna(subset=["gene_id"]).drop_duplicates("gene_id")


# RNA-seq tables
# -----------------------------------------------------
def read_rna(path):
    """The ID, baseMean, log2FoldChange and padj columns of one table"""
    table = read_table(path)
    logging.info(f"{len(table)} rows in {path}")
    id_column = columns["gene_id"] or table.columns[0]
    required = {
        "gene_id": id_column,
        "log2fc": columns["log2fc"],
        "padj": columns["padj"],
    }
    missing = [f"{key} ({col})" for key, col in required.items() if col not in table]
    if missing:
        raise ValueError(
            f"Columns not found in {path}: {', '.join(missing)}. "
            f"Columns of the table: {', '.join(map(str, table.columns))}. "
            "Set them in rnaseq: columns"
        )
    return pd.DataFrame(
        {
            "id": table[id_column].astype(str).str.strip(),
            "baseMean": (
                pd.to_numeric(table[columns["base_mean"]], errors="coerce")
                if columns["base_mean"] in table
                else np.nan
            ),
            "log2FoldChange": pd.to_numeric(table[columns["log2fc"]], errors="coerce"),
            "padj": pd.to_numeric(table[columns["padj"]], errors="coerce"),
        }
    )


rna = pd.concat([read_rna(path) for path in table_files], ignore_index=True)

# Genes
# -----------------------------------------------------
genes = read_genes(gtf_file)
logging.info(f"{len(genes)} genes in {gtf_file}")
ids = set(genes["gene_id"])
# Gene names that belong to one gene only
names = genes.dropna(subset=["gene_name"])
names = names[~names["gene_name"].duplicated(keep=False)]
name_to_id = dict(zip(names["gene_name"], names["gene_id"]))
id_to_name = dict(zip(genes["gene_id"], genes["gene_name"]))

# Ensembl IDs without version (also ..._PAR_Y of human chrY genes)
unversioned = rna["id"].str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
by_id = rna["id"].where(rna["id"].isin(ids), unversioned.where(unversioned.isin(ids)))
by_name = rna["id"].map(name_to_id)
n_id, n_name = by_id.notna().sum(), by_name.notna().sum()
logging.info(f"Rows matching a gene ID: {n_id}, a unique gene name: {n_name}")
if id_type == "auto":
    id_type = "ensembl" if n_id >= n_name else "symbol"
    logging.info(f"Gene identifiers: {id_type}")
rna["gene_id"] = by_id if id_type == "ensembl" else by_name
rna["gene_name"] = rna["gene_id"].map(id_to_name)

# TE subfamilies (TE module)
# -----------------------------------------------------
rna["subfamily"] = pd.Series(np.nan, index=rna.index, dtype=object)
if subfamilies_file:
    subfamilies = set(pd.read_csv(subfamilies_file, sep="\t")["subfamily"].astype(str))

    def subfamily_of(row_id):
        """The ID up to the first separator, if it is a subfamily"""
        first = row_id.split(separator)[0] if separator else row_id
        return first if first in subfamilies else None

    # Gene IDs take precedence, unless the ID has the separator
    has_separator = rna["id"].str.contains(re.escape(separator)) if separator else False
    candidates = rna["gene_id"].isna() | has_separator
    rna.loc[candidates, "subfamily"] = rna.loc[candidates, "id"].map(subfamily_of)
    rna.loc[rna["subfamily"].notna(), ["gene_id", "gene_name"]] = np.nan
    logging.info(f"Rows matching a repeat subfamily: {rna['subfamily'].notna().sum()}")

rna["type"] = np.where(
    rna["subfamily"].notna(),
    "te",
    np.where(rna["gene_id"].notna(), "gene", "unmatched"),
)
unmatched = rna.loc[rna["type"] == "unmatched", "id"]
# A table of TE rows only (TE module) is valid
if len(unmatched) == len(rna):
    raise ValueError(
        f"No rows of {table_file} match a gene of {gtf_file} by {id_type}"
        + (" or a repeat subfamily" if subfamilies_file else "")
        + ". Check rnaseq: columns: gene_id and rnaseq: id_type, and that the "
        "RNA-seq data used the same genome annotation"
    )
logging.info(
    f"Rows matching no gene or subfamily (left out): {len(unmatched)}"
    + (f", e.g. {', '.join(unmatched.head(10))}" if len(unmatched) else "")
)
rna = rna[rna["type"] != "unmatched"].copy()

# One row per gene or subfamily: the row with the lowest padj
rna["key"] = rna["gene_id"].where(rna["type"] == "gene", rna["subfamily"])
rna = rna.sort_values("padj", na_position="last")
duplicated = rna["key"].duplicated()
if duplicated.any():
    logging.warning(
        f"{duplicated.sum()} rows of genes or subfamilies that are in the table "
        "more than once were removed (the row with the lowest padj is kept)"
    )
rna = rna[~duplicated]

# Differential expression
# -----------------------------------------------------
significant = (rna["padj"] < fdr) & (rna["log2FoldChange"].abs() > lfc_threshold)
rna["rna_change"] = np.select(
    [
        rna["padj"].isna() | rna["log2FoldChange"].isna(),
        significant & (rna["log2FoldChange"] > 0),
        significant & (rna["log2FoldChange"] < 0),
    ],
    ["not tested", "up", "down"],
    default="unchanged",
)
for kind, group in rna.groupby("type"):
    counts = group["rna_change"].value_counts().to_dict()
    logging.info(f"{kind}: {len(group)} rows, {counts}")

rna = rna.sort_values(["type", "padj"], na_position="last")
rna[
    [
        "id",
        "type",
        "gene_id",
        "gene_name",
        "subfamily",
        "baseMean",
        "log2FoldChange",
        "padj",
        "rna_change",
    ]
].to_csv(output, sep="\t", index=False, na_rep="NA")
logging.info("Done!")
