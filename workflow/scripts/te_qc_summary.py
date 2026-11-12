"""
QC summary of the TE alignments per sample: read pairs kept, duplicates and
multimappers (te_library_stats.tsv), TSS enrichment and the fraction of reads
in nucleosome-free regions (NFR) with the NFR / mono-nucleosome ratio. Also
written as MultiQC custom content (general statistics table).
"""

import logging

import pandas as pd

# Load Snakemake variables
stats_file = snakemake.input["stats"]
tss_files = snakemake.input["tss"]
nucleosomal_files = snakemake.input["nucleosomal"]
samples = snakemake.params["samples"]
summary_out = snakemake.output["summary"]
mqc_out = snakemake.output["mqc"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

stats = pd.read_csv(stats_file, sep="\t", index_col="sample").loc[samples]
summary = pd.DataFrame(index=pd.Index(samples, name="sample"))
summary["proper_pairs"] = stats["proper_pairs"]
summary["kept_pairs"] = stats["kept_pairs"]
summary["frac_duplicate_pairs"] = stats["frac_duplicate_pairs"]
summary["frac_filter_chrs_pairs"] = stats["filter_chrs_pairs"] / stats["proper_pairs"]
summary["frac_multimapping_pairs"] = stats["frac_multimapping_pairs"]
summary["frac_insertions_in_repeats"] = stats["frac_insertions_in_repeats"]

for sample, tss, nucleosomal in zip(samples, tss_files, nucleosomal_files):
    with open(tss) as fh:
        summary.loc[sample, "tss_enrichment"] = float(fh.read().strip())
    metrics = pd.read_csv(nucleosomal, sep="\t", index_col="metric")["value"]
    summary.loc[sample, "frac_reads_in_nfr"] = float(
        metrics["Fraction of reads in NFR"]
    )
    summary.loc[sample, "nfr_mono_nuc_ratio"] = float(metrics["NFR / mono-nuc reads"])

summary = summary.round(4)
summary.to_csv(summary_out, sep="\t")

# MultiQC custom content: general statistics table
header = "\n".join(
    [
        "# plot_type: 'generalstats'",
        "# pconfig:",
        "#     frac_duplicate_pairs:",
        "#         title: 'TE dup'",
        "#         description: 'Fraction of duplicate read pairs (TE alignment)'",
        "#         format: '{:,.3f}'",
        "#     frac_multimapping_pairs:",
        "#         title: 'TE multi'",
        "#         description: 'Fraction of multimapping read pairs'",
        "#         format: '{:,.3f}'",
        "#     frac_insertions_in_repeats:",
        "#         title: 'In repeats'",
        "#         description: 'Fraction of Tn5 insertions in repeats'",
        "#         format: '{:,.3f}'",
        "#     tss_enrichment:",
        "#         title: 'TSS enrich.'",
        "#         description: 'TSS enrichment (TE alignment)'",
        "#         format: '{:,.2f}'",
        "#     frac_reads_in_nfr:",
        "#         title: 'NFR'",
        "#         description: 'Fraction of reads in nucleosome-free regions'",
        "#         format: '{:,.3f}'",
    ]
)
columns = [
    "frac_duplicate_pairs",
    "frac_multimapping_pairs",
    "frac_insertions_in_repeats",
    "tss_enrichment",
    "frac_reads_in_nfr",
]
with open(mqc_out, "w") as fh:
    fh.write(header + "\n")
    summary[columns].to_csv(fh, sep="\t")
logging.info("Done!")
