"""
QC summary of the TE alignments per sample: read pairs kept, duplicates,
multimappers and the fraction on filter_chrs, by default the mitochondrial
chromosome (te_library_stats.tsv), TSS enrichment, the fraction of reads in
nucleosome-free regions (NFR) with the NFR / mono-nucleosome ratio, and the
fraction of reads in the peaks of the sample (FRiP). Also written as MultiQC
custom content (general statistics table).
"""

import logging

import pandas as pd

# Load Snakemake variables
stats_file = snakemake.input["stats"]
tss_files = snakemake.input["tss"]
nucleosomal_files = snakemake.input["nucleosomal"]
frip_files = snakemake.input["frip"]
samples = snakemake.params["samples"]
filter_chrs = list(snakemake.params["filter_chrs"])
mito = snakemake.params["mito"]
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
# te_filter removes the duplicates first, so the reads on filter_chrs are
# counted among the deduplicated read pairs
summary["frac_filter_chrs_pairs"] = stats["filter_chrs_pairs"] / (
    stats["proper_pairs"] - stats["duplicate_pairs"]
)
summary["frac_multimapping_pairs"] = stats["frac_multimapping_pairs"]
summary["frac_insertions_in_repeats"] = stats["frac_insertions_in_repeats"]

for sample, tss, nucleosomal, frip in zip(
    samples, tss_files, nucleosomal_files, frip_files
):
    with open(tss) as fh:
        summary.loc[sample, "tss_enrichment"] = float(fh.read().strip())
    metrics = pd.read_csv(nucleosomal, sep="\t", index_col="metric")["value"]
    summary.loc[sample, "frac_reads_in_nfr"] = float(
        metrics["Fraction of reads in NFR"]
    )
    summary.loc[sample, "nfr_mono_nuc_ratio"] = float(metrics["NFR / mono-nuc reads"])
    with open(frip) as fh:
        summary.loc[sample, "frip"] = float(fh.read().strip())

summary = summary.round(4)
summary.to_csv(summary_out, sep="\t")

# MultiQC custom content: general statistics table. The reads on filter_chrs
# are the mitochondrial reads with the default filter_chrs (only the
# mitochondrial chromosome)
if filter_chrs == [mito]:
    filter_title = "TE mito"
    filter_description = (
        f"Fraction of deduplicated read pairs on the mitochondrial chromosome ({mito})"
    )
else:
    filter_title = "TE filter chrs"
    filter_description = (
        "Fraction of deduplicated read pairs on filter_chrs "
        f"({', '.join(filter_chrs)})"
    )
header = "\n".join(
    [
        "# plot_type: 'generalstats'",
        "# pconfig:",
        "#     frac_duplicate_pairs:",
        "#         title: 'TE dup'",
        "#         description: 'Fraction of duplicate read pairs (TE alignment)'",
        "#         format: '{:,.3f}'",
        "#     frac_filter_chrs_pairs:",
        f"#         title: '{filter_title}'",
        f"#         description: '{filter_description}'",
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
        "#     frip:",
        "#         title: 'TE FRiP'",
        "#         description: 'Fraction of reads in the peaks of the sample (TE peaks)'",
        "#         format: '{:,.3f}'",
    ]
)
columns = [
    "frac_duplicate_pairs",
    "frac_filter_chrs_pairs",
    "frac_multimapping_pairs",
    "frac_insertions_in_repeats",
    "tss_enrichment",
    "frac_reads_in_nfr",
    "frip",
]
with open(mqc_out, "w") as fh:
    fh.write(header + "\n")
    summary[columns].to_csv(fh, sep="\t")
logging.info("Done!")
