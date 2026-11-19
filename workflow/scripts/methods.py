"""
Methods section of the run (snakemake-methods). The prose is in
workflow/methods/methods.md.j2; this script derives the values of the
configuration that the template needs. Only the steps enabled in the
configuration are described (ENCODE peaks and QC unless te: only, IDR,
DESeq2, TE module, TE profile, TE correlation and RNA-seq integration).
"""

from snakemake_methods.snakemake import run
from snakemake_methods.text import code, join

config = snakemake.config
params = snakemake.params
reference = params["reference"]

# Settings of the run
# -----------------------------------------------------
te_only = params["te_only"]
te = params["te"]
deseq2 = params["deseq2"]
replicates = params["replicates"]
min_replicate_fraction = config["peaks"].get("min_replicate_fraction", 1)
lfc_threshold = config.get("deseq2", {}).get("lfc_threshold", 0)
macs2 = config["macs2"]
smooth_win = macs2["smooth_win"]
cutadapt = config["cutadapt"]
multimapping = config["bowtie2"]["multimapping"]
seed = config["pseudoreplication_random_seed"]
profile = params["profile"]

# Genome name: "mouse genome GRCm39", from the Ensembl FASTA file name
COMMON_NAMES = {
    "Homo_sapiens": "human",
    "Mus_musculus": "mouse",
    "Drosophila_melanogaster": "Drosophila melanogaster",
}
parts = reference["fasta"].split(".")
if len(parts) < 3:
    genome_name = "genome"
else:
    species = COMMON_NAMES.get(parts[0], parts[0].replace("_", " "))
    genome_name = f"{species} genome {parts[1]}"

# Extra command-line arguments (extra: in the config) of the tools that run
EXTRA_NAMES = {
    "bowtie2_build": "bowtie2-build",
    "fastqc": "FastQC",
    "markduplicates": "Picard MarkDuplicates",
    "ataqv": "ataqv",
    "featurecounts": "featureCounts",
    "multiqc": "MultiQC",
    "bamcoverage": "bamCoverage",
    "bigwigaverage": "bigwigAverage",
    "computematrix": "computeMatrix",
    "plotprofile": "plotProfile",
}
profiles = bool(profile["subfamilies"])
# deepTools profiles: TE profile, TSS profile of the RNA-seq integration (TE
# module), H3K9me3 and DNA methylation profiles (signal around the peaks, with
# BigWigs; ATAC-seq over the regions with the TE module)
deeptools_profiles = (
    (te and profiles)
    or (te and bool(params["rnaseq"]))
    or (te and any(params["epi"].values()))
    or any(s["bigwigs"] for c in params["epi"].values() for s in c.values())
)
EXTRA_RUNS = {
    "markduplicates": not te_only,
    "ataqv": not te_only,
    "featurecounts": (deseq2 and not te_only) or te,
    "bamcoverage": te,
    "bigwigaverage": te,
    "computematrix": deeptools_profiles,
    "plotprofile": deeptools_profiles,
}
extras = [
    f"{EXTRA_NAMES.get(key, key)} {code(value)}"
    for key, value in (config.get("extra") or {}).items()
    if value and EXTRA_RUNS.get(key, True)
]

# R scripts: GTF annotation and ChIPseeker (ENCODE peaks, or TE peaks with
# DESeq2), the AnnotationHub blacklist, DESeq2 and eulerr
annotation = not te_only or deseq2

context = {
    "idr": config["idr"]["enable"] and not te_only,
    "replicated": [c for c, n in replicates.items() if n > 1],
    "samples_text": join(f"{c} (n = {n})" for c, n in replicates.items()),
    "n_samples": sum(replicates.values()),
    "genome_name": genome_name,
    "filter_chrs": config["filter"]["filter_chrs"],
    "all_replicates": min_replicate_fraction >= 1,
    "min_fraction": f"{min_replicate_fraction:.0%}",
    # :g keeps small thresholds (0.004 -> 0.4%) and drops float noise
    "min_overlap": f"{config['peaks'].get('min_overlap', 0.5) * 100:g}%",
    "fdr": config.get("deseq2", {}).get("fdr", 0.05),
    "lfc_text": (
        f"|log2 fold change| > {lfc_threshold}"
        if lfc_threshold
        else "log2 fold change ≠ 0"
    ),
    "count_options": "-p -F SAF --read2pos 5 --largestOverlap -s 0",
    "tss_region": config["peaks"].get("tss_region", 3000),
    "min_count": config.get("deseq2", {}).get("min_count", 10),
    "bin_size": config.get("te", {}).get("bigwig_bin_size", 10),
    "mapq_thresh": config["filter"]["mapq_thresh"],
    "pca_ntop": config.get("deseq2", {}).get("pca_ntop", 500),
    "contrast_text": join(f"{t} vs {c}" for t, c in params["contrasts"].values()),
    "cutadapt_options": (
        f"{cutadapt['extra']} -a {cutadapt['adapter_r1']} "
        f"-A {cutadapt['adapter_r2']}"
    ),
    "multimapping": multimapping,
    "bowtie2_options": (
        (f"-k {multimapping + 1} " if multimapping else "") + config["bowtie2"]["extra"]
    ),
    "macs2_options": (
        f"callpeak -f BED --nomodel --shift {-int(round(smooth_win / 2))} "
        f"--extsize {smooth_win} -p {macs2['pval_thresh']} --keep-dup all "
        "--call-summits -B --SPMR"
    ),
    "keep_chr": "^(" + config["peaks"]["keep_chr_regex"] + ")$",
    "seed_text": (
        "seeded with the tagAlign file size" if seed == 0 else f"random seed {seed}"
    ),
    "extras": extras,
    "rnaseq": params["rnaseq"],
    "rnaseq_bins": join(
        f"{b / 1000:g} kb" if b >= 1000 else f"{b} bp"
        for b in params["rnaseq_config"]["distance_bins"]
    ),
    "rnaseq_max": f"{max(params['rnaseq_config']['distance_bins']) / 1000:g} kb",
    "rnaseq_fdr": params["rnaseq_config"]["fdr"],
    "rnaseq_pair_lfc": params["rnaseq_config"]["pair_lfc_threshold"],
    "rnaseq_tss_flank": params["rnaseq_config"]["tss_profile_flank"],
    "rnaseq_lfc_text": (
        f"|log2 fold change| > {params['rnaseq_config']['lfc_threshold']}"
        if params["rnaseq_config"]["lfc_threshold"]
        else "log2 fold change ≠ 0"
    ),
    **{
        f"{mark}{suffix}": value
        for mark, contrasts in params["epi"].items()
        for suffix, value in {
            "": list(contrasts),
            "_table": any(c["table"] for c in contrasts.values()),
            "_bed": any(not c["table"] for c in contrasts.values()),
            "_bigwigs": any(c["bigwigs"] for c in contrasts.values()),
            "_fdr": params["epi_config"][mark]["fdr"],
            "_min_effect": params["epi_config"][mark]["min_effect"],
            "_flank": params["epi_config"][mark]["profile_flank"],
            "_body": params["epi_config"][mark]["profile_body_length"],
            "_min_copies": params["epi_config"][mark]["min_copies"],
        }.items()
    },
    "epi_bin": params["computematrix_bin_size"],
    "annotation": annotation,
    "r_used": annotation or reference["blacklist"]["source"] == "annotationhub",
}

run(snakemake, context=context)
