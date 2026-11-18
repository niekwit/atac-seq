# Integration with RNA-seq (optional, rnaseq: tables, with deseq2: enable):
# the differential accessibility of a DESeq2 contrast is compared with an
# RNA-seq DESeq2 results table of the same contrast. Genes are compared with
# the ENCODE consensus peaks and with the TE peaks (nearest TSS), and, with
# the TE module, repeat subfamilies (TE rows of the RNA-seq table) with the
# subfamily counts and TE peaks. Not part of the ENCODE pipeline.


# Prepared RNA-seq table, next to the figures of the contrast: with te: only,
# in the TE module results, so that all output is in results/te/integration/
RNASEQ_TABLE = (
    "results/te/integration/rna-seq/{contrast}/rnaseq_table.tsv"
    if te_only()
    else "results/integration/rna-seq/{contrast}/rnaseq_table.tsv"
)
RNASEQ_LOG = (
    "logs/te/integration/rnaseq_prepare_{contrast}.log"
    if te_only()
    else "logs/integration/rnaseq_prepare_{contrast}.log"
)


# Genes and TE subfamilies of the RNA-seq table, with their differential
# expression
# -----------------------------------------------------
rule rnaseq_prepare:
    input:
        rnaseq=lambda w: rnaseq_contrasts()[w.contrast],
        gtf=resources.gtf,
        # Subfamilies only with the TE module
        subfamilies=lambda w: (["resources/te/subfamilies.tsv"] if te_enabled() else []),
    output:
        table=RNASEQ_TABLE,
    params:
        columns=rnaseq_config()["columns"],
        id_type=rnaseq_config()["id_type"],
        fdr=rnaseq_config()["fdr"],
        lfc_threshold=rnaseq_config()["lfc_threshold"],
        te_separator=rnaseq_config()["te_separator"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        RNASEQ_LOG,
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/rnaseq_prepare.py"


# Differential expression of genes against the differential accessibility of
# the ENCODE consensus peaks
# -----------------------------------------------------
rule rnaseq_genes:
    input:
        atac="results/deseq2/{contrast}/differential_peaks.tsv",
        rnaseq=RNASEQ_TABLE,
    output:
        genes="results/integration/rna-seq/{contrast}/genes.tsv",
        concordant="results/integration/rna-seq/{contrast}/concordant_genes.tsv",
        enrichment="results/integration/rna-seq/{contrast}/enrichment.tsv",
        enrichment_plot="results/integration/rna-seq/{contrast}/enrichment.pdf",
        correlation="results/integration/rna-seq/{contrast}/correlation.tsv",
        pairs="results/integration/rna-seq/{contrast}/pairs.tsv",
        scatter="results/integration/rna-seq/{contrast}/scatter.pdf",
        distance="results/integration/rna-seq/{contrast}/distance.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        features="peaks",
        peak_set="ENCODE consensus peaks (uniquely mapped reads)",
        distance_bins=rnaseq_config()["distance_bins"],
        pair_lfc_threshold=rnaseq_config()["pair_lfc_threshold"],
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=6000,
        cpus_per_task=cpus_for_memory(6000),
    log:
        "logs/integration/rnaseq_genes_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/rnaseq_genes.R"


# Differential expression of genes (non-TE rows of the RNA-seq table) against
# the differential accessibility of the peaks of the TE module (genome-wide
# peaks of the alignment that keeps multimapping reads, not only on TEs)
# -----------------------------------------------------
use rule rnaseq_genes as te_rnaseq_genes with:
    input:
        atac="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        rnaseq=RNASEQ_TABLE,
    output:
        genes="results/te/integration/rna-seq/{contrast}/genes/genes.tsv",
        concordant="results/te/integration/rna-seq/{contrast}/genes/concordant_genes.tsv",
        enrichment="results/te/integration/rna-seq/{contrast}/genes/enrichment.tsv",
        enrichment_plot="results/te/integration/rna-seq/{contrast}/genes/enrichment.pdf",
        correlation="results/te/integration/rna-seq/{contrast}/genes/correlation.tsv",
        pairs="results/te/integration/rna-seq/{contrast}/genes/pairs.tsv",
        scatter="results/te/integration/rna-seq/{contrast}/genes/scatter.pdf",
        distance="results/te/integration/rna-seq/{contrast}/genes/distance.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        features="peaks",
        peak_set="TE module peaks (alignment keeping multimapping reads)",
        distance_bins=rnaseq_config()["distance_bins"],
        pair_lfc_threshold=rnaseq_config()["pair_lfc_threshold"],
    log:
        "logs/te/integration/rnaseq_genes_{contrast}.log",


# Differential expression of repeat subfamilies (TE rows of the RNA-seq
# table) against their differential accessibility (TE module)
# -----------------------------------------------------
rule te_rnaseq_subfamilies:
    input:
        atac="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
        peaks="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        rnaseq=RNASEQ_TABLE,
    output:
        subfamilies="results/te/integration/rna-seq/{contrast}/subfamilies/subfamilies.tsv",
        correlation="results/te/integration/rna-seq/{contrast}/subfamilies/correlation.tsv",
        scatter="results/te/integration/rna-seq/{contrast}/subfamilies/scatter.pdf",
        cdf="results/te/integration/rna-seq/{contrast}/subfamilies/cdf.pdf",
        cdf_tests="results/te/integration/rna-seq/{contrast}/subfamilies/cdf_tests.tsv",
        enrichment="results/te/integration/rna-seq/{contrast}/subfamilies/enrichment.tsv",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=3000,
    log:
        "logs/te/integration/rnaseq_subfamilies_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/rnaseq_subfamilies.R"
