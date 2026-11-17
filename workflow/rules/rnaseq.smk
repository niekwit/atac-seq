# Integration with RNA-seq (optional, rnaseq: tables, with deseq2: enable):
# the differential accessibility of a DESeq2 contrast is compared with an
# RNA-seq DESeq2 results table of the same contrast. Genes are compared with
# the ENCODE consensus peaks and with the TE peaks (nearest TSS), and, with
# the TE module, repeat subfamilies (TE rows of the RNA-seq table) with the
# subfamily counts and TE peaks. Not part of the ENCODE pipeline.


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
        table="results/integration/rnaseq/{contrast}.tsv",
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
        "logs/integration/rnaseq_prepare_{contrast}.log",
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
        rnaseq="results/integration/rnaseq/{contrast}.tsv",
    output:
        genes="results/integration/{contrast}/genes.tsv",
        concordant="results/integration/{contrast}/concordant_genes.tsv",
        cdf="results/integration/{contrast}/cdf.pdf",
        cdf_tests="results/integration/{contrast}/cdf_tests.tsv",
        enrichment="results/integration/{contrast}/enrichment.tsv",
        enrichment_plot="results/integration/{contrast}/enrichment.pdf",
        correlation="results/integration/{contrast}/correlation.tsv",
        scatter="results/integration/{contrast}/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        features="peaks",
        distance_bins=rnaseq_config()["distance_bins"],
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=6000,
        cpus_per_task=cpus_for_memory(6000),
    log:
        "logs/integration/rnaseq_genes_{contrast}.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/rnaseq_genes.R"


# Differential expression of genes against the differential accessibility of
# the TE peaks (TE module)
# -----------------------------------------------------
use rule rnaseq_genes as te_rnaseq_genes with:
    input:
        atac="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        rnaseq="results/integration/rnaseq/{contrast}.tsv",
    output:
        genes="results/te/integration/{contrast}/peaks/genes.tsv",
        concordant="results/te/integration/{contrast}/peaks/concordant_genes.tsv",
        cdf="results/te/integration/{contrast}/peaks/cdf.pdf",
        cdf_tests="results/te/integration/{contrast}/peaks/cdf_tests.tsv",
        enrichment="results/te/integration/{contrast}/peaks/enrichment.tsv",
        enrichment_plot="results/te/integration/{contrast}/peaks/enrichment.pdf",
        correlation="results/te/integration/{contrast}/peaks/correlation.tsv",
        scatter="results/te/integration/{contrast}/peaks/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        features="TE peaks",
        distance_bins=rnaseq_config()["distance_bins"],
    log:
        "logs/te/integration/rnaseq_genes_{contrast}.log",


# Differential expression of repeat subfamilies (TE rows of the RNA-seq
# table) against their differential accessibility (TE module)
# -----------------------------------------------------
rule te_rnaseq_subfamilies:
    input:
        atac="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
        peaks="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        rnaseq="results/integration/rnaseq/{contrast}.tsv",
    output:
        subfamilies="results/te/integration/{contrast}/subfamilies/subfamilies.tsv",
        correlation="results/te/integration/{contrast}/subfamilies/correlation.tsv",
        scatter="results/te/integration/{contrast}/subfamilies/scatter.pdf",
        cdf="results/te/integration/{contrast}/subfamilies/cdf.pdf",
        cdf_tests="results/te/integration/{contrast}/subfamilies/cdf_tests.tsv",
        enrichment="results/te/integration/{contrast}/subfamilies/enrichment.tsv",
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
        "../envs/deseq2.yaml"
    script:
        "../scripts/rnaseq_subfamilies.R"
