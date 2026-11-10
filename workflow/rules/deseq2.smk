# Comparison of all samples with DESeq2 (optional, deseq2: enable; not part of
# the ENCODE pipeline, which analyses each condition on its own)


# Consensus peaks: the merged optimal peak sets of all conditions
# -----------------------------------------------------
rule consensus_peaks:
    input:
        peaks=expand(
            "results/{method}/{condition}/optimal_peak.narrowPeak",
            method=deseq2_peak_method(),
            condition=CONDITIONS,
        ),
    output:
        bed="results/deseq2/consensus_peaks.bed",
        saf="results/deseq2/consensus_peaks.saf",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=2000,
    log:
        "logs/deseq2/consensus_peaks.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(cat {input.peaks} | cut -f 1-3 | LC_COLLATE=C sort -S 1G -k1,1 -k2,2n | "
        "bedtools merge -i stdin | "
        'awk \'BEGIN{{OFS="\\t"}} {{print $1, $2, $3, "peak_" NR}}\' > {output.bed} && '
        'awk \'BEGIN{{OFS="\\t"; print "GeneID", "Chr", "Start", "End", "Strand"}} '
        '{{print $4, $1, $2 + 1, $3, "."}}\' {output.bed} > {output.saf}) 2> {log}'


# Count Tn5 insertions (5' end of each read) in the consensus peaks
# -----------------------------------------------------
rule count_peak_insertions:
    input:
        bams=expand("results/filtered/{sample}.nodup.bam", sample=SAMPLES),
        saf="results/deseq2/consensus_peaks.saf",
    output:
        counts="results/deseq2/consensus_peak_counts.txt",
    params:
        # -p without --countReadPairs: each read (Tn5 insertion) is counted
        options="-p -F SAF --read2pos 5 --largestOverlap -s 0",
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=8000,
    log:
        "logs/deseq2/count_peak_insertions.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "featureCounts {params.options} -T {threads} -a {input.saf} "
        "-o {output.counts} {input.bams} > {log} 2>&1"


# PCA and correlation of all samples
# -----------------------------------------------------
rule deseq2_qc:
    input:
        counts="results/deseq2/consensus_peak_counts.txt",
    output:
        vst="results/deseq2/vst_counts.tsv",
        pca="results/deseq2/pca.png",
        pca_data="results/deseq2/pca.tsv",
        correlation="results/deseq2/sample_correlation.png",
        correlation_data="results/deseq2/sample_correlation.tsv",
    params:
        samples=SAMPLES,
        conditions=SAMPLE_SHEET["condition"].tolist(),
        levels=deseq2_levels(),
        ntop=config.get("deseq2", {}).get("pca_ntop", 500),
    threads: 2
    resources:
        runtime=runtime(60),
        mem_mb=6000,
    log:
        "logs/deseq2/deseq2_qc.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_qc.R"


# Annotate the consensus peaks with ChIPseeker
# -----------------------------------------------------
rule annotate_consensus_peaks:
    input:
        peaks="results/deseq2/consensus_peaks.bed",
        edb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        txdb=f"resources/{resources.genome}_{resources.build}_txdb.Rdata",
    output:
        txt="results/deseq2/consensus_peaks_annotated.txt",
    threads: 2
    resources:
        runtime=runtime(30),
        mem_mb=4000,
    log:
        "logs/deseq2/annotate_consensus_peaks.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/annotate_peaks.R"


# Fit the DESeq2 model (~ condition) to all samples
# -----------------------------------------------------
rule deseq2_fit:
    input:
        counts="results/deseq2/consensus_peak_counts.txt",
    output:
        dds="results/deseq2/dds.rds",
    params:
        samples=SAMPLES,
        conditions=SAMPLE_SHEET["condition"].tolist(),
        levels=deseq2_levels(),
    threads: 4
    resources:
        runtime=runtime(120),
        mem_mb=12000,
    log:
        "logs/deseq2/deseq2_fit.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_fit.R"


# Differential accessibility between two conditions
# -----------------------------------------------------
rule deseq2_contrast:
    input:
        dds="results/deseq2/dds.rds",
        annotation="results/deseq2/consensus_peaks_annotated.txt",
        peaks="results/deseq2/consensus_peaks.bed",
    output:
        table="results/deseq2/{contrast}/differential_peaks.tsv",
        ma="results/deseq2/{contrast}/ma_plot.png",
        volcano="results/deseq2/{contrast}/volcano_plot.png",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        fdr=config.get("deseq2", {}).get("fdr", 0.05),
        lfc_threshold=config.get("deseq2", {}).get("lfc_threshold", 0),
    threads: 2
    resources:
        runtime=runtime(60),
        mem_mb=6000,
    log:
        "logs/deseq2/{contrast}.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_contrast.R"
