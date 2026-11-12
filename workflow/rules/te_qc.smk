# QC of the TE alignments (TE module), mainly for te: only, where the ENCODE
# alignment and its QC are not run. Same tools as the ENCODE QC, but on the
# alignments that keep multimappers, so the values are not directly
# comparable with the ENCODE standards.


use rule tss_enrichment as te_tss_enrichment with:
    input:
        bam="results/te/{sample}.sorted.bam",
        bai="results/te/{sample}.sorted.bam.bai",
        tss="resources/tss.bed",
        fastq="reads/{sample}_R1_001.fastq.gz",
    output:
        score="results/te/qc/{sample}.tss_enrich.qc",
        profile="results/te/qc/{sample}.tss_enrich_profile.tsv",
        plot="results/te/qc/{sample}.tss_enrich.png",
    log:
        "logs/te/qc/tss_enrichment_{sample}.log",


use rule fraglen_stat as te_fraglen_stat with:
    input:
        bam="results/te/{sample}.sorted.bam",
    output:
        metrics="results/te/qc/{sample}.insert_size_metrics.txt",
        histogram_pdf="results/te/qc/{sample}.insert_size_histogram.pdf",
        qc="results/te/qc/{sample}.nucleosomal.qc",
        plot="results/te/qc/{sample}.fraglen_dist.png",
    log:
        "logs/te/qc/fraglen_stat_{sample}.log",


use rule gc_bias as te_gc_bias with:
    input:
        bam="results/te/{sample}.sorted.bam",
        fasta=resources.fasta,
    output:
        metrics="results/te/qc/{sample}.gc_bias_metrics.txt",
        summary="results/te/qc/{sample}.gc_bias_summary.txt",
        chart_pdf="results/te/qc/{sample}.gc_bias.pdf",
        plot="results/te/qc/{sample}.gc_plot.png",
    log:
        "logs/te/qc/gc_bias_{sample}.log",


rule te_samtools_stats:
    input:
        "results/te/{sample}.sorted.bam",
    output:
        "results/te/samtools_stats/{sample}.txt",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=1000,
    log:
        "logs/te/qc/samtools_stats_{sample}.log",
    conda:
        "../envs/te.yaml"
    shell:
        "samtools stats {input} > {output} 2> {log}"


# Summary per sample: read pairs, duplicates, multimappers, TSS enrichment and
# nucleosome-free fraction, also as MultiQC custom content
# -----------------------------------------------------
rule te_qc_summary:
    input:
        stats="results/te/te_library_stats.tsv",
        tss=expand("results/te/qc/{sample}.tss_enrich.qc", sample=SAMPLES),
        nucleosomal=expand("results/te/qc/{sample}.nucleosomal.qc", sample=SAMPLES),
    output:
        summary="results/te/qc/te_qc_summary.tsv",
        mqc="results/te/qc/te_qc_summary_mqc.tsv",
    params:
        samples=SAMPLES,
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        "logs/te/qc/te_qc_summary.log",
    conda:
        "../envs/te.yaml"
    script:
        "../scripts/te_qc_summary.py"


rule te_multiqc:
    input:
        expand(
            "results/fastqc/{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read=["R1", "R2"],
        ),
        expand("logs/cutadapt/{sample}.log", sample=SAMPLES),
        expand("logs/te/te_align_{sample}.log", sample=SAMPLES),
        expand("results/te/samtools_stats/{sample}.txt", sample=SAMPLES),
        expand("results/te/qc/{sample}.insert_size_metrics.txt", sample=SAMPLES),
        expand("results/te/qc/{sample}.gc_bias_metrics.txt", sample=SAMPLES),
        "results/te/qc/te_qc_summary_mqc.tsv",
    output:
        report="results/te/multiqc/multiqc_report.html",
    params:
        # te_align_{sample}.log: sample name without te_align_
        extra=(
            "--verbose --no-data-dir "
            "--cl-config \"extra_fn_clean_trim: ['te_align_']\""
        ),
        # MultiQC searches the directories of the input files; logs/te/ also
        # has other logs, which MultiQC ignores
        dirs=lambda w, input: sorted({os.path.dirname(f) for f in input}),
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=2000,
    log:
        "logs/te/multiqc.log",
    conda:
        "../envs/reports.yaml"
    shell:
        "multiqc {params.extra} --force --outdir results/te/multiqc "
        "--filename multiqc_report {params.dirs} > {log} 2>&1"
