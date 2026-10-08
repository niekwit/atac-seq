# Library complexity (NRF, PBC1, PBC2) of the filtered, duplicate marked
# BAM files, excluding mitochondrial reads (ENCODE: pbc_qc_pe)
# -----------------------------------------------------
rule library_complexity:
    input:
        "results/filtered/{sample}.dupmark.bam",
    output:
        "results/qc/{sample}.lib_complexity.qc",
    params:
        mito=config["mito_chr_name"],
        awk=(
            "BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} "
            "{m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; "
            "m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; "
            'printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}'
        ),
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=6000,
    log:
        "logs/library_complexity/{sample}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(printf 'TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\t"
        "NRF\\tPBC1\\tPBC2\\n' > {output}; "
        "samtools sort -@ {threads} -n -T {output}.tmp -o - {input} | "
        "bedtools bamtobed -bedpe -i stdin | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$4,$6,$9,$10}}' | "
        "grep -v '^{params.mito}\\s' | sort -S 2G -T $(dirname {output}) | uniq -c | "
        "awk '{params.awk}' >> {output}) 2> {log}"


# Fraction of mitochondrial reads among the mapped reads (ENCODE: frac_mito)
# -----------------------------------------------------
rule frac_mito:
    input:
        bam="results/bowtie2/{sample}.bam",
        bai="results/bowtie2/{sample}.bam.bai",
    output:
        "results/qc/{sample}.frac_mito.qc",
    params:
        mito=config["mito_chr_name"],
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=1000,
    log:
        "logs/frac_mito/{sample}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(total=$(samtools view -c -F 0x904 {input.bam}); "
        "if samtools idxstats {input.bam} | cut -f 1 | grep -qx '{params.mito}'; then "
        "mito=$(samtools view -c -F 0x904 {input.bam} {params.mito}); else mito=0; fi; "
        "printf 'non_mito_reads\\tmito_reads\\tfrac_mito_reads\\n' > {output}; "
        "awk -v t=$total -v m=$mito "
        "'BEGIN{{printf \"%d\\t%d\\t%f\\n\", t - m, m, (t > 0 ? m / t : 0)}}' >> {output}"
        ") 2> {log}"


# TSS enrichment as calculated by ENCODE (encode_task_tss_enrich.py)
# -----------------------------------------------------
rule tss_enrichment:
    input:
        bam="results/filtered/{sample}.nodup.bam",
        bai="results/filtered/{sample}.nodup.bam.bai",
        tss="resources/tss.bed",
        fastq="reads/{sample}_R1_001.fastq.gz",
    output:
        score="results/qc/{sample}.tss_enrich.qc",
        profile="results/qc/{sample}.tss_enrich_profile.tsv",
        plot="results/qc/{sample}.tss_enrich.png",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=2000,
    log:
        "logs/tss_enrichment/{sample}.log",
    conda:
        "../envs/tss.yaml"
    script:
        "../scripts/tss_enrichment.py"


# Fragment length distribution and nucleosomal QC (ENCODE: fraglen_stat_pe)
# -----------------------------------------------------
rule fraglen_stat:
    input:
        bam="results/filtered/{sample}.nodup.bam",
    output:
        metrics="results/qc/{sample}.insert_size_metrics.txt",
        histogram_pdf="results/qc/{sample}.insert_size_histogram.pdf",
        qc="results/qc/{sample}.nucleosomal.qc",
        plot="results/qc/{sample}.fraglen_dist.png",
    params:
        java_heap="3g",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=4000,
    log:
        "logs/fraglen_stat/{sample}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/fraglen_stat.py"


# GC bias (ENCODE: gc_bias)
# -----------------------------------------------------
rule gc_bias:
    input:
        bam="results/filtered/{sample}.nodup.bam",
        fasta=resources.fasta,
    output:
        metrics="results/qc/{sample}.gc_bias_metrics.txt",
        summary="results/qc/{sample}.gc_bias_summary.txt",
        chart_pdf="results/qc/{sample}.gc_bias.pdf",
        plot="results/qc/{sample}.gc_plot.png",
    params:
        java_heap="8g",
    threads: 1
    resources:
        runtime=runtime(180),
        mem_mb=10000,
    log:
        "logs/gc_bias/{sample}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/gc_bias.py"


# Fraction of reads in annotated regions (ENCODE: annot_enrich)
# -----------------------------------------------------
rule annot_enrich:
    input:
        unpack(lambda w: annotation_regions()),
        ta="results/tagalign/{sample}.tagAlign.gz",
        blacklist="resources/blacklist.bed",
    output:
        "results/qc/{sample}.annot_enrich.qc",
    threads: 1
    resources:
        runtime=runtime(60),
        mem_mb=4000,
    log:
        "logs/annot_enrich/{sample}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/annot_enrich.py"


# Fingerprint and Jensen-Shannon distance of the replicates (ENCODE: jsd)
# -----------------------------------------------------
rule jsd:
    input:
        bams=lambda w: expand(
            "results/filtered/{sample}.nodup.bam", sample=replicates(w.condition)
        ),
        bais=lambda w: expand(
            "results/filtered/{sample}.nodup.bam.bai", sample=replicates(w.condition)
        ),
        blacklist="resources/blacklist.bed",
    output:
        plot="results/qc/{condition}.jsd_plot.png",
        metrics="results/qc/{condition}.jsd.qc",
    params:
        samples=lambda w: replicates(w.condition),
        mapq=30,
    threads: 8
    resources:
        runtime=runtime(180),
        mem_mb=8000,
    log:
        "logs/jsd/{condition}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/jsd.py"


# Collate FastQC, cutadapt, bowtie2, Picard, deepTools and samtools stats output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read=["R1", "R2"],
        ),
        expand("logs/cutadapt/{sample}.log", sample=SAMPLES),
        expand("logs/bowtie2/{sample}.log", sample=SAMPLES),
        expand("results/filtered/{sample}.dup.qc", sample=SAMPLES),
        expand("results/qc/{sample}.insert_size_metrics.txt", sample=SAMPLES),
        expand("results/qc/{sample}.gc_bias_metrics.txt", sample=SAMPLES),
        expand("results/qc/{condition}.jsd.qc", condition=CONDITIONS),
        expand(
            "results/samtools_stats/{stage}.txt",
            stage=[f"bowtie2/{s}" for s in SAMPLES]
            + [f"filtered/{s}.nodup" for s in SAMPLES],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra="--verbose --dirs",
    resources:
        runtime=runtime(30),
        mem_mb=2000,
    log:
        "logs/multiqc.log",
    wrapper:
        "v8.1.1/bio/multiqc"


# Summary of ENCODE QC metrics with ENCODE standards
# -----------------------------------------------------
rule encode_qc_summary:
    input:
        bowtie2=expand("logs/bowtie2/{sample}.log", sample=SAMPLES),
        frac_mito=expand("results/qc/{sample}.frac_mito.qc", sample=SAMPLES),
        lib_complexity=expand("results/qc/{sample}.lib_complexity.qc", sample=SAMPLES),
        nodup_stats=expand(
            "results/samtools_stats/filtered/{sample}.nodup.txt", sample=SAMPLES
        ),
        frip=expand("results/macs2/{sample}.frip.qc", sample=SAMPLES),
        peaks=expand("results/macs2/{sample}.bfilt.narrowPeak", sample=SAMPLES),
        ataqv=expand("results/ataqv/{sample}.json.gz", sample=SAMPLES),
        tss_enrich=expand("results/qc/{sample}.tss_enrich.qc", sample=SAMPLES),
        nucleosomal=expand("results/qc/{sample}.nucleosomal.qc", sample=SAMPLES),
        annot_enrich=expand("results/qc/{sample}.annot_enrich.qc", sample=SAMPLES),
        jsd=expand("results/qc/{condition}.jsd.qc", condition=CONDITIONS),
        gc_bias=expand("results/qc/{sample}.gc_bias_summary.txt", sample=SAMPLES),
        reproducibility=expand(
            "results/{method}/{condition}/reproducibility.qc",
            method=reproducibility_methods(),
            condition=CONDITIONS,
        ),
        frip_pairs=[
            f"results/{method}/{condition}/{pair}.frip.qc"
            for method in reproducibility_methods()
            for condition in CONDITIONS
            for pair in peak_pairs(condition)
        ],
    output:
        samples="results/qc/encode_qc_summary.tsv",
        conditions="results/qc/encode_reproducibility_summary.tsv",
    params:
        samples=SAMPLES,
        conditions=CONDITIONS,
        replicates={c: replicates(c) for c in CONDITIONS},
        methods=reproducibility_methods(),
        genome=config["genome"]["ensembl"],
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        "logs/encode_qc_summary.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/encode_qc_summary.py"
