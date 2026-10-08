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
        "grep -v '^{params.mito}\\s' | sort | uniq -c | "
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
        mem_mb=12000,
    log:
        "logs/tss_enrichment/{sample}.log",
    conda:
        "../envs/tss.yaml"
    script:
        "../scripts/tss_enrichment.py"


# Collate FastQC, cutadapt, bowtie2, Picard and samtools stats output
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
    log:
        "logs/encode_qc_summary.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/encode_qc_summary.py"
