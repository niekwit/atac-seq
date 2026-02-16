# Download fasta file from Ensembl
# -----------------------------------------------------
rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    threads: 1
    conda:
        "../envs/mapping.yaml"
    script:
        "../scripts/get_resource.sh"


# Download gtf file from Ensembl
# -----------------------------------------------------
use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


rule generate_tss_file:
    input:
        gtf=resources.gtf,
    output:
        "resources/tss.bed",
    log:
        "logs/resources/generate_tss_file.log",
    conda:
        "../envs/mapping.yaml"
    shell:
        "python workflow/scripts/generate_tss_file.py {input.gtf} {output}"


# Index genome with bwa_mem2
# -----------------------------------------------------
rule bwa_mem2_index:
    input:
        resources.fasta,
    output:
        "resources/bwa_mem2_index/index.0123",
        "resources/bwa_mem2_index/index.amb",
        "resources/bwa_mem2_index/index.ann",
        "resources/bwa_mem2_index/index.bwt.2bit.64",
        "resources/bwa_mem2_index/index.pac",
    threads: 24
    log:
        "logs/bwa-mem2/index.log",
    wrapper:
        "v7.2.0/bio/bwa-mem2/index"


# Create BED file of blacklisted regions
# -----------------------------------------------------
rule create_blacklist_bed:
    output:
        bed="resources/blacklist.bed",
    params:
        genome=resources.genome,
    log:
        "logs/resources/create_blacklist_bed.log"
    threads: 1
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/create_blacklist_bed.R"


# Make QC report
# -----------------------------------------------------
rule fastqc:
    input:
        fastq="reads/{sample}_{read}_001.fastq.gz",
    output:
        html="results/fastqc/{sample}.bwa.{read}_fastqc.html",
        zip="results/fastqc/{sample}.bwa.{read}_fastqc.zip",
    params:
        extra="--quiet --memory 1024",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "results/fastqc/{sample}.bwa.{read}.log",
    threads: 4
    wrapper:
        "v6.0.0/bio/fastqc"


# Run multiQC on tool output
# -----------------------------------------------------
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}.bwa.{read}_fastqc.{ext}",
            sample=SAMPLES,
            read=["R1", "R2"],
            ext=["html", "zip"],
        ),
    output:
        report="results/multiqc/multiqc_report.html",
    params:
        extra="--verbose --dirs",
    message:
        """--- Generating MultiQC report for seq data."""
    log:
        "results/multiqc/multiqc.log",
    wrapper:
        "v6.0.0/bio/multiqc"


# Adapter trimming
# -----------------------------------------------------
rule trim_galore_pe:
    input:
        ["reads/{sample}_R1_001.fastq.gz", "reads/{sample}_R2_001.fastq.gz"],
    output:
        fasta_fwd="trimmed/{sample}_R1.fq.gz",
        report_fwd="trimmed/reports/{sample}_R1_trimming_report.txt",
        fasta_rev="trimmed/{sample}_R2.fq.gz",
        report_rev="trimmed/reports/{sample}_R2_trimming_report.txt",
    threads: 4
    params:
        extra=config["trim_galore"]["extra"],
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v3.14.1/bio/trim_galore/pe"


# Mapping with bwa_mem2
# -----------------------------------------------------
rule bwa_mem2_mem:
    input:
        reads=["trimmed/{sample}_R1.fq.gz", "trimmed/{sample}_R2.fq.gz"],
        # Index needs to be a list of all index files created by bwa
        idx=multiext("resources/bwa_mem2_index/index", ".amb", ".ann", ".bwt.2bit.64", ".pac", ".0123"),
    output:
        "results/mapped/{sample}.bam",
    log:
        "logs/bwa_mem2/{sample}.log",
    params:
        extra= config["bwa_mem2"]["extra"] + r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="picard",  # Can be 'none', 'samtools', or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard sorts.
    threads: 8
    wrapper:
        "v8.1.1/bio/bwa-mem2/mem"


# Remove reads mapping to mitochondrial genome and
# those with low mapping quality (MAPQ < 30)
# -----------------------------------------------------
rule filter_bam:
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/filtered/{sample}.bam",
    log:
        "logs/filter_bam/{sample}.log",
    params:
        q="30",
    threads: 4
    shell:
        "samtools view -h -q {params.q} {input} | "
        """awk '$3 != "MT"' | """
        "samtools view -Sb - > {output}"


# Mark duplicates with Picard
# -----------------------------------------------------
rule markduplicates_bam:
    input:
        bams="results/filtered/{sample}.bam",
    output:
        bam="results/dedup/{sample}.bam",
        metrics="results/dedup/{sample}.metrics.txt",
    log:
        "logs/mark_duplicates/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES false",
    resources:
        mem_mb=2048,
    wrapper:
        "v9.0.0/bio/picard/markduplicates"


# Index BAM files with samtools
# -----------------------------------------------------
rule samtools_index:
    input:
        "results/dedup/{sample}.bam",
    output:
        "results/dedup/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 3  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"
