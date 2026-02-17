# Download fasta file from Ensembl
# -----------------------------------------------------
rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log",
    threads: 1
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/get_resource.sh"


# Create chromosome sizes file from fasta file for use in downstream tools
# -----------------------------------------------------
rule create_chrom_sizes:
    input:
        resources.fasta,
    output:
        "resources/chrom_sizes.txt",
    log:
        "logs/resources/create_chrom_sizes.log",
    threads: 1
    conda:
        "../envs/atac.yaml"
    shell:
        "samtools faidx {input}; "
        "cut -f1,2 {input}.fai > {output}"


# Download gtf file from Ensembl
# -----------------------------------------------------
use rule get_fasta as get_gtf with:
    output:
        resources.gtf,
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",


# Generate BED file of TSS regions from GTF file for ataqv
# -----------------------------------------------------
rule generate_tss_file:
    input:
        gtf=resources.gtf,
    output:
        bed="resources/tss.bed",
    log:
        "logs/resources/generate_tss_file.log",
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/generate_tss_file.py"


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
        "logs/resources/create_blacklist_bed.log",
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
        report="results/multiqc/fastq.html",
    params:
        extra="--verbose",
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
        fasta_fwd="results/trimmed/{sample}_R1.fq.gz",
        report_fwd="results/trimmed/reports/{sample}_R1_trimming_report.txt",
        fasta_rev="results/trimmed/{sample}_R2.fq.gz",
        report_rev="results/trimmed/reports/{sample}_R2_trimming_report.txt",
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
        reads=["results/trimmed/{sample}_R1.fq.gz", "results/trimmed/{sample}_R2.fq.gz"],
        idx=multiext(
            "resources/bwa_mem2_index/index",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
            ".0123",
        ),
    output:
        "results/mapped/{sample}.bam",
    log:
        "logs/bwa_mem2/{sample}.log",
    params:
        sort=config["bwa_mem2"]["sort"],
        sort_order=config["bwa_mem2"]["sort_order"],
        sort_extra=config["bwa_mem2"]["sort_extra"],
        extra=config["bwa_mem2"]["extra"] + " -R '@RG\\tID:{sample}\\tSM:{sample}\\tLB:lib1\\tPL:ILLUMINA'",
    threads: 8
    wrapper:
        "v8.1.1/bio/bwa-mem2/mem"


# Mark duplicates with Picard
# -----------------------------------------------------
rule markduplicates:
    input:
        bams="results/mapped/{sample}.bam",
    output:
        bam="results/dedup/{sample}.bam",
        metrics="results/dedup/{sample}.metrics.txt",
    log:
        "logs/mark_duplicates/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES false",
    resources:
        mem_mb=4096,
    wrapper:
        "v9.0.0/bio/picard/markduplicates"


# Index dedup marked BAM files with samtools
# -----------------------------------------------------
rule index_dedup_marked:
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


# Remove reads mapping to mitochondrial genome and
# those with low mapping quality (MAPQ < 30)
# -----------------------------------------------------
rule filter_bam:
    input:
        "results/dedup/{sample}.bam",
    output:
        "results/filtered/{sample}.bam",
    log:
        "logs/filter_bam/{sample}.log",
    params:
        q="30",
    threads: 4
    conda:
        "../envs/atac.yaml"
    shell:
        "samtools view -h -q {params.q} {input} | "
        """awk '$3 != "MT"' | """
        "samtools view -Sb - > {output}"


# Index filtered BAM files with samtools
# -----------------------------------------------------
# Index filteredBAM files with samtools
# -----------------------------------------------------
rule index_filtered:
    input:
        "results/filtered/{sample}.bam",
    output:
        "results/filtered/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 3  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"


# Get BAM file statistics with samtools stats
# -----------------------------------------------------
rule samtools_stats:
    input:
        bam="results/dedup/{sample}.bam",
    output:
        "results/samtools_stats/{sample}.txt",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "logs/samtools_stats/{sample}.log",
    wrapper:
        "v8.1.1/bio/samtools/stats"


# Collate BAM file statistics with MultiQC
# -----------------------------------------------------
rule multiqc_samtools_stats:
    input:
        expand("results/samtools_stats/{sample}.txt", sample=SAMPLES),
    output:
        report="results/multiqc/samtools_stats.html",
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    wrapper:
        "v8.1.1/bio/multiqc"
