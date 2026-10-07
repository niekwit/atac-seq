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
    shell:
        # Non-nucleotide characters in sequence lines (which break bowtie2-build)
        # are replaced by N, keeping the coordinates (| as sed delimiter, as
        # snakemake --lint mistakes a quoted string starting with / for a path)
        "(wget -q {params.url} -O - | pigz -dc | "
        "sed '\\|^>|!s|[^ACGTNacgtn]|N|g' > {output}) 2> {log}"


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


# BED file of all chromosomes except those whose reads are removed
# after deduplication (ENCODE filter_chrs)
# -----------------------------------------------------
rule create_keep_chroms_bed:
    input:
        "resources/chrom_sizes.txt",
    output:
        "resources/keep_chroms.bed",
    params:
        filter_chrs=" ".join(config["filter"]["filter_chrs"]),
    log:
        "logs/resources/create_keep_chroms_bed.log",
    threads: 1
    conda:
        "../envs/encode.yaml"
    shell:
        "awk -v chrs='{params.filter_chrs}' "
        '\'BEGIN{{OFS="\\t"; n=split(chrs, a, " "); for (i=1; i<=n; i++) skip[a[i]]=1}} '
        "!($1 in skip){{print $1, 0, $2}}' {input} > {output} 2> {log}"


# Create annotation files for use in downstream tools
# -----------------------------------------------------
rule create_annotation_file:
    input:
        gtf=resources.gtf,
    output:
        rdata=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        txdb=f"resources/{resources.genome}_{resources.build}_txdb.Rdata",
    log:
        "logs/resources/create_annotation_file.log",
    threads: 2
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/create_annotation_file.R"


# Download gtf file from Ensembl
# -----------------------------------------------------
rule get_gtf:
    output:
        resources.gtf,
    retries: 3
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log",
    threads: 1
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/get_resource.sh"


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


# Index genome with bowtie2
# -----------------------------------------------------
rule bowtie2_index:
    input:
        resources.fasta,
    output:
        multiext(
            "resources/bowtie2_index/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        prefix=lambda w, output: output[0].removesuffix(".1.bt2"),
    threads: 24
    log:
        "logs/bowtie2/index.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "bowtie2-build --threads {threads} {input} {params.prefix} > {log} 2>&1"


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
        html="results/fastqc/{sample}_{read}_fastqc.html",
        zip="results/fastqc/{sample}_{read}_fastqc.zip",
    params:
        extra="--quiet --memory 1024",
    message:
        """--- Checking fastq files with FastQC."""
    log:
        "logs/fastqc/{sample}_{read}.log",
    threads: 4
    wrapper:
        "v6.0.0/bio/fastqc"


# Adapter trimming (ENCODE: cutadapt -e 0.1 -m 5)
# -----------------------------------------------------
rule cutadapt:
    input:
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz",
    output:
        r1=temp("results/trimmed/{sample}_R1.fastq.gz"),
        r2=temp("results/trimmed/{sample}_R2.fastq.gz"),
    params:
        adapter_r1=config["cutadapt"]["adapter_r1"],
        adapter_r2=config["cutadapt"]["adapter_r2"],
        extra=config["cutadapt"]["extra"],
    threads: 4
    log:
        "logs/cutadapt/{sample}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "cutadapt -j {threads} {params.extra} "
        "-a {params.adapter_r1} -A {params.adapter_r2} "
        "-o {output.r1} -p {output.r2} "
        "{input.r1} {input.r2} > {log}"


# Mapping with bowtie2 (ENCODE: -k multimapping+1 -X2000 --mm)
# -----------------------------------------------------
rule bowtie2_align:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        idx=rules.bowtie2_index.output,
    output:
        "results/bowtie2/{sample}.bam",
    params:
        prefix=lambda w, input: input.idx[0].removesuffix(".1.bt2"),
        multimapping=(
            "-k {}".format(config["bowtie2"]["multimapping"] + 1)
            if config["bowtie2"]["multimapping"]
            else ""
        ),
        extra=config["bowtie2"]["extra"],
    threads: 8
    log:
        bowtie2="logs/bowtie2/{sample}.log",
        sort="logs/bowtie2/{sample}.sort.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "bowtie2 {params.multimapping} {params.extra} --threads {threads} "
        "--rg-id {wildcards.sample} --rg SM:{wildcards.sample} --rg LB:{wildcards.sample} --rg PL:ILLUMINA "
        "-x {params.prefix} -1 {input.r1} -2 {input.r2} 2> {log.bowtie2} | "
        "samtools sort -@ {threads} -T {output}.tmp -o {output} - 2> {log.sort}"


# Remove unmapped, low quality, orphan and multimapping reads
# -----------------------------------------------------
rule filter_bam:
    input:
        "results/bowtie2/{sample}.bam",
    output:
        temp("results/filtered/{sample}.filt.bam"),
    params:
        multimapping=config["bowtie2"]["multimapping"],
        mapq_thresh=config["filter"]["mapq_thresh"],
        assign_multimappers=f"{workflow.basedir}/scripts/assign_multimappers.py",
    threads: 4
    log:
        "logs/filter_bam/{sample}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/filter_bam.sh"


# Mark duplicates with Picard
# -----------------------------------------------------
rule markduplicates:
    input:
        bams="results/filtered/{sample}.filt.bam",
    output:
        bam="results/filtered/{sample}.dupmark.bam",
        metrics="results/filtered/{sample}.dup.qc",
    log:
        "logs/mark_duplicates/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT",
    resources:
        mem_mb=4096,
    wrapper:
        "v9.0.0/bio/picard/markduplicates"


# Remove duplicates and reads on filter_chrs (e.g. mitochondrial reads)
# -----------------------------------------------------
rule remove_duplicates:
    input:
        bam="results/filtered/{sample}.dupmark.bam",
        bed="resources/keep_chroms.bed",
    output:
        "results/filtered/{sample}.nodup.bam",
    threads: 4
    log:
        "logs/remove_duplicates/{sample}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "samtools view -@ {threads} -F 1804 -f 2 -L {input.bed} -b "
        "-o {output} {input.bam} 2> {log}"


# Index BAM files with samtools
# -----------------------------------------------------
rule samtools_index:
    input:
        "results/{dir}/{sample}.bam",
    output:
        "results/{dir}/{sample}.bam.bai",
    log:
        "logs/samtools_index/{dir}/{sample}.log",
    wildcard_constraints:
        dir="bowtie2|filtered",
        sample=r"[a-zA-Z0-9_]+(\.dupmark|\.nodup)?",
    params:
        extra="",  # optional params string
    threads: 3  # This value - 1 will be sent to -@
    wrapper:
        "v8.1.1/bio/samtools/index"


# Get BAM file statistics with samtools stats
# (raw alignments and final filtered alignments)
# -----------------------------------------------------
rule samtools_stats:
    input:
        bam="results/{dir}/{sample}.bam",
    output:
        "results/samtools_stats/{dir}/{sample}.txt",
    wildcard_constraints:
        dir="bowtie2|filtered",
        sample=r"[a-zA-Z0-9_]+(\.nodup)?",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "logs/samtools_stats/{dir}/{sample}.log",
    wrapper:
        "v8.1.1/bio/samtools/stats"
