# Generate bigwig files from filtered BAM files using deepTools bamCoverage
# -----------------------------------------------------
rule bigwig:
    input:
        bam="results/filtered/{sample}.bam",
        bai="results/filtered/{sample}.bam.bai",
    output:
        "results/bigwig/samples/{sample}.bw",
    params:
        genome=bamcoverage_genome(),
        read_length=config["read_length"],
        extra="",
    log:
        "logs/deeptools/bamcoverage_{sample}.log",
    wrapper:
        "v5.6.0/bio/deeptools/bamcoverage"


# Create intermediate wig files by averaging bigwig files from samples in the same condition
# -----------------------------------------------------
rule average_wig:
    input:
        expand("results/bigwig/samples/{sample}.bw", sample=SAMPLES),
    output:
        wig=temp("results/bigwig/{condition}.wig"),
    threads: 1
    log:
        "logs/wiggletools/wig_average_{condition}.log",
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/average_wig.py"


# Convert wig files to bigwig files
# -----------------------------------------------------
rule wig2bigwig:
    input:
        wig="results/bigwig/{condition}.wig",
        cs="resources/chrom_sizes.txt",
    output:
        "results/bigwig/average/{condition}.bw",
    params:
        extra="",
    threads: 1
    log:
        "logs/wigToBigWig/{condition}.log",
    conda:
        "../envs/atac.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"