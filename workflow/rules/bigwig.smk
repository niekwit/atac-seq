# Generate bigwig files from filtered BAM files using deepTools bamCoverage
# -----------------------------------------------------
rule bigwig:
    input:
        bam="results/filtered/{sample}.bam",
        bai="results/filtered/{sample}.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        genome=bamcoverage_genome(),
        read_length=config["read_length"],
        extra="",
    log:
        "logs/deeptools/bamcoverage_{sample}.log",
    wrapper:
        "v5.6.0/bio/deeptools/bamcoverage"


rule average_wig:
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
    output:
        wig=temp("results/bigwig/{condition}.wig"),
    threads: 2
    log:
        "logs/wiggletools/wig_average_{condition}.log",
    conda:
        "../envs/atac.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig="results/bigwig/{condition}.wig",
        cs="resources/chrom.sizes.txt",
    output:
        "results/bigwig/{condition}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/wigToBigWig/{condition}.log",
    conda:
        "../envs/atac.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"