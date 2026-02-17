# Identifying peaks in ATAC-seq data using MACS2
# -----------------------------------------------------
rule callpeak:
    input:
        treatment="results/dedup/{sample}.bam",
        bai="results/dedup/{sample}.bam.bai",
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/macs2/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    threads: 1
    log:
        "logs/macs2/{sample}.log"
    params:
        config["macs2"]["extra"],
    wrapper:
        "v2.9.1/bio/macs2/callpeak"


# Remove blacklisted regions from peak files
# -----------------------------------------------------
rule remove_blacklisted_regions:
    input:
        left="results/macs2/{sample}_peaks.narrowPeak",
        right="resources/blacklist.bed"
    output:
        "results/macs2/{sample}.no_blacklist.narrowPeak"
    params:
        extra="-v"
    log:
        "logs/blacklist_removal/{sample}.log"
    wrapper:
        "v8.1.1/bio/bedtools/intersect"


# QC of ATAC-seq data using ataqv
# -----------------------------------------------------
rule ataqv:
    input:
        bam="results/dedup/{sample}.bam", # use deduplicated BAM
        bai="results/dedup/{sample}.bam.bai",
        macs2="results/macs2/{sample}_peaks.narrowPeak", # raw peaks
        tss="resources/tss.bed",
        blacklist="resources/blacklist.bed"
    output:
        json="results/ataqv/{sample}.json.gz",
        out="results/ataqv/{sample}.out"
    params:
        organism=ataqv_organism(),
        extra=""
    log:
        "logs/ataqv/{sample}.log"
    threads: 4
    conda:
        "../envs/atac.yaml"
    shell:
        "ataqv {params.organism} "
        "--threads {threads} "
        "--mitochondrial-reference-name MT "
        "--name {wildcards.sample} "
        "--peak-file {input.macs2} "
        "--tss-file {input.tss} "
        "--metrics-file {output.json} "
        "--excluded-region-file {input.blacklist} "
        "{params.extra} "
        "{input.bam} > {output.out} 2> {log}"


# Create HTML report of ATAC-seq QC metrics using ataqv mkarv
# -----------------------------------------------------
rule ataqv_report:
    input:
        json=expand("results/ataqv/{sample}.json.gz", sample=SAMPLES),
    output:
        html="results/ataqv_report/index.html",
    params:
        dir=lambda w, output: os.path.dirname(output["html"]),
        extra=""
    log:
        "logs/ataqv/ataqv_report.log"
    threads: 1
    conda:
        "../envs/atac.yaml"
    shell:
        "mkarv --force {params.dir} {input.json} 2> {log}"


# Generate consensus peak set for each condition
# -----------------------------------------------------
rule consensus_peaks:
    input:
        peaks=expand("results/macs2/{sample}.no_blacklist.narrowPeak", sample=SAMPLES),
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        bed="results/macs2/{condition}_consensus_peaks.bed",
    params:
        method=config["consensus_peaks"]["method"],
        keep=config["consensus_peaks"]["multiinter_options"]["keep"],
    log:
        "logs/consensus_peaks/{condition}.log"
    script:
        "../scripts/consensus_peaks.py"