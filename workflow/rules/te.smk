# Transposable element (TE) module (optional, te: enable)
# Accessibility of repeat subfamilies, families and classes. Young repeats
# have many near-identical copies, so their reads are removed by the ENCODE
# multimapper filter. Here, reads are aligned again (bowtie2 default mode:
# one of the best alignments, picked at random for multimappers), and Tn5
# insertions in RepeatMasker repeats are counted per subfamily, family and
# class. See docs/te_alignment_comparison for the comparison with -k 100.
# The ENCODE outputs are not affected.


# Download the RepeatMasker annotation of the genome from UCSC
# -----------------------------------------------------
rule get_rmsk:
    output:
        "resources/te/rmsk.txt.gz",
    retries: 3
    params:
        url=(
            "https://hgdownload.soe.ucsc.edu/goldenPath/"
            f"{UCSC_GENOMES[config['genome']['ensembl']]}/database/rmsk.txt.gz"
        ),
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=1000,
    log:
        "logs/te/get_rmsk.log",
    conda:
        "../envs/te.yaml"
    shell:
        "wget -q {params.url} -O {output} 2> {log}"


# Repeat annotation (SAF) and subfamily table, with Ensembl chromosome names
# -----------------------------------------------------
rule prepare_te_annotation:
    input:
        rmsk="resources/te/rmsk.txt.gz",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        saf="resources/te/repeats.saf",
        subfamilies="resources/te/subfamilies.tsv",
    params:
        classes=config.get("te", {}).get("classes", []),
        filter_chrs=config["filter"]["filter_chrs"],
    threads: 2
    resources:
        runtime=runtime(30),
        mem_mb=4000,
    log:
        "logs/te/prepare_te_annotation.log",
    conda:
        "../envs/te.yaml"
    script:
        "../scripts/prepare_te_annotation.py"


# Alignment keeping multimappers (bowtie2 default mode, reproducible random
# choice among equally good alignments), in read name order
# -----------------------------------------------------
rule te_align:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        idx=rules.bowtie2_index.output,
    output:
        temp("results/te/{sample}.bowtie2.bam"),
    params:
        prefix=lambda w, input: input.idx[0].removesuffix(".1.bt2"),
        extra=config["bowtie2"]["extra"],
    threads: 8
    resources:
        runtime=runtime(360),
        mem_mb=8000,
    log:
        "logs/te/te_align_{sample}.log",
    conda:
        "../envs/te.yaml"
    shell:
        "(bowtie2 --seed 42 {params.extra} --threads {threads} "
        "-x {params.prefix} -1 {input.r1} -2 {input.r2} | "
        "samtools view -b -@ 2 -o {output} -) 2> {log}"


# Keep properly paired alignments, remove duplicates and reads on filter_chrs
# and mark multimappers (NH:i:2)
# -----------------------------------------------------
rule te_filter:
    input:
        "results/te/{sample}.bowtie2.bam",
    output:
        bam=temp("results/te/{sample}.filtered.bam"),
        stats="results/te/{sample}.filter_stats.tsv",
    params:
        filter_chrs=config["filter"]["filter_chrs"],
    threads: 4
    resources:
        runtime=runtime(240),
        mem_mb=8000,
    log:
        "logs/te/te_filter_{sample}.log",
    conda:
        "../envs/te.yaml"
    script:
        "../scripts/te_filter.py"


# Count Tn5 insertions (5' end of each read) in repeats: all reads (-M) and
# uniquely mapping read pairs only (NH:i:1)
# -----------------------------------------------------
rule te_count:
    input:
        bams=expand("results/te/{sample}.filtered.bam", sample=SAMPLES),
        saf="resources/te/repeats.saf",
    output:
        multi="results/te/featurecounts_multi.txt",
        unique="results/te/featurecounts_unique.txt",
    params:
        # -p without --countReadPairs: each read (Tn5 insertion) is counted
        options="-p -F SAF --read2pos 5 --largestOverlap -s 0",
    threads: 8
    resources:
        runtime=runtime(120),
        mem_mb=8000,
    log:
        "logs/te/te_count.log",
    conda:
        "../envs/te.yaml"
    shell:
        "(featureCounts {params.options} -T {threads} -M "
        "-a {input.saf} -o {output.multi} {input.bams} && "
        "featureCounts {params.options} -T {threads} "
        "-a {input.saf} -o {output.unique} {input.bams}) > {log} 2>&1"


# Counts, CPM and enrichment of subfamilies, families and classes
# -----------------------------------------------------
rule te_summary:
    input:
        multi="results/te/featurecounts_multi.txt",
        unique="results/te/featurecounts_unique.txt",
        subfamilies="resources/te/subfamilies.tsv",
        stats=expand("results/te/{sample}.filter_stats.tsv", sample=SAMPLES),
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        subfamily="results/te/te_subfamily.tsv",
        family="results/te/te_family.tsv",
        repeat_class="results/te/te_class.tsv",
        stats="results/te/te_library_stats.tsv",
    params:
        samples=SAMPLES,
        filter_chrs=config["filter"]["filter_chrs"],
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=2000,
    log:
        "logs/te/te_summary.log",
    conda:
        "../envs/te.yaml"
    script:
        "../scripts/te_summary.py"


# Fit the DESeq2 model (~ condition) to the subfamily counts of all samples
# (with deseq2: enable)
# -----------------------------------------------------
rule te_deseq2_fit:
    input:
        counts="results/te/featurecounts_multi.txt",
        stats="results/te/te_library_stats.tsv",
    output:
        dds="results/te/deseq2/dds.rds",
    params:
        samples=SAMPLES,
        conditions=SAMPLE_SHEET["condition"].tolist(),
        levels=deseq2_levels(),
        normalisation=config.get("te", {}).get("normalisation", "deseq2"),
    threads: 2
    resources:
        runtime=runtime(30),
        mem_mb=4000,
    log:
        "logs/te/te_deseq2_fit.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/te_deseq2_fit.R"


# Differential accessibility of subfamilies between two conditions
# -----------------------------------------------------
rule te_deseq2_contrast:
    input:
        dds="results/te/deseq2/dds.rds",
        subfamilies="resources/te/subfamilies.tsv",
    output:
        table="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
        ma="results/te/deseq2/{contrast}/ma_plot.png",
        volcano="results/te/deseq2/{contrast}/volcano_plot.png",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        fdr=config.get("deseq2", {}).get("fdr", 0.05),
        lfc_threshold=config.get("deseq2", {}).get("lfc_threshold", 0),
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=3000,
    log:
        "logs/te/te_deseq2_{contrast}.log",
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/te_deseq2_contrast.R"
