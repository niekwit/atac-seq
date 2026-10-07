# Peak files that are blacklist filtered and used for FRiP
STEM_REGEX = (
    r"results/macs2/(pooled/)?[a-zA-Z0-9_]+(\.pr[12])?"
    r"|results/(idr|overlap)/[a-zA-Z0-9_]+/"
    r"(rep\d+_vs_rep\d+|rep\d+-pr1_vs_rep\d+-pr2|pooled-pr1_vs_pooled-pr2)"
)


# Convert filtered BAM files to Tn5 shifted tagAlign files
# -----------------------------------------------------
rule bam2ta:
    input:
        bam="results/filtered/{sample}.nodup.bam",
    output:
        "results/tagalign/{sample}.tagAlign.gz",
    params:
        subsample=config["subsample_reads"],
    threads: 4
    log:
        "logs/bam2ta/{sample}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/bam2ta.sh"


# Split each replicate into two pseudoreplicates
# -----------------------------------------------------
rule pseudoreplicates:
    input:
        "results/tagalign/{sample}.tagAlign.gz",
    output:
        pr1="results/tagalign/{sample}.pr1.tagAlign.gz",
        pr2="results/tagalign/{sample}.pr2.tagAlign.gz",
    params:
        seed=config["pseudoreplication_random_seed"],
    threads: 1
    log:
        "logs/pseudoreplicates/{sample}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/spr.sh"


# Pool the replicates (and their pseudoreplicates) of each condition
# -----------------------------------------------------
rule pool_tagalign:
    input:
        lambda w: expand(
            "results/tagalign/{sample}.tagAlign.gz", sample=replicates(w.condition)
        ),
    output:
        "results/tagalign/pooled/{condition}.tagAlign.gz",
    threads: 1
    log:
        "logs/pool_tagalign/{condition}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "zcat -f {input} | gzip -nc > {output} 2> {log}"


use rule pool_tagalign as pool_pseudoreplicates with:
    input:
        lambda w: expand(
            "results/tagalign/{sample}.{pr}.tagAlign.gz",
            sample=replicates(w.condition),
            pr=w.pr,
        ),
    output:
        "results/tagalign/pooled/{condition}.{pr}.tagAlign.gz",
    wildcard_constraints:
        pr="pr1|pr2",
    log:
        "logs/pool_tagalign/{condition}.{pr}.log",


# Call peaks with MACS2
# -----------------------------------------------------
rule macs2_callpeak:
    input:
        ta="results/tagalign/{prefix}.tagAlign.gz",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        narrowpeak="results/macs2/{prefix}.narrowPeak",
        treat=temp("results/macs2/{prefix}_treat_pileup.bdg"),
        control=temp("results/macs2/{prefix}_control_lambda.bdg"),
    params:
        gsize=config["macs2"]["gsize"],
        pval_thresh=config["macs2"]["pval_thresh"],
        smooth_win=config["macs2"]["smooth_win"],
        cap_num_peak=config["macs2"]["cap_num_peak"],
    threads: 1
    log:
        "logs/macs2/{prefix}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/macs2_callpeak.sh"


# IDR on pairs of replicates/pseudoreplicates
# -----------------------------------------------------
rule idr:
    input:
        unpack(pair_peaks),
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        narrowpeak="results/idr/{condition}/{pair}.narrowPeak",
        unthresholded="results/idr/{condition}/{pair}.unthresholded-peaks.txt.gz",
        plot="results/idr/{condition}/{pair}.png",
        idr_log="results/idr/{condition}/{pair}.idr.log",
    params:
        threshold=config["idr"]["threshold"],
    threads: 1
    log:
        "logs/idr/{condition}/{pair}.log",
    conda:
        "../envs/idr.yaml"
    script:
        "../scripts/idr.sh"


# Naive overlap of pairs of replicates/pseudoreplicates: pooled peaks
# that overlap a peak in both sets by >= 50% (of either peak)
# -----------------------------------------------------
rule overlap:
    input:
        unpack(pair_peaks),
    output:
        "results/overlap/{condition}/{pair}.narrowPeak",
    params:
        awk=r"{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}",
    threads: 1
    log:
        "logs/overlap/{condition}/{pair}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(intersectBed -nonamecheck -wo -a {input.pooled} -b {input.peak1} | "
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {params.awk}\' | '
        "cut -f 1-10 | sort | uniq | "
        "intersectBed -nonamecheck -wo -a stdin -b {input.peak2} | "
        'awk \'BEGIN{{FS="\\t";OFS="\\t"}} {params.awk}\' | '
        "cut -f 1-10 | sort | uniq > {output}) 2> {log}"


# Remove peaks in blacklisted regions and on non-standard chromosomes
# -----------------------------------------------------
rule blacklist_filter:
    input:
        peak="{stem}.narrowPeak",
        blacklist="resources/blacklist.bed",
    output:
        "{stem}.bfilt.narrowPeak",
    wildcard_constraints:
        stem=STEM_REGEX,
    params:
        regex=config["peaks"]["keep_chr_regex"],
    threads: 1
    log:
        "logs/blacklist_filter/{stem}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(bedtools intersect -nonamecheck -v -a {input.peak} -b {input.blacklist} | "
        "awk 'BEGIN{{OFS=\"\\t\"}} {{if ($5>1000) $5=1000; print $0}}' | "
        "awk -v re='^({params.regex})$' '$1 ~ re' > {output}) 2> {log}"


# Fraction of reads in peaks
# -----------------------------------------------------
rule frip:
    input:
        peak="{stem}.bfilt.narrowPeak",
        ta=frip_tagalign,
    output:
        "{stem}.frip.qc",
    wildcard_constraints:
        stem=STEM_REGEX,
    threads: 1
    log:
        "logs/frip/{stem}.log",
    conda:
        "../envs/encode.yaml"
    shell:
        "(in_peaks=$(bedtools intersect -nonamecheck -a {input.ta} -b {input.peak} -wa -u | wc -l); "
        "total=$(zcat -f {input.ta} | wc -l); "
        "awk -v a=$in_peaks -v b=$total 'BEGIN{{print a / b}}' > {output}) 2> {log}"


# Select optimal and conservative peak sets and check reproducibility
# -----------------------------------------------------
rule reproducibility:
    input:
        peaks=lambda w: expand(
            "results/{method}/{condition}/{pair}.bfilt.narrowPeak",
            method=w.method,
            condition=w.condition,
            pair=[p for p in peak_pairs(w.condition) if re.match(r"rep\d+_vs_", p)],
        ),
        peaks_pr=lambda w: expand(
            "results/{method}/{condition}/{pair}.bfilt.narrowPeak",
            method=w.method,
            condition=w.condition,
            pair=[p for p in peak_pairs(w.condition) if re.match(r"rep\d+-pr1_", p)],
        ),
        peak_ppr=lambda w: expand(
            "results/{method}/{condition}/{pair}.bfilt.narrowPeak",
            method=w.method,
            condition=w.condition,
            pair=[p for p in peak_pairs(w.condition) if p.startswith("pooled")],
        ),
    output:
        optimal="results/{method}/{condition}/optimal_peak.narrowPeak",
        conservative="results/{method}/{condition}/conservative_peak.narrowPeak",
        qc="results/{method}/{condition}/reproducibility.qc",
    params:
        pairs=lambda w: peak_pairs(w.condition),
    threads: 1
    log:
        "logs/reproducibility/{method}/{condition}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/reproducibility.py"


# QC of ATAC-seq data using ataqv
# -----------------------------------------------------
rule ataqv:
    input:
        bam="results/filtered/{sample}.dupmark.bam",  # duplicates marked, not removed
        bai="results/filtered/{sample}.dupmark.bam.bai",
        peaks="results/macs2/{sample}.bfilt.narrowPeak",
        tss="resources/tss.bed",
        blacklist="resources/blacklist.bed",
    output:
        json="results/ataqv/{sample}.json.gz",
        out="results/ataqv/{sample}.out",
    params:
        organism=ataqv_organism(),
        mito=config["mito_chr_name"],
        extra="",
    log:
        "logs/ataqv/{sample}.log",
    threads: 4
    conda:
        "../envs/atac.yaml"
    shell:
        "ataqv {params.organism} "
        "--threads {threads} "
        "--mitochondrial-reference-name {params.mito} "
        "--name {wildcards.sample} "
        "--peak-file {input.peaks} "
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
        extra="",
    log:
        "logs/ataqv/ataqv_report.log",
    threads: 1
    conda:
        "../envs/atac.yaml"
    shell:
        "mkarv --force {params.dir} {input.json} 2> {log}"


# Annotate the optimal reproducible peak set with ChIPseeker
# -----------------------------------------------------
rule annotate_peaks:
    input:
        peaks="results/{method}/{condition}/optimal_peak.narrowPeak",
        edb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        txdb=f"resources/{resources.genome}_{resources.build}_txdb.Rdata",
    output:
        txt="results/{method}/{condition}/optimal_peak_annotated.txt",
    threads: 2
    log:
        "logs/annotate_peaks/{method}/{condition}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/annotate_peaks.R"
