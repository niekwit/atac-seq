# Peaks of the TE module (te: enable): peaks called from the alignments that
# keep multimappers, so that open regions in young repeats are included.
# Same steps as the main peaks where possible (Tn5-shifted tagAlign, MACS2,
# peaks replicated in min_replicate_fraction of the replicates, consensus
# peaks, DESeq2), and with te: correlation: datasets, GAT enrichment of the
# differential peaks in other data sets (e.g. EM-seq DMRs).


# Tn5 shifted tagAlign of the TE alignments
# -----------------------------------------------------
use rule bam2ta as te_bam2ta with:
    input:
        bam="results/te/{sample}.filtered.bam",
    output:
        "results/te/tagalign/{sample}.tagAlign.gz",
    log:
        "logs/te/bam2ta/{sample}.log",


use rule pool_tagalign as te_pool_tagalign with:
    input:
        lambda w: expand(
            "results/te/tagalign/{sample}.tagAlign.gz", sample=replicates(w.condition)
        ),
    output:
        "results/te/tagalign/pooled/{condition}.tagAlign.gz",
    log:
        "logs/te/pool_tagalign/{condition}.log",


# MACS2 as for the main peaks; multimappers can align anywhere, so the
# effective genome size is the total length of the chromosomes
# -----------------------------------------------------
use rule macs2_callpeak as te_macs2_callpeak with:
    input:
        ta="results/te/tagalign/{prefix}.tagAlign.gz",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        narrowpeak="results/te/macs2/{prefix}.narrowPeak",
        treat=temp("results/te/macs2/{prefix}_treat_pileup.bdg"),
        control=temp("results/te/macs2/{prefix}_control_lambda.bdg"),
    params:
        gsize="chrom_sizes",
        pval_thresh=config["macs2"]["pval_thresh"],
        smooth_win=config["macs2"]["smooth_win"],
        cap_num_peak=config["macs2"]["cap_num_peak"],
    log:
        "logs/te/macs2/{prefix}.log",


use rule blacklist_filter as te_blacklist_filter with:
    wildcard_constraints:
        stem=r"results/te/macs2/pooled/[a-zA-Z0-9_]+",


# Peaks of the pooled replicates that overlap a peak of at least
# min_replicate_fraction of the replicates by >= 50%
# -----------------------------------------------------
use rule replicated_peaks as te_replicated_peaks with:
    input:
        pooled="results/te/macs2/pooled/{condition}.bfilt.narrowPeak",
        replicates=lambda w: expand(
            "results/te/macs2/{sample}.narrowPeak", sample=replicates(w.condition)
        ),
        pairs=[],
    output:
        "results/te/peaks/{method}/{condition}/replicated_peak.narrowPeak",
    params:
        pairs=[],
        min_fraction=config["peaks"].get("min_replicate_fraction", 1),
    wildcard_constraints:
        method="overlap",
    log:
        "logs/te/replicated_peaks/{method}/{condition}.log",


# Consensus peaks: the merged replicated peaks of all conditions
# -----------------------------------------------------
use rule consensus_peaks as te_consensus_peaks with:
    input:
        peaks=expand(
            "results/te/peaks/overlap/{condition}/replicated_peak.narrowPeak",
            condition=CONDITIONS,
        ),
    output:
        bed="results/te/peaks/consensus_peaks.bed",
        saf="results/te/peaks/consensus_peaks.saf",
    log:
        "logs/te/peaks/consensus_peaks.log",


# Tn5 insertions in the consensus peaks: all read pairs (-M) and uniquely
# mapping read pairs only (NH:i:1)
# -----------------------------------------------------
rule te_count_peak_insertions:
    input:
        bams=expand("results/te/{sample}.filtered.bam", sample=SAMPLES),
        saf="results/te/peaks/consensus_peaks.saf",
    output:
        multi="results/te/peaks/consensus_peak_counts.txt",
        unique="results/te/peaks/consensus_peak_counts_unique.txt",
    params:
        # -p without --countReadPairs: each read (Tn5 insertion) is counted
        options="-p -F SAF --read2pos 5 --largestOverlap -s 0",
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=8000,
    log:
        "logs/te/peaks/count_peak_insertions.log",
    conda:
        "../envs/deseq2.yaml"
    shell:
        "(featureCounts {params.options} -M -T {threads} -a {input.saf} "
        "-o {output.multi} {input.bams} && "
        "featureCounts {params.options} -T {threads} -a {input.saf} "
        "-o {output.unique} {input.bams}) > {log} 2>&1"


# Repeat with the largest overlap and fraction of insertions from uniquely
# mapping reads of each consensus peak
# -----------------------------------------------------
rule te_peak_annotation:
    input:
        peaks="results/te/peaks/consensus_peaks.bed",
        multi="results/te/peaks/consensus_peak_counts.txt",
        unique="results/te/peaks/consensus_peak_counts_unique.txt",
        saf="resources/te/repeats.saf",
        subfamilies="resources/te/subfamilies.tsv",
    output:
        table="results/te/peaks/consensus_peaks_repeats.tsv",
        bed="results/te/peaks/consensus_peaks_frac_unique.bed",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=3000,
    log:
        "logs/te/peaks/peak_annotation.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/te_peak_annotation.py"


use rule annotate_consensus_peaks as te_annotate_consensus_peaks with:
    input:
        peaks="results/te/peaks/consensus_peaks.bed",
        edb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        txdb=f"resources/{resources.genome}_{resources.build}_txdb.Rdata",
    output:
        txt="results/te/peaks/consensus_peaks_annotated.txt",
    log:
        "logs/te/peaks/annotate_consensus_peaks.log",


# Differential accessibility of the TE peaks (with deseq2: enable)
# -----------------------------------------------------
use rule deseq2_fit as te_peaks_deseq2_fit with:
    input:
        counts="results/te/peaks/consensus_peak_counts.txt",
    output:
        dds="results/te/peaks/deseq2/dds.rds",
    log:
        "logs/te/peaks/deseq2_fit.log",


use rule deseq2_contrast as te_peaks_deseq2_contrast with:
    input:
        dds="results/te/peaks/deseq2/dds.rds",
        annotation="results/te/peaks/consensus_peaks_annotated.txt",
        peaks="results/te/peaks/consensus_peaks.bed",
        extra="results/te/peaks/consensus_peaks_repeats.tsv",
    output:
        table="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        ma="results/te/peaks/deseq2/{contrast}/ma_plot.png",
        volcano="results/te/peaks/deseq2/{contrast}/volcano_plot.png",
    log:
        "logs/te/peaks/deseq2_{contrast}.log",


# Differential peaks as BED6 (name: peak_id, score: 0; the fold changes are in
# differential_peaks.tsv)
# -----------------------------------------------------
rule te_differential_peak_beds:
    input:
        "results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
    output:
        more="results/te/peaks/deseq2/{contrast}/more_accessible.bed",
        less="results/te/peaks/deseq2/{contrast}/less_accessible.bed",
    params:
        awk=(
            'BEGIN{FS=OFS="\\t"} NR==1{for (i=1; i<=NF; i++) c[$i]=i; next} '
            '$c["change"]==change{print $c["chr"], $c["start"], $c["end"], '
            '$c["peak_id"], 0, "."}'
        ),
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        "logs/te/peaks/differential_peak_beds_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    shell:
        "(awk -v change='more accessible' '{params.awk}' {input} | "
        "sort -k1,1 -k2,2n > {output.more} && "
        "awk -v change='less accessible' '{params.awk}' {input} | "
        "sort -k1,1 -k2,2n > {output.less}) 2> {log}"


# GAT enrichment of the differential peaks in the regions of other data sets
# (te: correlation: datasets), within the TE consensus peaks
# -----------------------------------------------------
rule te_gat:
    input:
        more="results/te/peaks/deseq2/{contrast}/more_accessible.bed",
        less="results/te/peaks/deseq2/{contrast}/less_accessible.bed",
        workspace="results/te/peaks/consensus_peaks.bed",
        chrom_sizes="resources/chrom_sizes.txt",
        datasets=lambda w: list(te_correlation_config()["datasets"].values()),
    output:
        table="results/te/correlation/{contrast}/gat.tsv",
    params:
        # Functions, so the data sets are checked only when GAT runs
        datasets=lambda w: te_correlation_config()["datasets"],
        samples=lambda w: te_correlation_config()["gat_samples"],
    threads: 4
    resources:
        runtime=runtime(120),
        mem_mb=12000,
    log:
        "logs/te/correlation/gat_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/te_gat.py"


# Area-proportional Venn diagrams of the differential peaks and each data set,
# with the GAT results
# -----------------------------------------------------
rule te_venn:
    input:
        table="results/te/correlation/{contrast}/gat.tsv",
    output:
        venn="results/te/correlation/{contrast}/venn.png",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=2000,
    log:
        "logs/te/correlation/venn_{contrast}.log",
    conda:
        "../envs/eulerr.yaml"
    script:
        "../scripts/te_venn.R"
