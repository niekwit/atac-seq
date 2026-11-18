# Integration with differential H3K9me3 (optional, h3k9me3: contrasts, with
# deseq2: enable): the differential accessibility of a DESeq2 contrast is
# compared with the H3K9me3 regions that change in the same contrast, from a
# DiffBind results table (log2 fold change and FDR per site) or a BED file
# of regions that change in one direction. The consensus peaks and, with the
# TE module, the TE peaks are compared with the regions; with the TE module
# also the ATAC-seq signal over the regions and the repeat subfamilies. With
# H3K9me3 BigWig files, the H3K9me3 signal around the differential peaks.
# Not part of the ENCODE pipeline.


# Prepared H3K9me3 regions of a contrast: with te: only, in the TE module
# results, so that all output is in results/te/integration/
H3K9ME3_REGIONS = (
    "results/te/integration/h3k9me3/{contrast}/regions.tsv"
    if te_only()
    else "results/integration/h3k9me3/{contrast}/regions.tsv"
)
H3K9ME3_LOST = (
    "results/te/integration/h3k9me3/{contrast}/lost.bed"
    if te_only()
    else "results/integration/h3k9me3/{contrast}/lost.bed"
)
H3K9ME3_GAINED = (
    "results/te/integration/h3k9me3/{contrast}/gained.bed"
    if te_only()
    else "results/integration/h3k9me3/{contrast}/gained.bed"
)
# Copy of the regions file of the contrast (h3k9me3: contrasts: regions)
H3K9ME3_INPUTS = (
    "results/te/integration/h3k9me3/{contrast}/input_tables"
    if te_only()
    else "results/integration/h3k9me3/{contrast}/input_tables"
)
H3K9ME3_PREPARE_LOG = (
    "logs/te/integration/h3k9me3/prepare_{contrast}.log"
    if te_only()
    else "logs/integration/h3k9me3/prepare_{contrast}.log"
)
H3K9ME3_INPUTS_LOG = (
    "logs/te/integration/h3k9me3/input_tables_{contrast}.log"
    if te_only()
    else "logs/integration/h3k9me3/input_tables_{contrast}.log"
)


# Lost, gained and unchanged H3K9me3 regions of the contrast
# -----------------------------------------------------
rule h3k9me3_prepare:
    input:
        regions=lambda w: h3k9me3_contrasts()[w.contrast]["regions"],
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        table=H3K9ME3_REGIONS,
        lost=H3K9ME3_LOST,
        gained=H3K9ME3_GAINED,
    params:
        diffbind=lambda w: h3k9me3_contrasts()[w.contrast]["diffbind"],
        direction=lambda w: h3k9me3_contrasts()[w.contrast]["direction"],
        columns=h3k9me3_config()["columns"],
        flip_fold=h3k9me3_config()["flip_fold"],
        fdr=h3k9me3_config()["fdr"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        H3K9ME3_PREPARE_LOG,
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/h3k9me3_prepare.py"


# Copy of the H3K9me3 regions file of the contrast, unchanged, so that the
# results keep the input they were made from (BigWig files are not copied)
# -----------------------------------------------------
rule h3k9me3_inputs:
    input:
        lambda w: h3k9me3_contrasts()[w.contrast]["regions"],
    output:
        directory(H3K9ME3_INPUTS),
    params:
        files=lambda w: {"regions": h3k9me3_contrasts()[w.contrast]["regions"]},
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        H3K9ME3_INPUTS_LOG,
    script:
        "../scripts/copy_inputs.py"


# The H3K9me3 region each tested consensus peak overlaps most
# -----------------------------------------------------
rule h3k9me3_annotate_peaks:
    input:
        atac="results/deseq2/{contrast}/differential_peaks.tsv",
        regions=H3K9ME3_REGIONS,
    output:
        peaks="results/integration/h3k9me3/{contrast}/peaks.tsv",
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/h3k9me3/annotate_peaks_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/h3k9me3_annotate_peaks.py"


use rule h3k9me3_annotate_peaks as te_h3k9me3_annotate_peaks with:
    input:
        atac="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        regions=H3K9ME3_REGIONS,
    output:
        peaks="results/te/integration/h3k9me3/{contrast}/peaks/peaks.tsv",
    log:
        "logs/te/integration/h3k9me3/annotate_peaks_{contrast}.log",


# Differential peaks in the H3K9me3 regions, and the accessibility changes by
# H3K9me3 change
# -----------------------------------------------------
rule h3k9me3_peaks:
    input:
        peaks="results/integration/h3k9me3/{contrast}/peaks.tsv",
    output:
        enrichment="results/integration/h3k9me3/{contrast}/enrichment.tsv",
        enrichment_plot="results/integration/h3k9me3/{contrast}/enrichment.pdf",
        distribution="results/integration/h3k9me3/{contrast}/distribution.pdf",
        distribution_tests="results/integration/h3k9me3/{contrast}/distribution_tests.tsv",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        changes=lambda w: h3k9me3_changes(w.contrast),
        peak_set="ENCODE consensus peaks (uniquely mapped reads)",
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/h3k9me3/peaks_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/h3k9me3_peaks.R"


use rule h3k9me3_peaks as te_h3k9me3_peaks with:
    input:
        peaks="results/te/integration/h3k9me3/{contrast}/peaks/peaks.tsv",
    output:
        enrichment="results/te/integration/h3k9me3/{contrast}/peaks/enrichment.tsv",
        enrichment_plot="results/te/integration/h3k9me3/{contrast}/peaks/enrichment.pdf",
        distribution="results/te/integration/h3k9me3/{contrast}/peaks/distribution.pdf",
        distribution_tests="results/te/integration/h3k9me3/{contrast}/peaks/distribution_tests.tsv",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        changes=lambda w: h3k9me3_changes(w.contrast),
        peak_set="TE module peaks (alignment keeping multimapping reads)",
    log:
        "logs/te/integration/h3k9me3/peaks_{contrast}.log",


# DiffBind data: ATAC-seq against H3K9me3 log2 fold changes of the peaks in
# DiffBind sites
# -----------------------------------------------------
rule h3k9me3_scatter:
    input:
        peaks="results/integration/h3k9me3/{contrast}/peaks.tsv",
    output:
        correlation="results/integration/h3k9me3/{contrast}/correlation.tsv",
        scatter="results/integration/h3k9me3/{contrast}/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        peak_set="ENCODE consensus peaks",
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/h3k9me3/scatter_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/h3k9me3_scatter.R"


use rule h3k9me3_scatter as te_h3k9me3_scatter with:
    input:
        peaks="results/te/integration/h3k9me3/{contrast}/peaks/peaks.tsv",
    output:
        correlation="results/te/integration/h3k9me3/{contrast}/peaks/correlation.tsv",
        scatter="results/te/integration/h3k9me3/{contrast}/peaks/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        peak_set="TE module peaks",
    log:
        "logs/te/integration/h3k9me3/scatter_{contrast}.log",


# H3K9me3 signal (h3k9me3: contrasts: bigwigs) of the control and the
# treatment around the centre of the more, less and unchanged accessible peaks
# -----------------------------------------------------
rule h3k9me3_signal_profile:
    input:
        table="results/deseq2/{contrast}/differential_peaks.tsv",
        # Control first, as the reference (grey)
        bigwigs=lambda w: [
            h3k9me3_contrasts()[w.contrast]["bigwigs"]["control"],
            h3k9me3_contrasts()[w.contrast]["bigwigs"]["treatment"],
        ],
    output:
        matrix="results/integration/h3k9me3/{contrast}/signal_profile.matrix.gz",
        profile="results/integration/h3k9me3/{contrast}/signal_profile.pdf",
        heatmap="results/integration/h3k9me3/{contrast}/signal_heatmap.pdf",
    params:
        group_column="change",
        groups={
            "more accessible": "more accessible",
            "less accessible": "less accessible",
            "unchanged": "unchanged",
        },
        mode="reference-point",
        flank=h3k9me3_config()["profile_flank"],
        labels=lambda w: h3k9me3_conditions(w.contrast),
        colours=["#999999", "#dd3b3b"],
        y_label="H3K9me3 signal",
        region_name="peaks",
        extra_computematrix=extra("computematrix"),
        extra_plotprofile=extra("plotprofile"),
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=8000,
    log:
        "logs/integration/h3k9me3/signal_profile_{contrast}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/h3k9me3_profile.py"


use rule h3k9me3_signal_profile as te_h3k9me3_signal_profile with:
    input:
        table="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        bigwigs=lambda w: [
            h3k9me3_contrasts()[w.contrast]["bigwigs"]["control"],
            h3k9me3_contrasts()[w.contrast]["bigwigs"]["treatment"],
        ],
    output:
        matrix="results/te/integration/h3k9me3/{contrast}/peaks/signal_profile.matrix.gz",
        profile="results/te/integration/h3k9me3/{contrast}/peaks/signal_profile.pdf",
        heatmap="results/te/integration/h3k9me3/{contrast}/peaks/signal_heatmap.pdf",
    log:
        "logs/te/integration/h3k9me3/signal_profile_{contrast}.log",


# Tn5 insertions (CPM, mean of the replicates) of the control and the
# treatment over the lost, gained and unchanged H3K9me3 regions (TE module,
# which has the CPM signal tracks)
# -----------------------------------------------------
use rule h3k9me3_signal_profile as te_h3k9me3_atac_profile with:
    input:
        table=H3K9ME3_REGIONS,
        bigwigs=lambda w: expand(
            "results/te/bigwig/mean/{condition}.cpm.bigwig",
            condition=h3k9me3_conditions(w.contrast),
        ),
    output:
        matrix="results/te/integration/h3k9me3/{contrast}/atac_profile.matrix.gz",
        profile="results/te/integration/h3k9me3/{contrast}/atac_profile.pdf",
        heatmap="results/te/integration/h3k9me3/{contrast}/atac_heatmap.pdf",
    params:
        group_column="change",
        groups=lambda w: {
            "lost": "H3K9me3 lost",
            "gained": "H3K9me3 gained",
            "unchanged": "H3K9me3 unchanged",
        },
        mode="scale-regions",
        flank=h3k9me3_config()["profile_flank"],
        body_length=h3k9me3_config()["profile_body_length"],
        labels=lambda w: h3k9me3_conditions(w.contrast),
        colours=["#999999", "#dd3b3b"],
        y_label="Tn5 insertions per million",
        region_name="H3K9me3 regions",
        extra_computematrix=extra("computematrix"),
        extra_plotprofile=extra("plotprofile"),
    log:
        "logs/te/integration/h3k9me3/atac_profile_{contrast}.log",


# Repeat copies in the lost and gained H3K9me3 regions
# -----------------------------------------------------
rule te_h3k9me3_copies:
    input:
        saf="resources/te/repeats.saf",
        regions=H3K9ME3_REGIONS,
        # Only with te: profile: subfamilies (otherwise no full-length copies)
        full_length=lambda w: (
            ["resources/te/full_length_copies.bed"]
            if te_profile_config()["subfamilies"]
            else []
        ),
        peaks="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        copies="results/te/integration/h3k9me3/{contrast}/subfamilies/copies.tsv",
        full_length="results/te/integration/h3k9me3/{contrast}/subfamilies/full_length_copies.tsv",
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=8000,
        cpus_per_task=cpus_for_memory(8000),
    log:
        "logs/te/integration/h3k9me3/copies_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/h3k9me3_te_copies.py"


# Repeat subfamilies in the H3K9me3 regions against their differential
# accessibility, and the full-length copies of the profile subfamilies
# -----------------------------------------------------
rule te_h3k9me3_subfamilies:
    input:
        copies="results/te/integration/h3k9me3/{contrast}/subfamilies/copies.tsv",
        full_length="results/te/integration/h3k9me3/{contrast}/subfamilies/full_length_copies.tsv",
        atac="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
    output:
        subfamilies="results/te/integration/h3k9me3/{contrast}/subfamilies/subfamilies.tsv",
        correlation="results/te/integration/h3k9me3/{contrast}/subfamilies/correlation.tsv",
        scatter="results/te/integration/h3k9me3/{contrast}/subfamilies/scatter.pdf",
        full_length_tests="results/te/integration/h3k9me3/{contrast}/subfamilies/full_length.tsv",
        full_length="results/te/integration/h3k9me3/{contrast}/subfamilies/full_length.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        changes=lambda w: h3k9me3_changes(w.contrast),
        min_copies=h3k9me3_config()["min_copies"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=3000,
    log:
        "logs/te/integration/h3k9me3/subfamilies_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/h3k9me3_subfamilies.R"
