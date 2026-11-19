# Integration with differential epigenomic data (optional, with deseq2:
# enable): the differential accessibility of a DESeq2 contrast is compared
# with the regions where a mark ({mark}: h3k9me3 or methylation, sections of
# the config) changes in the same contrast, from a results table (DiffBind,
# methylKit: effect size and FDR per region) or a BED file of regions that
# change in one direction. The consensus peaks and, with the TE module, the
# TE peaks are compared with the regions; with the TE module also the
# ATAC-seq signal over the regions and the repeat subfamilies. With BigWig
# files of the mark, its signal around the differential peaks. Not part of
# the ENCODE pipeline.


wildcard_constraints:
    mark="|".join(EPI_MARKS),


# Prepared regions of a contrast: with te: only, in the TE module results, so
# that all output is in results/te/integration/
EPI_REGIONS = (
    "results/te/integration/{mark}/{contrast}/regions.tsv"
    if te_only()
    else "results/integration/{mark}/{contrast}/regions.tsv"
)
EPI_DECREASED = (
    "results/te/integration/{mark}/{contrast}/decreased.bed"
    if te_only()
    else "results/integration/{mark}/{contrast}/decreased.bed"
)
EPI_INCREASED = (
    "results/te/integration/{mark}/{contrast}/increased.bed"
    if te_only()
    else "results/integration/{mark}/{contrast}/increased.bed"
)
# Copy of the regions file of the contrast ({mark}: contrasts: regions)
EPI_INPUTS = (
    "results/te/integration/{mark}/{contrast}/input_tables"
    if te_only()
    else "results/integration/{mark}/{contrast}/input_tables"
)
EPI_PREPARE_LOG = (
    "logs/te/integration/{mark}/prepare_{contrast}.log"
    if te_only()
    else "logs/integration/{mark}/prepare_{contrast}.log"
)
EPI_INPUTS_LOG = (
    "logs/te/integration/{mark}/input_tables_{contrast}.log"
    if te_only()
    else "logs/integration/{mark}/input_tables_{contrast}.log"
)


# Regions where the mark decreases, increases or does not change
# -----------------------------------------------------
rule epi_prepare:
    input:
        regions=lambda w: epi_contrasts(w.mark)[w.contrast]["regions"],
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        table=EPI_REGIONS,
        decreased=EPI_DECREASED,
        increased=EPI_INCREASED,
    params:
        name=lambda w: EPI_MARKS[w.mark]["name"],
        prefix=lambda w: w.mark,
        table=lambda w: epi_contrasts(w.mark)[w.contrast]["table"],
        direction=lambda w: epi_contrasts(w.mark)[w.contrast]["direction"],
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        effect_name=lambda w: EPI_MARKS[w.mark]["effect_name"],
        columns=lambda w: epi_config(w.mark)["columns"],
        flip=lambda w: epi_config(w.mark)["flip"],
        fdr=lambda w: epi_config(w.mark)["fdr"],
        min_effect=lambda w: epi_config(w.mark)["min_effect"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        EPI_PREPARE_LOG,
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/epi_prepare.py"


# Copy of the regions file of the contrast, unchanged, so that the results
# keep the input they were made from (BigWig files are not copied)
# -----------------------------------------------------
rule epi_inputs:
    input:
        lambda w: epi_contrasts(w.mark)[w.contrast]["regions"],
    output:
        directory(EPI_INPUTS),
    params:
        files=lambda w: {"regions": epi_contrasts(w.mark)[w.contrast]["regions"]},
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        EPI_INPUTS_LOG,
    script:
        "../scripts/copy_inputs.py"


# The region each tested consensus peak overlaps most
# -----------------------------------------------------
rule epi_annotate_peaks:
    input:
        atac="results/deseq2/{contrast}/differential_peaks.tsv",
        regions=EPI_REGIONS,
    output:
        peaks="results/integration/{mark}/{contrast}/peaks.tsv",
    params:
        prefix=lambda w: w.mark,
        effect_name=lambda w: EPI_MARKS[w.mark]["effect_name"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/{mark}/annotate_peaks_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/epi_annotate_peaks.py"


use rule epi_annotate_peaks as te_epi_annotate_peaks with:
    input:
        atac="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        regions=EPI_REGIONS,
    output:
        peaks="results/te/integration/{mark}/{contrast}/peaks/peaks.tsv",
    params:
        prefix=lambda w: w.mark,
        effect_name=lambda w: EPI_MARKS[w.mark]["effect_name"],
    log:
        "logs/te/integration/{mark}/annotate_peaks_{contrast}.log",


# Differential peaks in the regions, and the accessibility changes by change
# of the mark
# -----------------------------------------------------
rule epi_peaks:
    input:
        peaks="results/integration/{mark}/{contrast}/peaks.tsv",
    output:
        enrichment="results/integration/{mark}/{contrast}/enrichment.tsv",
        enrichment_plot="results/integration/{mark}/{contrast}/enrichment.pdf",
        distribution="results/integration/{mark}/{contrast}/distribution.pdf",
        distribution_tests="results/integration/{mark}/{contrast}/distribution_tests.tsv",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        prefix=lambda w: w.mark,
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        changes=lambda w: epi_changes(w.mark, w.contrast),
        labels=lambda w: epi_labels(w.mark, w.contrast),
        peak_set="ENCODE consensus peaks (uniquely mapped reads)",
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/{mark}/peaks_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/epi_peaks.R"


use rule epi_peaks as te_epi_peaks with:
    input:
        peaks="results/te/integration/{mark}/{contrast}/peaks/peaks.tsv",
    output:
        enrichment="results/te/integration/{mark}/{contrast}/peaks/enrichment.tsv",
        enrichment_plot="results/te/integration/{mark}/{contrast}/peaks/enrichment.pdf",
        distribution="results/te/integration/{mark}/{contrast}/peaks/distribution.pdf",
        distribution_tests="results/te/integration/{mark}/{contrast}/peaks/distribution_tests.tsv",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        prefix=lambda w: w.mark,
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        changes=lambda w: epi_changes(w.mark, w.contrast),
        labels=lambda w: epi_labels(w.mark, w.contrast),
        peak_set="TE module peaks (alignment keeping multimapping reads)",
    log:
        "logs/te/integration/{mark}/peaks_{contrast}.log",


# Results tables only: ATAC-seq against the effect sizes of the mark for the
# peaks in tested regions
# -----------------------------------------------------
rule epi_scatter:
    input:
        peaks="results/integration/{mark}/{contrast}/peaks.tsv",
    output:
        correlation="results/integration/{mark}/{contrast}/correlation.tsv",
        scatter="results/integration/{mark}/{contrast}/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        prefix=lambda w: w.mark,
        effect_name=lambda w: EPI_MARKS[w.mark]["effect_name"],
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        labels=lambda w: epi_labels(w.mark, w.contrast),
        peak_set="ENCODE consensus peaks",
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=4000,
    log:
        "logs/integration/{mark}/scatter_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/epi_scatter.R"


use rule epi_scatter as te_epi_scatter with:
    input:
        peaks="results/te/integration/{mark}/{contrast}/peaks/peaks.tsv",
    output:
        correlation="results/te/integration/{mark}/{contrast}/peaks/correlation.tsv",
        scatter="results/te/integration/{mark}/{contrast}/peaks/scatter.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        prefix=lambda w: w.mark,
        effect_name=lambda w: EPI_MARKS[w.mark]["effect_name"],
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        labels=lambda w: epi_labels(w.mark, w.contrast),
        peak_set="TE module peaks",
    log:
        "logs/te/integration/{mark}/scatter_{contrast}.log",


# Signal of the mark ({mark}: contrasts: bigwigs) of the control and the
# treatment around the centre of the more, less and unchanged accessible peaks
# -----------------------------------------------------
rule epi_signal_profile:
    input:
        table="results/deseq2/{contrast}/differential_peaks.tsv",
        # Control first, as the reference (grey)
        bigwigs=lambda w: [
            epi_contrasts(w.mark)[w.contrast]["bigwigs"]["control"],
            epi_contrasts(w.mark)[w.contrast]["bigwigs"]["treatment"],
        ],
    output:
        matrix="results/integration/{mark}/{contrast}/signal_profile.matrix.gz",
        profile="results/integration/{mark}/{contrast}/signal_profile.pdf",
        heatmap="results/integration/{mark}/{contrast}/signal_heatmap.pdf",
    params:
        group_column="change",
        groups={
            "more accessible": "more accessible",
            "less accessible": "less accessible",
            "unchanged": "unchanged",
        },
        mode="reference-point",
        flank=lambda w: epi_config(w.mark)["profile_flank"],
        labels=lambda w: epi_conditions(w.contrast),
        colours=["#999999", "#dd3b3b"],
        y_label=lambda w: EPI_MARKS[w.mark]["labels"]["signal"],
        missing_as_zero=lambda w: EPI_MARKS[w.mark]["missing_as_zero"],
        region_name="peaks",
        extra_computematrix=extra("computematrix"),
        extra_plotprofile=extra("plotprofile"),
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=8000,
    log:
        "logs/integration/{mark}/signal_profile_{contrast}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/epi_profile.py"


use rule epi_signal_profile as te_epi_signal_profile with:
    input:
        table="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        bigwigs=lambda w: [
            epi_contrasts(w.mark)[w.contrast]["bigwigs"]["control"],
            epi_contrasts(w.mark)[w.contrast]["bigwigs"]["treatment"],
        ],
    output:
        matrix="results/te/integration/{mark}/{contrast}/peaks/signal_profile.matrix.gz",
        profile="results/te/integration/{mark}/{contrast}/peaks/signal_profile.pdf",
        heatmap="results/te/integration/{mark}/{contrast}/peaks/signal_heatmap.pdf",
    log:
        "logs/te/integration/{mark}/signal_profile_{contrast}.log",


# Tn5 insertions (signal tracks of the TE module, mean of the replicates) of
# the control and the treatment over the regions where the mark decreases,
# increases or does not change
# -----------------------------------------------------
use rule epi_signal_profile as te_epi_atac_profile with:
    input:
        table=EPI_REGIONS,
        bigwigs=lambda w: expand(
            f"results/te/bigwig/mean/{{condition}}.{TE_BIGWIG}.bigwig",
            condition=epi_conditions(w.contrast),
        ),
    output:
        matrix="results/te/integration/{mark}/{contrast}/atac_profile.matrix.gz",
        profile="results/te/integration/{mark}/{contrast}/atac_profile.pdf",
        heatmap="results/te/integration/{mark}/{contrast}/atac_heatmap.pdf",
    params:
        group_column="change",
        groups=lambda w: {
            EPI_MARKS[w.mark]["down"]: EPI_MARKS[w.mark]["labels"]["down_short"],
            EPI_MARKS[w.mark]["up"]: EPI_MARKS[w.mark]["labels"]["up_short"],
            "unchanged": EPI_MARKS[w.mark]["labels"]["unchanged"],
        },
        mode="scale-regions",
        flank=lambda w: epi_config(w.mark)["profile_flank"],
        body_length=lambda w: epi_config(w.mark)["profile_body_length"],
        labels=lambda w: epi_conditions(w.contrast),
        colours=["#999999", "#dd3b3b"],
        y_label=TE_BIGWIG_LABEL,
        missing_as_zero=True,
        region_name=lambda w: EPI_MARKS[w.mark]["labels"]["regions"],
        extra_computematrix=extra("computematrix"),
        extra_plotprofile=extra("plotprofile"),
    log:
        "logs/te/integration/{mark}/atac_profile_{contrast}.log",


# Repeat copies in the regions where the mark decreases and increases
# -----------------------------------------------------
rule te_epi_copies:
    input:
        saf="resources/te/repeats.saf",
        regions=EPI_REGIONS,
        # Only with te: profile: subfamilies (otherwise no full-length copies)
        full_length=lambda w: (
            ["resources/te/full_length_copies.bed"]
            if te_profile_config()["subfamilies"]
            else []
        ),
        peaks="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        copies="results/te/integration/{mark}/{contrast}/subfamilies/copies.tsv",
        full_length="results/te/integration/{mark}/{contrast}/subfamilies/full_length_copies.tsv",
    params:
        prefix=lambda w: w.mark,
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
    threads: 1
    resources:
        runtime=runtime(30),
        mem_mb=8000,
        cpus_per_task=cpus_for_memory(8000),
    log:
        "logs/te/integration/{mark}/copies_{contrast}.log",
    conda:
        "../envs/gat.yaml"
    script:
        "../scripts/epi_te_copies.py"


# Repeat subfamilies in the regions against their differential accessibility,
# and the full-length copies of the profile subfamilies
# -----------------------------------------------------
rule te_epi_subfamilies:
    input:
        copies="results/te/integration/{mark}/{contrast}/subfamilies/copies.tsv",
        full_length="results/te/integration/{mark}/{contrast}/subfamilies/full_length_copies.tsv",
        atac="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
    output:
        subfamilies="results/te/integration/{mark}/{contrast}/subfamilies/subfamilies.tsv",
        correlation="results/te/integration/{mark}/{contrast}/subfamilies/correlation.tsv",
        scatter="results/te/integration/{mark}/{contrast}/subfamilies/scatter.pdf",
        full_length_tests="results/te/integration/{mark}/{contrast}/subfamilies/full_length.tsv",
        full_length="results/te/integration/{mark}/{contrast}/subfamilies/full_length.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        prefix=lambda w: w.mark,
        down=lambda w: EPI_MARKS[w.mark]["down"],
        up=lambda w: EPI_MARKS[w.mark]["up"],
        changes=lambda w: epi_changes(w.mark, w.contrast),
        labels=lambda w: epi_labels(w.mark, w.contrast),
        min_copies=lambda w: epi_config(w.mark)["min_copies"],
        min_consensus_fraction=te_profile_config()["min_consensus_fraction"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=3000,
    log:
        "logs/te/integration/{mark}/subfamilies_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/epi_subfamilies.R"
