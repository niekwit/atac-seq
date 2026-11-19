# Integration of all data sets of a contrast (optional, multiomics: enable,
# with the TE module): ATAC-seq, RNA-seq (rnaseq: tables), H3K9me3 and DNA
# methylation (h3k9me3: and methylation: contrasts, with BigWig files), for
# the contrasts that have all of them. Per repeat subfamily (heatmap of the
# four changes), per TE peak (UpSet plot and enrichment by combination of the
# marks) and per full-length copy of the te: profile: subfamilies entries
# (heatmap and profile of the signal of all data sets). Output in
# results/te/integration/multiomics/{contrast}/. Not part of the ENCODE
# pipeline.


# H3K9me3 and DNA methylation of the control and the treatment over all
# copies of each repeat subfamily
# -----------------------------------------------------
rule te_multiomics_signal:
    input:
        saf="resources/te/repeats.saf",
        bigwigs=lambda w: [
            epi_contrasts(mark)[w.contrast]["bigwigs"][condition]
            for mark in ("h3k9me3", "methylation")
            for condition in ("control", "treatment")
        ],
    output:
        "results/te/integration/multiomics/{contrast}/subfamilies/signal.tsv",
    params:
        names=[
            "h3k9me3_control",
            "h3k9me3_treatment",
            "methylation_control",
            "methylation_treatment",
        ],
        missing_as_zero=[
            EPI_MARKS["h3k9me3"]["missing_as_zero"],
            EPI_MARKS["h3k9me3"]["missing_as_zero"],
            EPI_MARKS["methylation"]["missing_as_zero"],
            EPI_MARKS["methylation"]["missing_as_zero"],
        ],
    threads: 1
    resources:
        runtime=runtime(120),
        mem_mb=16000,
        cpus_per_task=cpus_for_memory(16000),
    log:
        "logs/te/integration/multiomics/signal_{contrast}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/multiomics_signal.py"


# Heatmap of the changes in accessibility, expression, H3K9me3 and DNA
# methylation of the repeat subfamilies
# -----------------------------------------------------
rule te_multiomics_subfamilies:
    input:
        atac="results/te/deseq2/{contrast}/differential_subfamilies.tsv",
        rnaseq=RNASEQ_TABLE,
        signal="results/te/integration/multiomics/{contrast}/subfamilies/signal.tsv",
        h3k9me3="results/te/integration/h3k9me3/{contrast}/subfamilies/copies.tsv",
        methylation="results/te/integration/methylation/{contrast}/subfamilies/copies.tsv",
        annotation="resources/te/subfamilies.tsv",
    output:
        table="results/te/integration/multiomics/{contrast}/subfamilies/subfamilies.tsv",
        heatmap="results/te/integration/multiomics/{contrast}/subfamilies/heatmap.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        marks=lambda w: multiomics_marks(w.contrast),
        min_copies=multiomics_config()["min_copies"],
        max_subfamilies=multiomics_config()["max_subfamilies"],
        # Subfamilies of these entries are always in the heatmap
        profile_entries=te_profile_config()["subfamilies"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=3000,
    log:
        "logs/te/integration/multiomics/subfamilies_{contrast}.log",
    conda:
        "../envs/multiomics.yaml"
    script:
        "../scripts/multiomics_subfamilies.R"


# TE peaks by change in accessibility, H3K9me3, DNA methylation and
# expression of the nearest gene: UpSet plot, and enrichment by combination
# of the marks
# -----------------------------------------------------
rule te_multiomics_peaks:
    input:
        atac="results/te/peaks/deseq2/{contrast}/differential_peaks.tsv",
        h3k9me3="results/te/integration/h3k9me3/{contrast}/peaks/peaks.tsv",
        methylation="results/te/integration/methylation/{contrast}/peaks/peaks.tsv",
        rnaseq=RNASEQ_TABLE,
    output:
        table="results/te/integration/multiomics/{contrast}/peaks/peaks.tsv",
        upset="results/te/integration/multiomics/{contrast}/peaks/upset.pdf",
        enrichment="results/te/integration/multiomics/{contrast}/peaks/enrichment.tsv",
        enrichment_plot="results/te/integration/multiomics/{contrast}/peaks/enrichment.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        marks=lambda w: multiomics_marks(w.contrast),
        max_distance=max(rnaseq_config()["distance_bins"]),
        max_combinations=multiomics_config()["max_combinations"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=3000,
    log:
        "logs/te/integration/multiomics/peaks_{contrast}.log",
    conda:
        "../envs/multiomics.yaml"
    script:
        "../scripts/multiomics_peaks.R"


# ATAC-seq, H3K9me3 and CpG methylation of the control and the treatment over
# the full-length copies of the te: profile: subfamilies entries: matrix,
# heatmap (copies sorted by the change in accessibility) and mean profiles
# -----------------------------------------------------
rule te_multiomics_copies:
    input:
        bed="resources/te/full_length_copies.bed",
        # Control first, then the treatment, of each data set
        bigwigs=lambda w: [
            *expand(
                f"results/te/bigwig/mean/{{condition}}.{TE_BIGWIG}.bigwig",
                condition=epi_conditions(w.contrast),
            ),
            *[
                epi_contrasts(mark)[w.contrast]["bigwigs"][condition]
                for mark in ("h3k9me3", "methylation")
                for condition in ("control", "treatment")
            ],
        ],
        h3k9me3="results/te/integration/h3k9me3/{contrast}/subfamilies/full_length_copies.tsv",
        methylation="results/te/integration/methylation/{contrast}/subfamilies/full_length_copies.tsv",
    output:
        matrix="results/te/integration/multiomics/{contrast}/copies/matrix.gz",
        copies="results/te/integration/multiomics/{contrast}/copies/copies.tsv",
        heatmap="results/te/integration/multiomics/{contrast}/copies/heatmap.pdf",
        profile="results/te/integration/multiomics/{contrast}/copies/profile.tsv",
    params:
        datasets=["ATAC-seq"] * 2 + ["H3K9me3"] * 2 + ["Methylation"] * 2,
        columns=["atac"] * 2 + ["h3k9me3"] * 2 + ["methylation"] * 2,
        roles=["control", "treatment"] * 3,
        conditions=lambda w: epi_conditions(w.contrast) * 3,
        missing_as_zero=[True, True]
        + [EPI_MARKS["h3k9me3"]["missing_as_zero"]] * 2
        + [EPI_MARKS["methylation"]["missing_as_zero"]] * 2,
        colour_maps=["Reds"] * 2 + ["Blues"] * 2 + ["Greens"] * 2,
        flank=te_profile_config()["flank"],
        body_length=te_profile_config()["body_length"],
        bin_size=computematrix_bin_size(),
    threads: 4
    resources:
        runtime=runtime(60),
        mem_mb=8000,
    log:
        "logs/te/integration/multiomics/copies_{contrast}.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/multiomics_copies.py"


rule te_multiomics_profile:
    input:
        "results/te/integration/multiomics/{contrast}/copies/profile.tsv",
    output:
        "results/te/integration/multiomics/{contrast}/copies/profile.pdf",
    params:
        treatment=lambda w: deseq2_contrasts()[w.contrast][0],
        control=lambda w: deseq2_contrasts()[w.contrast][1],
        flank=te_profile_config()["flank"],
        body_length=te_profile_config()["body_length"],
    threads: 1
    resources:
        runtime=runtime(15),
        mem_mb=2000,
    log:
        "logs/te/integration/multiomics/profile_{contrast}.log",
    conda:
        "../envs/integration.yaml"
    script:
        "../scripts/multiomics_profile.R"
