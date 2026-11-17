from snakemake import __version__ as SNAKEMAKE_VERSION


# Packages installed in a conda environment of the workflow (name, version),
# from the conda-meta records (name-version-build.json) of the environment
# -----------------------------------------------------
rule software_versions:
    output:
        "methods/software_versions/{env}.tsv",
    wildcard_constraints:
        env="|".join(methods_envs()),
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        "logs/methods/software_versions_{env}.log",
    conda:
        "../envs/{env}.yaml"
    shell:
        '(ls "$CONDA_PREFIX/conda-meta" | '
        "sed -nE 's/^(.+)-([^-]+)-([^-]+)\\.json$/\\1\\t\\2/p' | sort > {output}) "
        "2> {log}"


# Methods section (Markdown) of the parts of the workflow that were run, with
# the configuration values, software versions and references. Runs after all
# other targets, so that it describes the finished analysis
# -----------------------------------------------------
rule methods:
    input:
        versions=expand("methods/software_versions/{env}.tsv", env=methods_envs()),
        targets=lambda w: ENCODE_TARGETS + TE_TARGETS,
    output:
        "methods/methods.md",
    params:
        snakemake_version=SNAKEMAKE_VERSION,
        workflow_version=workflow_version(),
        replicates={c: len(replicates(c)) for c in CONDITIONS},
        reference=reference_files(),
        te_only=te_only(),
        te=te_enabled(),
        deseq2=deseq2_enabled(),
        peak_method=deseq2_peak_method(),
        contrasts=lambda w: deseq2_contrasts() if deseq2_enabled() else {},
        profile=te_profile_config(),
        heatmap=te_heatmap_config(),
        # Functions, so the data sets are checked only when the TE module runs
        datasets=lambda w: (
            list(te_correlation_config()["datasets"])
            if te_enabled() and deseq2_enabled()
            else []
        ),
        gat_samples=lambda w: te_correlation_config()["gat_samples"],
    threads: 1
    resources:
        runtime=runtime(10),
        mem_mb=1000,
    log:
        "logs/methods/methods.log",
    conda:
        "../envs/reports.yaml"
    script:
        "../scripts/methods.py"
