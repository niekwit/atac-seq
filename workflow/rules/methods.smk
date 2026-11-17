from snakemake import __version__ as SNAKEMAKE_VERSION


# Versions (name, version) of the packages listed in an environment file of
# the workflow, as installed in its conda environment (conda-meta records:
# name-version-build.json), without their dependencies. The Python and R
# interpreters are kept even when the file does not list them, as the
# workflow's scripts run on them.
# -----------------------------------------------------
rule software_versions:
    output:
        "methods/software_versions/{env}.tsv",
    wildcard_constraints:
        env="|".join(methods_envs()),
    params:
        env_file=lambda w: os.path.join(workflow.basedir, "envs", f"{w.env}.yaml"),
        # Package names from the dependency lines of the environment file
        # (bioconda::macs2=2.2.9.1=py310h7c593f9_1 -> macs2), then the
        # installed packages with these names
        awk=(
            "NR == FNR { if ($0 ~ /^[[:space:]]*-[[:space:]]/) { "
            's = $0; sub(/^[[:space:]]*-[[:space:]]*/, "", s); '
            'sub(/^[^:]*::/, "", s); sub(/[=<>! #].*/, "", s); keep[s] = 1 } '
            "next } "
            '($1 in keep) || $1 == "python" || $1 == "r-base"'
        ),
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
        "sed -nE 's/^(.+)-([^-]+)-([^-]+)\\.json$/\\1\\t\\2/p' | "
        "awk -F '\\t' '{params.awk}' {params.env_file} - | sort > {output}) "
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
        # From the config, so the data set files are not checked without TE
        gat_samples=config.get("te", {})
        .get("correlation", {})
        .get("gat_samples", 10000),
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
