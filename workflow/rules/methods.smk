from snakemake import __version__ as SNAKEMAKE_VERSION


# Package names of the environment files (read while parsing)
ENV_PACKAGES = {env: env_packages(env) for env in methods_envs()}


# Versions (name, version) of the packages listed in an environment file of
# the workflow, as installed in its conda environment (conda-meta records:
# name-version-build.json), without their dependencies. The Python and R
# interpreters are kept even when the file does not list them, as the
# workflow's scripts run on them.
# -----------------------------------------------------
rule software_versions:
    output:
        # Combined into methods/software_versions.tsv by the methods rule
        temp("results/software_versions/{env}.tsv"),
    wildcard_constraints:
        env="|".join(methods_envs()),
    params:
        packages=lambda w: " ".join(ENV_PACKAGES[w.env]),
        # The installed packages with these names
        awk=(
            'BEGIN { n = split(packages, a, " "); for (i = 1; i <= n; i++) '
            "keep[a[i]] = 1 } "
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
        "awk -F '\\t' -v packages='{params.packages}' '{params.awk}' | "
        "sort > {output}) "
        "2> {log}"


# Methods section (Markdown) of the parts of the workflow that were run, with
# the configuration values, software versions and references. Runs after all
# other targets, so that it describes the finished analysis
# -----------------------------------------------------
rule methods:
    input:
        versions=expand("results/software_versions/{env}.tsv", env=methods_envs()),
        targets=lambda w: ENCODE_TARGETS + TE_TARGETS,
    output:
        methods="methods/methods.md",
        # Package versions of all conda environments (environment, package,
        # version)
        versions="methods/software_versions.tsv",
    params:
        snakemake_version=SNAKEMAKE_VERSION,
        workflow_version=VERSION,
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
