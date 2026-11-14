# Installation

## Requirements

- Linux with [conda](https://docs.conda.io/en/latest/) (or mamba) installed
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) ≥8.25.5
  (tested with 8.25.5 and 9.27.1)

Install Snakemake and snakedeploy in a new environment, so that conda picks a Python version
that both support:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=9.27.1 snakedeploy
conda activate snakemake
```

On a cluster, also install the executor plugin of its scheduler, e.g.
`snakemake-executor-plugin-slurm` for SLURM.

All other software is installed automatically by Snakemake in conda environments.

## Deploy the workflow

Create an analysis directory and deploy the workflow into it with
[snakedeploy](https://snakedeploy.readthedocs.io):

```bash
mkdir -p path/to/analysis && cd path/to/analysis
snakedeploy deploy-workflow https://github.com/niekwit/atac-seq . --branch main
```

This creates `workflow/Snakefile`, which loads the workflow from GitHub, and copies the
default configuration to `config/`. To use a fixed version, deploy a release instead of `main`,
e.g. `--tag v0.3.0`.

Alternatively, clone the repository and run the workflow from the analysis directory with
`--snakefile path/to/atac-seq/workflow/Snakefile`.

Next: [configure the analysis](configuration.md).
