# Snakemake workflow: `atac-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/atac-seq/actions/workflows/main.yaml/badge.svg)](https://github.com/niekwit/atac-seq/actions/workflows/main.yaml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/niekwit/atac-seq)

A Snakemake workflow for `ATAC-seq` data analysis that follows the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
See [config/README.md](config/README.md) for the workflow steps, parameters and output.

- [Snakemake workflow: `atac-seq`](#snakemake-workflow-atac-seq)
  - [Usage](#usage)
    - [Requirements](#requirements)
    - [Deploy the workflow](#deploy-the-workflow)
    - [Configure the analysis](#configure-the-analysis)
    - [Run the workflow](#run-the-workflow)
  - [Authors](#authors)
  - [References](#references)

## Usage

### Requirements

- Linux with [conda](https://docs.conda.io/en/latest/) (or mamba) installed
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) ≥8.25.5:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=8.25.5 snakedeploy
conda activate snakemake
```

All other software is installed automatically by Snakemake in conda environments.

### Deploy the workflow

Create an analysis directory and deploy the workflow into it with
[snakedeploy](https://snakedeploy.readthedocs.io):

```bash
mkdir -p path/to/analysis && cd path/to/analysis
snakedeploy deploy-workflow https://github.com/niekwit/atac-seq . --branch main
```

This creates `workflow/Snakefile`, which loads the workflow from GitHub, and copies the
default configuration to `config/`. Alternatively, clone the repository and run the workflow
from the analysis directory with `--snakefile path/to/atac-seq/workflow/Snakefile`.

### Configure the analysis

The analysis directory should look like this:

```
path/to/analysis
├── config
│   ├── config.yaml
│   └── samples.csv
└── reads
    ├── WT_1_R1_001.fastq.gz
    ├── WT_1_R2_001.fastq.gz
    ├── ...
```

1. Put the paired-end reads in `reads/` as `{sample}_R1_001.fastq.gz` and `{sample}_R2_001.fastq.gz`.
2. List the samples and their condition in `config/samples.csv`. Samples of the same condition are
   analysed as biological replicates.
3. Set the genome (`genome: ensembl` and `release`) and the effective genome size
   (`macs2: gsize`) in `config/config.yaml`. The other defaults follow the ENCODE pipeline.

See [config/README.md](config/README.md) for the sample sheet rules, all parameters and the output files.

### Run the workflow

Check what will be run with a dry run:

```bash
snakemake -n
```

Then run the workflow, letting Snakemake create the conda environments:

```bash
snakemake --sdm conda --cores 16
```

Genome files and the blacklist are downloaded on the first run. The main results are the
reproducible peak sets in `results/idr/{condition}/` and the QC summary in
`results/qc/encode_qc_summary.tsv`.

To run the workflow on a compute cluster, use an
[executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/) with a profile,
for example `snakemake --sdm conda --executor slurm --jobs 50`.

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)
  - [Google Scholar](https://scholar.google.co.uk/citations?user=USu7NNcAAAAJ&hl=en)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
