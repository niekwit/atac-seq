# Snakemake workflow: `atac-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/atac-seq/actions/workflows/main.yaml/badge.svg)](https://github.com/niekwit/atac-seq/actions/workflows/main.yaml)
[![Documentation](https://readthedocs.org/projects/atac-seq/badge/?version=latest)](https://atac-seq.readthedocs.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/niekwit/atac-seq)

A Snakemake workflow for paired-end `ATAC-seq` data that follows the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) and grades the
results against the [ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/).

**Documentation: <https://atac-seq.readthedocs.io>**

## Features

- ENCODE peaks (IDR and naive overlap), signal tracks and QC, validated against ENCODE's own
  processing of human and mouse data
- Peaks reproducible in several replicates of a condition
- Comparison of all samples with DESeq2: PCA, sample correlation and differential
  accessibility between conditions
- Optional module for transposable elements, including young repeats whose reads align to many
  places: counts and differential accessibility per subfamily, signal tracks, profiles, TE
  peaks and comparison with other data sets (e.g. DMRs)

## Quick start

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=9.27.1 snakedeploy
conda activate snakemake
mkdir -p path/to/analysis && cd path/to/analysis
snakedeploy deploy-workflow https://github.com/niekwit/atac-seq . --branch main
```

Put the paired-end reads in `reads/` as `{sample}_R1_001.fastq.gz` and
`{sample}_R2_001.fastq.gz`, list the samples in `config/samples.csv` and set the genome in
`config/config.yaml` (see [config/README.md](config/README.md)). Then run:

```bash
snakemake --sdm conda --cores 16 --resources mem_mb=100000
```

See the [documentation](https://atac-seq.readthedocs.io) for running on a cluster, the output
files and the methods.

## Building the documentation

```bash
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)
  - [Google Scholar](https://scholar.google.co.uk/citations?user=USu7NNcAAAAJ&hl=en)

## Citation

If you use this workflow, please cite Snakemake and the tools it uses (see the
[references](https://atac-seq.readthedocs.io/en/latest/references.html)):

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
