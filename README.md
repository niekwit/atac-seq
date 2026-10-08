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
  - [QC summary](#qc-summary)
  - [Validation against ENCODE](#validation-against-encode)
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

Each rule has a runtime limit (Snakemake's `runtime` resource, in minutes), which cluster
executors pass on as the job's time limit. The limits are about 3-5 times the runtimes for two
replicates of ~50 million read pairs (human), and are multiplied by the attempt number, so jobs
that run out of time are resubmitted with more time when the workflow is run with `--retries`
(e.g. `--retries 2`). For much deeper data, or to change a limit, override it in a profile or on
the command line, e.g. `--set-resources bowtie2_align:runtime=1440`.

## QC summary

`results/qc/encode_qc_summary.tsv` has one row per sample with the QC metrics of the ENCODE
pipeline. Most are graded as `ideal`, `acceptable` or `concerning` against the
[ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/#standards):

| Column                 | Meaning                                                                                                   | Ideal        | Acceptable |
| ---------------------- | --------------------------------------------------------------------------------------------------------- | ------------ | ---------- |
| `alignment_rate`       | Fraction of read pairs aligned by bowtie2                                                                 | > 0.95       | > 0.80     |
| `frac_mito_reads`      | Fraction of aligned reads on the mitochondrial genome (not graded; lower is better)                       |              |            |
| `nodup_non_mito_reads` | Reads left after filtering and removal of duplicates and mitochondrial reads; these are used for peak calling | ≥ 50 million |            |
| `NRF`                  | Non-redundant fraction: distinct fragments / all fragments. Low values mean many PCR duplicates           | > 0.9        | > 0.8      |
| `PBC1`                 | PCR bottlenecking coefficient 1: positions with exactly one fragment / positions with at least one        | > 0.9        | > 0.8      |
| `PCR_bottlenecking`    | Bottlenecking by PBC1 as in the ENCODE QC report: none (> 0.9), mild (0.8-0.9), moderate (0.5-0.8), severe (< 0.5) |  |  |
| `PBC2`                 | PCR bottlenecking coefficient 2: positions with one fragment / positions with two                         | > 3          | > 1        |
| `num_peaks`            | Peaks of the replicate (MACS2 p < 0.01, top 300,000, blacklist filtered)                                  |              |            |
| `FRiP`                 | Fraction of Tn5 cutting sites in the replicate's peaks                                                    | > 0.3        | > 0.2      |
| `TSS_enrichment`       | Signal at transcription start sites relative to the background 2 kb away, calculated as by ENCODE        | genome dependent¹ |       |
| `TSS_enrichment_ataqv` | TSS enrichment calculated by ataqv (not graded, see below)                                                |              |            |

¹ hg38: > 7 ideal, > 5 acceptable; hg19: > 10 ideal, > 6 acceptable; mm10/mm39: > 15 ideal,
> 10 acceptable.

`TSS_enrichment` follows the ENCODE pipeline (`encode_task_tss_enrich.py`): reads are centred
on their 5' end, their coverage is averaged over ±2 kb windows around the TSSs of protein
coding genes, and the maximum is divided by the coverage at the window edges. As in ENCODE,
reads near the window edges are partly excluded, which raises the score by about a third. The
ENCODE thresholds only apply to scores calculated this way: ataqv calculates TSS enrichment on
a different scale (about 4 where ENCODE reports about 29), so `TSS_enrichment_ataqv` is not
graded.

`results/qc/encode_reproducibility_summary.tsv` has one row per condition and reproducibility
method (IDR and overlap): the optimal and conservative peak sets, their number of peaks
(ENCODE: > 150,000 ideal, > 100,000 acceptable) and FRiP, and the rescue and self-consistency
ratios. A ratio above 2 means poor agreement between replicates (rescue ratio) or within a
replicate (self-consistency ratio); with one ratio above 2 reproducibility is `borderline`,
with both `fail`.

## Validation against ENCODE

The workflow was run with default settings (hg38, Ensembl release 115) on the raw reads of
ENCODE experiment [ENCSR422SUG](https://www.encodeproject.org/experiments/ENCSR422SUG/)
(MCF-7, two paired-end replicates), and its output was compared with the files ENCODE
processed with its ATAC-seq pipeline (v1.9.1).

Signal tracks (Pearson correlation of the mean signal in 1 kb bins on chr1-22 and X, and in
the union of both peak sets):

| Signal track                    | 1 kb bins | Peaks  |
| ------------------------------- | --------- | ------ |
| Pooled replicates, p-value      | 0.9999    | 0.9999 |
| Pooled replicates, fold change  | 0.9999    | 0.9996 |
| Replicate 1, p-value            | 0.9999    | 0.9999 |
| Replicate 2, p-value            | 0.9999    | 0.9999 |

Peaks:

| Peak set                          | Peaks (workflow / ENCODE) | Overlapping ENCODE peaks | Identical to ENCODE peaks¹ |
| --------------------------------- | ------------------------- | ------------------------ | -------------------------- |
| Overlap, optimal (ENCODE default) | 245,849 / 248,668         | 95%                      | 92%                        |
| IDR, optimal                      | 171,837 / 173,574         | 92%                      | 89%                        |
| IDR, conservative                 | 167,086 / 168,691         | 99%                      | 96%                        |

¹ Same coordinates and summit. The signal values and p-values of these peaks are identical
to ENCODE's (Pearson r = 1.000).

The remaining differences are expected from the different genome assembly (Ensembl primary
assembly instead of ENCODE's GRCh38 no-alt analysis set) and the newer pipeline and tool
versions. Both replicates passed the ENCODE standards for alignment rate, read depth, FRiP
and reproducibility.

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)
  - [Google Scholar](https://scholar.google.co.uk/citations?user=USu7NNcAAAAJ&hl=en)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
