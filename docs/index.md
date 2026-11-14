# atac-seq

A [Snakemake](https://snakemake.github.io) workflow for paired-end ATAC-seq data that follows
the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) (v2.2.3) and
grades the results against the
[ENCODE ATAC-seq data standards](https://www.encodeproject.org/atac-seq/). Samples of the same
condition are analysed as biological replicates of one ENCODE experiment.

On top of the ENCODE steps, the workflow can:

- compare all samples with DESeq2: PCA, sample correlation and differential accessibility
  between conditions ([Comparing samples](comparing_samples.md));
- give the peaks reproducible in several replicates of a condition
  ([QC and reproducibility](qc.md#peaks-reproducible-in-several-replicates));
- measure the accessibility of transposable elements, including young repeats whose reads
  align to many places ([Transposable elements](transposable_elements.md)).

The peaks, signal tracks and QC metrics match ENCODE's own processing of the same data
([Validation](validation.md)).

## Quick start

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=9.27.1 snakedeploy
conda activate snakemake
mkdir -p path/to/analysis && cd path/to/analysis
snakedeploy deploy-workflow https://github.com/niekwit/atac-seq . --branch main
# put the reads in reads/, edit config/samples.csv and config/config.yaml
snakemake --sdm conda --cores 16 --resources mem_mb=100000
```

See [Installation](installation.md), [Configuration](configuration.md) and
[Running the workflow](running.md) for the details.

```{toctree}
:caption: User guide
:maxdepth: 2
:hidden:

installation
configuration
running
outputs
```

```{toctree}
:caption: Methods
:maxdepth: 2
:hidden:

workflow
qc
comparing_samples
transposable_elements
```

```{toctree}
:caption: Validation
:maxdepth: 2
:hidden:

validation
te_alignment_comparison/README
```

```{toctree}
:caption: About
:maxdepth: 1
:hidden:

references
GitHub repository <https://github.com/niekwit/atac-seq>
```
