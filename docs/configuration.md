# Configuration

## Analysis directory

After [deploying the workflow](installation.md), the analysis directory should look like this:

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

1. Put the paired-end reads in `reads/`.
2. List the samples and their condition in `config/samples.csv`. Samples of the same condition
   are analysed as biological replicates.
3. Set the genome and the effective genome size in `config/config.yaml`.

```{include} ../config/README.md
:start-after: <!-- docs-include-start -->
```

## Extra tool arguments

Tools that the workflow runs directly as a command take extra command-line arguments from
`extra:` in `config/config.yaml`. They are added after the workflow's own options, so for most
tools they override them, and the [methods](outputs.md#methods) list the ones that were used. For
example, to give the MultiQC reports a title and to plot the standard error in the TE profiles:

```yaml
extra:
  multiqc: "--title 'My project'"
  plotprofile: "--plotType se"
```

Tools run from the workflow's own scripts (e.g. MACS2, IDR and `plotFingerprint`) are not
included. cutadapt and bowtie2 take extra arguments from `cutadapt: extra` and `bowtie2: extra`.
