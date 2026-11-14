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
