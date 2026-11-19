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

## Subsampling

`subsample_reads` randomly subsamples the Tn5 shifted tagAlign files of each replicate before
pseudoreplication and peak calling, as in the ENCODE pipeline: a number of reads (twice the
number of read pairs), or `min` for the number of reads of the smallest library after
duplicate removal (`results/tagalign/library_sizes.tsv`). Libraries with fewer reads keep all
of them. The tagAlign files of the TE module (`te: enable`) are subsampled in the same way, with
`min` from the TE alignments (`results/te/tagalign/library_sizes.tsv`).

Only peak calling and what is made from the tagAlign files use the subsampled reads: the peaks
(numbers, reproducibility, consensus peaks), FRiP, and the MACS2 fold-enrichment and p-value
signal tracks (`results/bigwig/`). The insertions in the consensus peaks and in the repeats are
counted from all reads, and the TE signal tracks are made from all reads: DESeq2 corrects for
sequencing depth with its size factors, and removing reads from the counts would only reduce
its power. The differential accessibility of the peaks can still change, as the consensus peaks
that are tested come from the subsampled peak calls; the test of the repeat subfamilies does
not change.

**When to use it:** MACS2 finds more peaks in deeper libraries, as weak peaks become
significant with more reads. If the libraries differ a lot in depth after duplicate removal
(e.g. more than 2-fold), the peaks of the deeper samples or conditions dominate the
reproducible and consensus peak sets, and comparisons of peak numbers or of overlaps between
conditions (e.g. peaks unique to one condition) reflect depth rather than biology. `min` gives
every replicate the same depth. With similar depths, or when only the differential tests matter,
leave it at 0, as subsampling discards reads. The depths after duplicate removal are in the
QC summaries (`nodup_non_mito_reads` in `results/qc/encode_qc_summary.tsv`, `kept_pairs` in
`results/te/qc/te_qc_summary.tsv`); compare those, not the raw read numbers, since duplicate
rates can differ a lot between libraries.

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
