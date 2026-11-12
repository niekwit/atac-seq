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
  - [Comparing samples](#comparing-samples)
  - [Transposable elements](#transposable-elements)
  - [Validation against ENCODE](#validation-against-encode)
  - [Authors](#authors)
  - [References](#references)

## Usage

### Requirements

- Linux with [conda](https://docs.conda.io/en/latest/) (or mamba) installed
- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) ≥8.25.5
  (tested with 8.25.5 and 9.27.1). Install it in a new environment, so that conda picks a
  Python version that Snakemake and snakedeploy support:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake=9.27.1 snakedeploy
conda activate snakemake
```

On a cluster, also install the executor plugin of its scheduler, e.g.
`snakemake-executor-plugin-slurm` for SLURM.

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
   analysed as biological replicates. Mark the samples of the reference (control) condition with
   `yes` in the `reference` column.
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
snakemake --sdm conda --cores 16 --resources mem_mb=100000
```

Each rule declares how much memory it needs (`mem_mb`), so `--resources mem_mb=100000` makes
Snakemake run only as many jobs at the same time as fit in 100 GB of RAM. Set it to somewhat
less than the RAM of your computer. The largest jobs need about 10-16 GB (genome index,
alignment, signal tracks, annotation).

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
the command line, e.g. `--set-resources bowtie2_align:runtime=1440`. The memory limits
(`mem_mb`) are passed on as well, for clusters that schedule by memory. For clusters that
assign memory per CPU, each rule requests enough threads to cover its `mem_mb` at 3,370 MiB per
CPU; on clusters with less memory per CPU, raise the threads with `--set-threads`.

## QC summary

`results/qc/encode_qc_summary.tsv` has one row per sample with the QC metrics of the ENCODE
pipeline. Most are graded as `ideal`, `acceptable` or `concerning` against the
[ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/#standards):

| Column                 | Meaning                                                                                                   | Ideal        | Acceptable |
| ---------------------- | --------------------------------------------------------------------------------------------------------- | ------------ | ---------- |
| `alignment_rate`       | Fraction of read pairs aligned by bowtie2                                                                 | > 0.95       | > 0.80     |
| `frac_mito_reads`      | Fraction of aligned reads that align to the mitochondrial genome, as by ENCODE (not graded; lower is better) |              |            |
| `nodup_non_mito_reads` | Reads left after filtering and removal of duplicates and mitochondrial reads; these are used for peak calling | ≥ 50 million |            |
| `NRF`                  | Non-redundant fraction: distinct fragments / all fragments. Low values mean many PCR duplicates           | > 0.9        | > 0.8      |
| `PBC1`                 | PCR bottlenecking coefficient 1: positions with exactly one fragment / positions with at least one        | > 0.9        | > 0.8      |
| `PCR_bottlenecking`    | Bottlenecking by PBC1 as in the ENCODE QC report: none (> 0.9), mild (0.8-0.9), moderate (0.5-0.8), severe (< 0.5) |  |  |
| `PBC2`                 | PCR bottlenecking coefficient 2: positions with one fragment / positions with two                         | > 3          | > 1        |
| `num_peaks`            | Peaks of the replicate (MACS2 p < 0.01, top 300,000, blacklist filtered)                                  |              |            |
| `FRiP`                 | Fraction of Tn5 cutting sites in the replicate's peaks                                                    | > 0.3        | > 0.2      |
| `TSS_enrichment`       | Signal at transcription start sites relative to the background 2 kb away, calculated as by ENCODE        | genome dependent¹ |       |
| `TSS_enrichment_ataqv` | TSS enrichment calculated by ataqv (not graded, see below)                                                |              |            |
| `frac_fragments_in_NFR`, `NFR_mono_nuc_ratio` | Fraction of fragments < 150 bp (nucleosome free region, NFR), and NFR / mononucleosome (150-300 bp) fragments (ENCODE flags < 0.4 and < 2.5) | | |
| `NFR_peak`, `mono_nuc_peak`, `di_nuc_peak` | Peaks in the fragment length distribution at 20-90, 120-250 and 300-500 bp            |              |            |
| `nucleosome_pattern_standard` | ENCODE standard: a NFR and a mononucleosome peak are present                                       | both present |            |
| `synthetic_JS_distance` | Jensen-Shannon distance from uniform coverage (deepTools fingerprint); higher means more enriched signal |              |            |
| `frac_reads_in_DHS`, `_promoters`, `_enhancers`, `_blacklist` | Fraction of Tn5 cutting sites in ENCODE's DNase I hypersensitive sites, promoters, enhancers and exclusion list (DHS, promoters and enhancers for GRCh38 and mm10 only) | | |
| `AT_dropout`, `GC_dropout` | Picard GC bias: reads missing from AT- or GC-rich regions (%)                                         |              |            |

¹ hg38: > 7 ideal, > 5 acceptable; hg19: > 10 ideal, > 6 acceptable; mm10/mm39: > 15 ideal,
> 10 acceptable.

`TSS_enrichment` follows the ENCODE pipeline (`encode_task_tss_enrich.py`): reads are centred
on their 5' end, their coverage is averaged over ±2 kb windows around the TSSs of protein
coding genes, and the maximum is divided by the coverage at the window edges. The TSSs are
those of ENCODE's TSS files (GRCh38, mm10; lifted over for mm39), as the score depends on the
annotation: in recent Ensembl releases the 5' end of many genes lies upstream of the main TSS,
which lowers the score (by ~15% for GRCm39, Ensembl 116), whereas the TSSs of canonical
transcripts raise it (by ~50% for GRCh38) above the scale of the ENCODE thresholds. As in ENCODE,
reads near the window edges are partly excluded, which raises the score by about a third. The
ENCODE thresholds only apply to scores calculated this way: ataqv calculates TSS enrichment on
a different scale (about 5 where ENCODE reports about 29), so `TSS_enrichment_ataqv` is not
graded.

`results/qc/encode_reproducibility_summary.tsv` has one row per condition and reproducibility
method (IDR and overlap): the optimal and conservative peak sets, their number of peaks
(ENCODE: > 150,000 ideal, > 100,000 acceptable) and FRiP, and the rescue and self-consistency
ratios. A ratio above 2 means poor agreement between replicates (rescue ratio) or within a
replicate (self-consistency ratio); with one ratio above 2 reproducibility is `borderline`,
with both `fail`.

The optimal and conservative peak sets come from a single pair of replicates (ENCODE picks the
pair with the most peaks), so with three or more replicates a peak only has to be reproducible
in two of them. `results/{idr,overlap}/{condition}/replicated_peak.narrowPeak` (not part of the
ENCODE pipeline) has the peaks of the pooled sample that are reproducible in at least *k* of the
*n* replicates of a condition, with *k* = `peaks: min_replicate_fraction` × *n*, rounded up and at
least 2 (default 1: all replicates):

- naive overlap: the peak overlaps a peak of at least *k* replicates by at least 50% (of either
  peak), the criterion of the ENCODE naive overlap (`bedtools intersect -f 0.5 -F 0.5 -e`)
- IDR: IDR only compares pairs, so there must be *k* replicates of which every pair passed the
  IDR threshold for the peak; its score (column 5) is the lowest scaled IDR of these pairs

With two replicates, or with *k* = *n*, the peaks are in the peak sets of all pairs of
replicates; with two replicates this is the conservative peak set. Requiring all replicates
makes the set smaller with every replicate added, and a single weaker replicate limits it: the
log (`logs/replicated_peaks/`) lists the peaks of each replicate (overlap) or pair (IDR). With
different numbers of replicates per condition the sets are not equally strict, so compare
conditions with the DESeq2 analysis, not with these sets. As in the other peak sets, a region
can be listed more than once, with different summits (`--call-summits`).

## Comparing samples

The ENCODE pipeline analyses each condition on its own. To compare all samples, the workflow
also uses [DESeq2](https://doi.org/10.1186/s13059-014-0550-8) (`deseq2: enable: true`):

1. Consensus peaks: the optimal peak sets of all conditions (`deseq2: peak_method`, overlap or
   IDR) are merged with `bedtools merge`
2. Tn5 insertions (the 5' end of each read of `{sample}.nodup.bam`) in the consensus peaks are
   counted with `featureCounts`
3. Variance stabilising transformation (blind to the conditions) of the peaks with at least 10
   insertions, which is used for the PCA and the sample correlations
4. Differential accessibility between conditions, for each `[treatment, control]` pair in
   `deseq2: contrasts` or, if that is empty, each condition with at least 2 replicates against
   the reference condition (`reference` column of the sample sheet, see
   [config/README.md](config/README.md#input-data)): DESeq2 model `~ condition` fitted to all
   samples, Wald test, and log2 fold changes shrunken with
   [ashr](https://doi.org/10.1093/biostatistics/kxw041)

Output in `results/deseq2/`:

| file                                                  | content                                                     |
| ----------------------------------------------------- | ----------------------------------------------------------- |
| `consensus_peaks.bed`                                 | consensus peaks                                             |
| `consensus_peak_counts.txt`                           | Tn5 insertions per peak and sample                          |
| `vst_counts.tsv`                                      | variance stabilised counts                                  |
| `pca.png`, `pca.tsv`                                  | PCA of the most variable peaks (`deseq2: pca_ntop`)         |
| `sample_correlation.png`, `sample_correlation.tsv`    | Pearson correlation of the samples over all peaks           |
| `consensus_peaks_annotated.txt`                       | consensus peaks annotated with `ChIPseeker`                 |
| `dds.rds`                                             | fitted DESeq2 object                                        |
| `{treatment}_vs_{control}/differential_peaks.tsv`     | peak coordinates, annotation and test results, by p-value   |
| `{treatment}_vs_{control}/ma_plot.png`, `volcano_plot.png` | plots of the shrunken log2 fold changes                |

In `differential_peaks.tsv`, `change` is `more accessible` or `less accessible` in the treatment
if `padj` is below `deseq2: fdr`. The p-values come from the unshrunken `log2FoldChange` (against
`deseq2: lfc_threshold`); use `log2FoldChange_shrunken` for ranking and plots.

DESeq2 scales the samples with median-of-ratios size factors, calculated over the peaks. This
assumes that most peaks do not change, so a global change in accessibility, or a large
difference in FRiP between samples, is not measured reliably: check `FRiP` in the QC summary
and the MA plot, in which the cloud of peaks should be centred on 0.

## Transposable elements

Young transposable elements (TEs), such as L1Md_T/A/Gf LINE-1s and IAPs in mouse or L1HS and
L1PA2/3 in human, have many nearly identical copies. Their reads align to many places, so the
ENCODE filter (at most 4 alignments, one of which is kept at random) removes them, and they
hardly show up in the peaks and signal tracks. Older, diverged copies map uniquely and are
analysed like the rest of the genome.

The optional TE module (`te: enable: true` in `config/config.yaml`) measures the accessibility
of repeats per subfamily, family and class instead of per copy. It does not change the ENCODE
outputs.

1. RepeatMasker annotation from UCSC (`rmsk`), with Ensembl chromosome names
2. A second alignment of the trimmed reads that keeps multimappers (`bowtie2 --seed 42 -X2000
   --mm`): bowtie2's default mode reports one of the best alignments of each read pair, picked
   at random (reproducibly, with the seed) among equally good ones
3. Properly paired read pairs are kept; PCR duplicates and reads on `filter_chrs` (MT) are
   removed. Read pairs in which a mate has an equally good alternative alignment (bowtie2
   `XS` ≥ `AS`) that the pairing does not resolve (MAPQ < 2) are marked as multimappers.
4. Counting of Tn5 insertions (5' end of each read) in repeats with `featureCounts`, once with
   all read pairs (`-M`) and once with uniquely mapping read pairs only
5. Signal tracks of Tn5 insertions (5' end of each read, 10 bp bins, per million insertions;
   `deepTools bamCoverage`) per sample, and the mean of the replicates of each condition
   (`deepTools bigwigAverage`). Unlike the ENCODE signal tracks, they include multimappers.
6. With `deseq2: enable`: differential accessibility of subfamilies with `DESeq2`, for the same
   contrasts, `fdr` and `lfc_threshold` as the peaks (Wald test, ashr-shrunken log2 fold
   changes for plots). Library sizes come from the median-of-ratios over the subfamilies
   (`te: normalisation: deseq2`), which assumes that most subfamilies do not change. If a
   condition changes the accessibility of many repeat families at once, use
   `te: normalisation: insertions` (all kept insertions of each sample) instead.
7. With subfamilies in `te: profile: subfamilies`: profile of the signal tracks over the
   full-length copies of these subfamilies (copies that cover at least 90% of the RepeatMasker
   consensus), scaled to the same length from 5' to 3' with 2 kb flanks (`deepTools
   computeMatrix scale-regions` and `plotProfile`), per condition and per sample. The profile
   is the mean over many copies, so it shows, for example, accessibility at the 5' UTR promoter
   of young LINE-1s, even though single copies can't be resolved.

Output in `results/te/`:

| file                     | content                                                                 |
| ------------------------ | ----------------------------------------------------------------------- |
| `te_subfamily.tsv`       | per subfamily (repName), with its family and class                      |
| `te_family.tsv`          | per family (repFamily, e.g. L1, ERVK, Alu)                              |
| `te_class.tsv`           | per class (repClass, e.g. LINE, LTR, SINE)                              |
| `te_library_stats.tsv`   | read pairs kept, duplicates, multimappers and fraction of insertions in repeats |
| `bigwig/{sample}.cpm.bigwig` | Tn5 insertions per million insertions, including multimappers (see the caveat below) |
| `bigwig/mean/{condition}.cpm.bigwig` | mean of the replicates of a condition |
| `deseq2/{contrast}/differential_subfamilies.tsv` | per subfamily: family, class, copies, `baseMean`, `log2FoldChange` (and shrunken), `pvalue`, `padj`, `change` |
| `deseq2/{contrast}/{ma,volcano}_plot.png` | MA and volcano plots, most significant subfamilies labelled |
| `profile/{conditions,samples}.png` | profile over the full-length copies, one panel per subfamily (with its number of copies), one line per condition or sample |
| `deseq2/{contrast}/heatmap_classes.png` | heatmap faceted by class and condition: row z-scores of variance-stabilised insertions in square cells. By default the significant subfamilies of LINE, SINE, LTR and DNA (at most 25 per class); set with `te: heatmap` |

For each sample the `te_{subfamily,family,class}.tsv` tables have:

| column                       | meaning                                                                    |
| ---------------------------- | -------------------------------------------------------------------------- |
| `{sample}_count`             | insertions of all read pairs, multimappers at one of their best alignments |
| `{sample}_count_unique`      | insertions of uniquely mapping read pairs                                  |
| `{sample}_frac_unique`       | `count_unique / count`: low for young repeats, whose signal comes from multimappers and is missed by the ENCODE outputs |
| `{sample}_CPM`               | insertions per million insertions of the sample                            |
| `{sample}_log2_enrichment`   | log2 of the fraction of insertions in the repeat over the fraction of the genome it covers; > 0: more accessible than the genome average |

Notes:

- Counts are per repeat subfamily, family or class, not per copy: a multimapper is assigned to
  a random copy among its best alignments, which is close to the right subfamily total but not a
  reliable per-copy count.
- **Signal tracks:** a multimapper is shown at one randomly chosen copy among its best
  alignments, so the signal at a single copy of a young repeat (e.g. one L1MdTf element) does
  not show the accessibility of that copy. Only the pattern over many copies, such as a higher
  signal over most L1MdTf copies in one condition, is meaningful. In uniquely mappable regions
  the tracks show the same insertions as the ENCODE alignments, but as insertions per million
  instead of MACS2 fold enrichment or p-values. Use the subfamily tables, not the tracks, for
  quantification.
- Random assignment of multimappers gives the same subfamily counts and differential
  accessibility as keeping up to 100 alignments per read pair (`bowtie2 -k 100`, multimappers
  weighted 1/number of alignments), at about 1/6 of the alignment time: see
  [docs/te_alignment_comparison](docs/te_alignment_comparison/README.md), which reproduces the
  increased accessibility of young LINE-1s in Morc1 knockout gonocytes (Uneme et al., 2024).
- Longer reads resolve more copies; with short reads (< 50 bp) more of the signal is in
  multimappers.
- The TE alignment takes about as long as the ENCODE alignment.

## Validation against ENCODE

### Human (MCF-7, GRCh38)

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
| Overlap, optimal (ENCODE default) | 248,218 / 248,668         | 96%                      | 93%                        |
| IDR, optimal                      | 173,525 / 173,574         | 93%                      | 90%                        |
| IDR, conservative                 | 168,761 / 168,691         | 99.8%                    | 96%                        |

¹ Same coordinates and summit. The signal values and p-values of these peaks are identical
to ENCODE's (Pearson r = 1.000).

The remaining differences are expected from the different genome assembly (Ensembl primary
assembly instead of ENCODE's GRCh38 no-alt analysis set), the different random split into
pseudoreplicates and the newer pipeline and tool versions.

QC metrics (workflow / ENCODE):

| Metric                       | Replicate 1             | Replicate 2             |
| ---------------------------- | ----------------------- | ----------------------- |
| Alignment rate               | 0.955 / 0.955           | 0.974 / 0.974           |
| Fraction mitochondrial reads | 0.042 / 0.042           | 0.021 / 0.021           |
| Non-duplicate reads          | 61,376,290 / 61,376,458 | 56,295,148 / 56,295,774 |
| NRF                          | 0.8219 / 0.8218         | 0.8361 / 0.8361         |
| PBC1                         | 0.8301 / 0.8300         | 0.8440 / 0.8440         |
| PBC2                         | 6.13 / 6.13             | 6.70 / 6.70             |
| TSS enrichment²              | 29.3 / 29.7             | 26.9 / 27.2             |
| Fraction of fragments in NFR | 0.3844 / 0.3844         | 0.3759 / 0.3759         |
| NFR / mononucleosome         | 1.201 / 1.201           | 1.209 / 1.209           |
| Synthetic JS distance        | 0.613 / 0.614           | 0.608 / 0.608           |
| Fraction of reads in DHS     | 0.619 / 0.619           | 0.604 / 0.605           |
| ... in promoters             | 0.332 / 0.332           | 0.301 / 0.301           |
| ... in enhancers             | 0.349 / 0.349           | 0.370 / 0.370           |
| ... in exclusion list        | 0.0024 / 0.0024         | 0.0020 / 0.0020         |

² The small difference comes from the alignments and from the edge normalisation of the
profile.

Both replicates met the ENCODE standards for alignment rate, read depth, FRiP, TSS enrichment,
nucleosomal pattern (NFR and mononucleosome peaks present) and reproducibility. Their library
complexity was acceptable (mild PCR bottlenecking).

### Mouse (dendritic cells, mm10 and GRCm39)

The workflow was also run on ENCODE experiment
[ENCSR984HFU](https://www.encodeproject.org/experiments/ENCSR984HFU/) (mouse dendritic cells,
two paired-end replicates, 42 bp), once on mm10 (`mm38`, Ensembl release 102) for comparison
with ENCODE's processed files (pipeline v1.9.2), and once on GRCm39 (`mm39`, Ensembl release
116).

mm10 compared with ENCODE:

| Peak set / signal                        | Workflow / ENCODE | Identical to ENCODE peaks¹ | Pearson r (1 kb bins) |
| ---------------------------------------- | ----------------- | -------------------------- | --------------------- |
| Overlap, optimal (ENCODE default)        | 163,718 / 163,939 | 95.5%                      |                       |
| Overlap, conservative                    | 104,145 / 104,153 | 99.7%                      |                       |
| IDR, optimal                             | 102,595 / 102,490 | 91.6%                      |                       |
| IDR, conservative                        | 38,987 / 38,986   | 99.7%                      |                       |
| Pooled replicates, fold change / p-value |                   |                            | 1.000 / 0.999         |

Peak sets of true replicates (conservative) are 99.7% identical; the optimal peak sets, which
come from pseudoreplicates, differ more because the reads are split into pseudoreplicates at
random.

QC metrics (workflow / ENCODE):

| Metric                       | Replicate 1             | Replicate 2             |
| ---------------------------- | ----------------------- | ----------------------- |
| Fraction mitochondrial reads | 0.0885 / 0.0885         | 0.2072 / 0.2072         |
| Non-duplicate reads          | 77,816,024 / 77,834,648 | 40,499,660 / 40,502,024 |
| NRF                          | 0.8793 / 0.8792         | 0.6027 / 0.6027         |
| PBC1                         | 0.9049 / 0.9048         | 0.5916 / 0.5916         |
| PBC2                         | 11.95 / 11.94           | 2.162 / 2.162           |
| TSS enrichment               | 39.9 / 40.0             | 14.2 / 14.2             |
| Fraction of fragments in NFR | 0.5436 / 0.5436         | 0.7367 / 0.7367         |
| Synthetic JS distance        | 0.639 / 0.638           | 0.298 / 0.297           |
| Fraction of reads in DHS     | 0.744 / 0.745           | 0.406 / 0.407           |

The rescue and self-consistency ratios (overlap: 1.57 / 1.57 and 2.76 / 2.77; IDR: 2.63 / 2.63
and 4.15 / 4.14) agree as well. Replicate 2 has a low library complexity (moderate PCR
bottlenecking) and FRiP (0.14), so reproducibility is borderline (overlap) or failing (IDR),
for both the workflow and ENCODE.

GRCm39 gives the same results as mm10. After lifting the mm10 peaks to GRCm39, 99.6-99.9% of
the GRCm39 replicate and conservative peaks overlap an mm10 peak and 97-98% are identical (same
coordinates and summit, identical p-values); the optimal peak sets differ more (93% identical),
as the random split into pseudoreplicates differs. The signal tracks correlate with r = 0.999
(1 kb bins) and the QC metrics are within 0.5% (e.g. TSS enrichment 39.9 and 14.2, FRiP 0.608
and 0.144, 163,773 optimal peaks). For GRCm39, ENCODE's mm10 exclusion list and TSS file are lifted over;
ENCODE has no DHS, promoter and enhancer regions for it. Peak annotation does differ:
Ensembl 116 has more and longer transcripts than Ensembl 102, so more peaks are annotated as
promoter (49% instead of 40%) and fewer as distal intergenic (15% instead of 25%).

## Authors

- Niek Wit
  - University of Cambridge
  - [ORCID profile](https://orcid.org/0009-0002-4330-5333)
  - [Google Scholar](https://scholar.google.co.uk/citations?user=USu7NNcAAAAJ&hl=en)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
