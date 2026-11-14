# Configuration

The full documentation, including the workflow steps and output files, is at
<https://atac-seq.readthedocs.io>.

<!-- docs-include-start -->

## Sample sheet

Paired-end reads go in `reads/{sample}_R1_001.fastq.gz` and `reads/{sample}_R2_001.fastq.gz`.
The sample sheet `config/samples.csv` has the following layout:

| sample | condition | reference |
| ------ | --------- | --------- |
| WT_1   | wt        | yes       |
| WT_2   | wt        | yes       |
| KO_1   | ko        | no        |
| KO_2   | ko        | no        |

Sample names may only contain letters, numbers and underscores, and must end with `_[0-9]`.
The replicates of a condition are numbered (rep1, rep2, ...) in the order of the sample sheet.

`reference` marks the samples of the reference (control) condition with `yes`, and all other
samples with `no`. It is the reference level of the DESeq2 model. Without `deseq2: contrasts`,
each other condition with at least 2 replicates is compared to it (positive log2 fold changes:
more accessible than in the reference); the reference condition then also needs at least 2
replicates. The column is needed for the DESeq2 analysis of more than one condition, unless all
comparisons are listed in `deseq2: contrasts`.

## Parameters

`config/config.yaml` sets the parameters below. Set at least the genome (`genome: ensembl` and
`release`) and the effective genome size (`macs2: gsize`); the other defaults are those of the
ENCODE pipeline.

| parameter                         | details                                                                    | default        |
| --------------------------------- | -------------------------------------------------------------------------- | -------------- |
| **genome**                        |                                                                            |                |
| ensembl                           | genome build: hg19, hg38, mm38, mm39 or dm6                                | mm39           |
| release                           | Ensembl release                                                            | 115            |
| **mito_chr_name**                 | name of the mitochondrial chromosome                                       | MT             |
| **cutadapt**                      |                                                                            |                |
| adapter_r1, adapter_r2            | adapter sequences of read 1 and read 2                                     | CTGTCTCTTATA   |
| extra                             | other cutadapt arguments                                                   | -e 0.1 -m 5    |
| **bowtie2**                       |                                                                            |                |
| multimapping                      | maximum number of alignments of a read; 0: filter on MAPQ instead          | 4              |
| extra                             | other bowtie2 arguments                                                    | -X2000 --mm    |
| **filter**                        |                                                                            |                |
| mapq_thresh                       | minimum MAPQ (only used when `multimapping` is 0)                          | 30             |
| filter_chrs                       | chromosomes whose reads are removed after deduplication                    | [MT]           |
| **subsample_reads**               | number of reads to subsample each replicate to (0: none)                   | 0              |
| **macs2**                         |                                                                            |                |
| gsize                             | effective genome size (hs, mm or a number)                                 | mm             |
| pval_thresh                       | p-value threshold                                                          | 0.01           |
| smooth_win                        | smoothing window (`--extsize`; `--shift` is -smooth_win/2)                 | 150            |
| cap_num_peak                      | maximum number of peaks                                                    | 300000         |
| **peaks**                         |                                                                            |                |
| keep_chr_regex                    | peaks on chromosomes not matching this (POSIX extended) regex are removed  | [0-9]+\|X\|Y   |
| min_replicate_fraction            | fraction of the replicates in which a peak of `replicated_peak.narrowPeak` is reproducible (rounded up, at least 2) | 1 |
| **idr**                           |                                                                            |                |
| enable                            | run IDR (naive overlap is always run)                                      | true           |
| threshold                         | IDR threshold                                                              | 0.05           |
| **pseudoreplication_random_seed** | seed for pseudoreplication (0: tagAlign file size, as ENCODE)              | 0              |
| **deseq2**                        |                                                                            |                |
| enable                            | compare all samples with DESeq2                                            | true           |
| peak_method                       | optimal peak sets merged into consensus peaks: overlap or idr              | overlap        |
| pca_ntop                          | number of most variable peaks used for the PCA                             | 500            |
| contrasts                         | [treatment, control] pairs; empty: each condition against the reference    | []             |
| fdr                               | adjusted p-value below which a peak is differential                        | 0.05           |
| lfc_threshold                     | test against this absolute log2 fold change (0: any change)                | 0              |
| **te**                            |                                                                            |                |
| enable                            | run the transposable element module                                        | false          |
| only                              | run only the TE module (with the genome, index and read trimming it needs), not the ENCODE peaks, signal tracks, QC and DESeq2 of the peaks | false |
| classes                           | RepeatMasker classes to count (empty: all repeats)                         | LINE, SINE, LTR, DNA, Retroposon, RC, Satellite |
| normalisation                     | library size for the DESeq2 comparison of subfamilies: `deseq2` (median-of-ratios over the subfamilies) or `insertions` (all kept insertions of a sample) | deseq2 |
| heatmap: classes                  | list of classes in the heatmap of each contrast; `["all"]` for every class | ["LINE", "SINE", "LTR", "DNA"] |
| heatmap: subfamilies              | list of subfamilies in the heatmap; `["all"]` for the significant subfamilies of the selected classes. Named subfamilies are plotted whether significant or not (significant ones marked with *) | ["all"] |
| heatmap: max_per_class            | most significant subfamilies per class in the heatmap (with `subfamilies: ["all"]`) | 25 |
| profile: subfamilies              | subfamilies for the profile over their full-length copies (empty: no profile) | [] |
| profile: min_consensus_fraction   | fraction of the RepeatMasker consensus a copy must cover to count as full-length | 0.9 |
| profile: flank                    | bp before and after the copies in the profile                              | 2000 |
| profile: body_length              | length (bp) to which the copies are scaled in the profile                  | 6000 |
| correlation: datasets             | data sets to compare with the differential TE peaks with GAT, as name: BED file (e.g. EM-seq DMRs; UCSC or Ensembl chromosome names) | {} |
| correlation: gat_samples          | number of random placements of the peaks in GAT                            | 10000 |
