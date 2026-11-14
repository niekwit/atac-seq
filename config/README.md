## Workflow overview

This workflow analyses paired-end ATAC-seq data following the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) (v2.2.3)
and checks the results against the [ENCODE ATAC-seq data standards](https://www.encodeproject.org/atac-seq/).
Samples of the same condition are treated as biological replicates of one ENCODE experiment.

1. Download genome sequence and annotation from Ensembl, and the ENCODE exclusion list
   (blacklist), TSS file and annotated regions (GRCh38, mm10; mm39: ENCODE mm10 exclusion list
   and TSS file lifted to mm39; other genomes: blacklist from AnnotationHub)
2. Quality control of reads (`FastQC`)
3. Adapter trimming (`cutadapt -e 0.1 -m 5`, Nextera adapter)
4. Alignment (`bowtie2 -k 5 -X2000 --mm`)
5. Removal of unmapped, unpaired and multimapping (> 4 alignments) reads (`samtools`, `pysam`)
6. Duplicate marking (`Picard MarkDuplicates`) and removal of duplicates and mitochondrial reads
7. Conversion to Tn5 shifted (+4/-5 bp) tagAlign files of both mates
8. Pseudoreplication of each replicate, and pooling of the replicates and of the pseudoreplicates per condition
9. Peak calling (`MACS2 callpeak -p 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits`),
   keeping the top 300,000 peaks
10. Removal of peaks in blacklisted regions and on non-standard chromosomes
11. Reproducibility: IDR (threshold 0.05) and naive overlap of all pairs of replicates,
    of the pseudoreplicates of each replicate and of the pooled pseudoreplicates,
    giving the optimal and conservative peak set of each condition
    and, not part of the ENCODE pipeline, the peaks reproducible in several replicates
12. Fold enrichment and p-value signal tracks (`MACS2 bdgcmp`) of each replicate and condition
13. QC: library complexity (NRF, PBC1, PBC2), fraction of mitochondrial reads, FRiP,
    TSS enrichment (calculated as by ENCODE), fragment length distribution and nucleosomal
    pattern (`Picard CollectInsertSizeMetrics`), GC bias (`Picard CollectGcBiasMetrics`),
    fingerprint and Jensen-Shannon distance (`deepTools plotFingerprint`), fraction of reads in
    DHS, promoters, enhancers and blacklisted regions, rescue and self-consistency ratios,
    `ataqv`, `MultiQC`
14. Annotation of the optimal peak set (`ChIPseeker`)
15. Comparison of all samples (not part of the ENCODE pipeline): PCA, sample correlation and
    differential accessibility between conditions from Tn5 insertions in the merged optimal peaks
    (`DESeq2`)
16. Optional (not part of the ENCODE pipeline): accessibility of transposable elements per
    RepeatMasker subfamily, family and class, keeping multimapping reads
    (see [Transposable elements](../README.md#transposable-elements))

Not included from the ENCODE pipeline: its optional (off by default) preseq, cross-correlation
and Roadmap comparison, its HTML/JSON QC report and the bigBed/starch/hammock peak formats.
The fraction of reads in DHS, promoters and enhancers is only calculated for GRCh38 and mm10
(`mm38`), for which ENCODE provides these regions.
TSS enrichment uses ENCODE's TSS files (one TSS per protein coding gene; GRCh38, mm10, and for
mm39 the mm10 file lifted to mm39), so that scores can be compared with the ENCODE thresholds.
For other genomes it uses the 5' end of each Ensembl protein coding gene.

## Running the workflow

### Input data

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

### Parameters

The defaults in `config/config.yaml` are those of the ENCODE pipeline.

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

### Output

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `results/filtered/{sample}.nodup.bam`                  | filtered, deduplicated alignments without mitochondrial reads   |
| `results/tagalign/{sample}.tagAlign.gz`                | Tn5 shifted reads                                               |
| `results/macs2/{sample}.bfilt.narrowPeak`              | blacklist-filtered peaks of each replicate                      |
| `results/{idr,overlap}/{condition}/optimal_peak.narrowPeak` | optimal reproducible peak set (ENCODE's main output)       |
| `results/{idr,overlap}/{condition}/conservative_peak.narrowPeak` | conservative reproducible peak set                    |
| `results/{idr,overlap}/{condition}/replicated_peak.narrowPeak` | peaks reproducible in at least `min_replicate_fraction` of the replicates (conditions with replicates) |
| `results/{idr,overlap}/{condition}/optimal_peak_annotated.txt` | annotated optimal peak set                              |
| `results/bigwig/{sample,pooled/condition}.{fc,pval}.signal.bigwig` | fold enrichment and -log10(p-value) signal tracks   |
| `results/qc/encode_qc_summary.tsv`                     | QC metrics of each sample, graded against the ENCODE standards  |
| `results/qc/encode_reproducibility_summary.tsv`        | reproducibility of each condition, graded against the standards |
| `results/multiqc/multiqc_report.html`                  | FastQC, cutadapt, bowtie2, Picard, deepTools and samtools stats report |
| `results/qc/{sample}.fraglen_dist.png`, `.nucleosomal.qc` | fragment length distribution and nucleosomal QC             |
| `results/qc/{sample}.tss_enrich.png`                   | TSS enrichment profile                                          |
| `results/qc/{sample}.gc_plot.png`                      | GC bias                                                         |
| `results/qc/{condition}.jsd_plot.png`                  | fingerprint plot of the replicates                              |
| `results/ataqv_report/index.html`                      | ataqv report                                                    |
| `results/deseq2/pca.png`, `sample_correlation.png`     | PCA and correlation of all samples                              |
| `results/deseq2/{treatment}_vs_{control}/differential_peaks.tsv` | differential accessibility between conditions         |
| `results/te/te_{subfamily,family,class}.tsv`           | Tn5 insertions in repeats (TE module)                           |
| `results/te/bigwig/{sample,mean/{condition}}.cpm.bigwig` | Tn5 insertion tracks including multimappers, per sample and mean per condition (TE module; signal at single copies of young repeats is not reliable) |
| `results/te/deseq2/pca.png`, `results/te/deseq2/sample_correlation.png` | PCA and sample correlation from the repeat subfamily counts (TE module with `deseq2: enable`) |
| `results/te/deseq2/{contrast}/differential_subfamilies.tsv` | differential accessibility of repeat subfamilies (TE module with `deseq2: enable`) |
| `results/te/deseq2/{contrast}/heatmap_classes.png` | heatmap of repeat subfamilies, faceted by class (TE module with `deseq2: enable`) |
| `results/te/profile/{conditions,samples}.png` | profile of the TE signal tracks over full-length copies of `te: profile: subfamilies` |
| `results/te/multiqc/multiqc_report.html`, `results/te/qc/` | QC of the TE alignments (also with `te: only`) |
| `results/te/peaks/consensus_peaks_repeats.tsv` | TE consensus peaks with their repeat and fraction of uniquely mapping insertions |
| `results/te/peaks/deseq2/pca.png`, `results/te/peaks/deseq2/sample_correlation.png` | PCA and sample correlation from the TE peak counts (with `deseq2: enable`) |
| `results/te/peaks/deseq2/{contrast}/differential_peaks.tsv` | differential accessibility of the TE peaks (with `deseq2: enable`) |
| `results/te/correlation/{contrast}/venn.png` | GAT enrichment and Venn diagrams of the differential TE peaks and `te: correlation: datasets` |
