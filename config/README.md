## Workflow overview

This workflow analyses paired-end ATAC-seq data following the
[ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) (v2.2.3)
and checks the results against the [ENCODE ATAC-seq data standards](https://www.encodeproject.org/atac-seq/).
Samples of the same condition are treated as biological replicates of one ENCODE experiment.

1. Download genome sequence and annotation from Ensembl, and the ENCODE exclusion list
   (blacklist) and annotated regions (GRCh38, mm10; other genomes: blacklist from AnnotationHub)
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
12. Fold enrichment and p-value signal tracks (`MACS2 bdgcmp`) of each replicate and condition
13. QC: library complexity (NRF, PBC1, PBC2), fraction of mitochondrial reads, FRiP,
    TSS enrichment (calculated as by ENCODE), fragment length distribution and nucleosomal
    pattern (`Picard CollectInsertSizeMetrics`), GC bias (`Picard CollectGcBiasMetrics`),
    fingerprint and Jensen-Shannon distance (`deepTools plotFingerprint`), fraction of reads in
    DHS, promoters, enhancers and blacklisted regions, rescue and self-consistency ratios,
    `ataqv`, `MultiQC`
14. Annotation of the optimal peak set (`ChIPseeker`)

Not included from the ENCODE pipeline: its optional (off by default) preseq, cross-correlation
and Roadmap comparison, its HTML/JSON QC report and the bigBed/starch/hammock peak formats.
The fraction of reads in DHS, promoters and enhancers is only calculated for GRCh38 and mm10
(`mm38`), for which ENCODE provides these regions.
As in ENCODE, TSS enrichment uses one TSS per protein coding gene (5' end of each Ensembl
protein coding gene).

## Running the workflow

### Input data

Paired-end reads go in `reads/{sample}_R1_001.fastq.gz` and `reads/{sample}_R2_001.fastq.gz`.
The sample sheet `config/samples.csv` has the following layout:

| sample | condition |
| ------ | --------- |
| WT_1   | wt        |
| WT_2   | wt        |
| KO_1   | ko        |
| KO_2   | ko        |

Sample names may only contain letters, numbers and underscores, and must end with `_[0-9]`.
The replicates of a condition are numbered (rep1, rep2, ...) in the order of the sample sheet.

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
| **idr**                           |                                                                            |                |
| enable                            | run IDR (naive overlap is always run)                                      | true           |
| threshold                         | IDR threshold                                                              | 0.05           |
| **pseudoreplication_random_seed** | seed for pseudoreplication (0: tagAlign file size, as ENCODE)              | 0              |

### Output

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `results/filtered/{sample}.nodup.bam`                  | filtered, deduplicated alignments without mitochondrial reads   |
| `results/tagalign/{sample}.tagAlign.gz`                | Tn5 shifted reads                                               |
| `results/macs2/{sample}.bfilt.narrowPeak`              | blacklist-filtered peaks of each replicate                      |
| `results/{idr,overlap}/{condition}/optimal_peak.narrowPeak` | optimal reproducible peak set (ENCODE's main output)       |
| `results/{idr,overlap}/{condition}/conservative_peak.narrowPeak` | conservative reproducible peak set                    |
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
