# Output files

The main output files, all in `results/`. The pages on [QC](qc.md),
[comparing samples](comparing_samples.md) and [transposable elements](transposable_elements.md)
describe their contents in more detail.

## ENCODE pipeline

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `filtered/{sample}.nodup.bam`                          | filtered, deduplicated alignments without mitochondrial reads   |
| `tagalign/{sample}.tagAlign.gz`                        | Tn5 shifted reads                                               |
| `macs2/{sample}.bfilt.narrowPeak`                      | blacklist-filtered peaks of each replicate                      |
| `{idr,overlap}/{condition}/optimal_peak.narrowPeak`    | optimal reproducible peak set (ENCODE's main output)            |
| `{idr,overlap}/{condition}/conservative_peak.narrowPeak` | conservative reproducible peak set                            |
| `{idr,overlap}/{condition}/replicated_peak.narrowPeak` | peaks reproducible in at least `min_replicate_fraction` of the replicates (conditions with replicates; not part of the ENCODE pipeline) |
| `{idr,overlap}/{condition}/optimal_peak_annotated.txt` | annotated optimal peak set                                      |
| `bigwig/{sample,pooled/condition}.{fc,pval}.signal.bigwig` | fold enrichment and -log10(p-value) signal tracks           |
| `qc/encode_qc_summary.tsv`                             | QC metrics of each sample, graded against the ENCODE standards  |
| `qc/encode_reproducibility_summary.tsv`                | reproducibility of each condition, graded against the standards |
| `multiqc/multiqc_report.html`                          | FastQC, cutadapt, bowtie2, Picard, deepTools and samtools stats report |
| `qc/{sample}.fraglen_dist.png`, `.nucleosomal.qc`      | fragment length distribution and nucleosomal QC                 |
| `qc/{sample}.tss_enrich.png`                           | TSS enrichment profile                                          |
| `qc/{sample}.gc_plot.png`                              | GC bias                                                         |
| `qc/{condition}.jsd_plot.png`                          | fingerprint plot of the replicates                              |
| `ataqv_report/index.html`                              | ataqv report                                                    |

## Comparing samples (`deseq2: enable`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `deseq2/pca.png`, `sample_correlation.png`             | PCA and correlation of all samples                              |
| `deseq2/{treatment}_vs_{control}/differential_peaks.tsv` | differential accessibility between conditions                 |

All files are listed in [Comparing samples](comparing_samples.md#output).

## Transposable elements (`te: enable`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `te/te_{subfamily,family,class}.tsv`                   | Tn5 insertions in repeats                                       |
| `te/bigwig/{sample,mean/{condition}}.cpm.bigwig`       | Tn5 insertion tracks including multimappers, per sample and mean per condition (signal at single copies of young repeats is not reliable) |
| `te/deseq2/pca.png`, `te/deseq2/sample_correlation.png` | PCA and sample correlation from the repeat subfamily counts (with `deseq2: enable`) |
| `te/deseq2/{contrast}/differential_subfamilies.tsv`    | differential accessibility of repeat subfamilies (with `deseq2: enable`) |
| `te/deseq2/{contrast}/heatmap_classes.png`             | heatmap of repeat subfamilies, faceted by class (with `deseq2: enable`) |
| `te/profile/{conditions,samples}.png`                  | profile of the TE signal tracks over full-length copies of `te: profile: subfamilies` |
| `te/multiqc/multiqc_report.html`, `te/qc/`             | QC of the TE alignments (also with `te: only`)                  |
| `te/peaks/consensus_peaks_repeats.tsv`                 | TE consensus peaks with their repeat and fraction of uniquely mapping insertions |
| `te/peaks/deseq2/pca.png`, `te/peaks/deseq2/sample_correlation.png` | PCA and sample correlation from the TE peak counts (with `deseq2: enable`) |
| `te/peaks/deseq2/{contrast}/differential_peaks.tsv`    | differential accessibility of the TE peaks (with `deseq2: enable`) |
| `te/correlation/{contrast}/venn.png`                   | GAT enrichment and Venn diagrams of the differential TE peaks and `te: correlation: datasets` |

All files are listed in [Transposable elements](transposable_elements.md#output).
