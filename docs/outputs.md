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
| `qc/{sample}.fraglen_dist.pdf`, `.nucleosomal.qc`      | fragment length distribution and nucleosomal QC                 |
| `qc/{sample}.tss_enrich.pdf`                           | TSS enrichment profile                                          |
| `qc/{sample}.gc_plot.pdf`                              | GC bias                                                         |
| `qc/{condition}.jsd_plot.pdf`                          | fingerprint plot of the replicates                              |
| `ataqv_report/index.html`                              | ataqv report                                                    |

## Comparing samples (`deseq2: enable`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `deseq2/pca.pdf`, `sample_correlation.pdf`             | PCA and correlation of all samples                              |
| `deseq2/{treatment}_vs_{control}/differential_peaks.tsv` | differential accessibility between conditions                 |

All files are listed in [Comparing samples](comparing_samples.md#output).

## Transposable elements (`te: enable`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `te/te_{subfamily,family,class}.tsv`                   | Tn5 insertions in repeats                                       |
| `te/bigwig/{sample,mean/{condition}}.{cpm,deseq2}.bigwig` | Tn5 insertion tracks including multimappers, per sample and mean per condition (signal at single copies of young repeats is not reliable) |
| `te/deseq2/pca.pdf`, `te/deseq2/sample_correlation.pdf` | PCA and sample correlation from the repeat subfamily counts (with `deseq2: enable`) |
| `te/deseq2/{contrast}/differential_subfamilies.tsv`    | differential accessibility of repeat subfamilies (with `deseq2: enable`) |
| `te/deseq2/{contrast}/heatmap_classes.pdf`             | heatmap of repeat subfamilies, faceted by class (with `deseq2: enable`) |
| `te/profile/{conditions,samples}.pdf`                  | profile of the TE signal tracks over full-length copies of `te: profile: subfamilies` |
| `te/profile/copy_lengths.pdf`                          | length distribution of the full-length copies in the profile and of all copies, per entry of `te: profile: subfamilies` |
| `te/multiqc/multiqc_report.html`, `te/qc/`             | QC of the TE alignments (also with `te: only`)                  |
| `te/peaks/consensus_peaks_repeats.tsv`                 | TE consensus peaks with their repeat and fraction of uniquely mapping insertions |
| `te/peaks/deseq2/pca.pdf`, `te/peaks/deseq2/sample_correlation.pdf` | PCA and sample correlation from the TE peak counts (with `deseq2: enable`) |
| `te/peaks/deseq2/{contrast}/differential_peaks.tsv`    | differential accessibility of the TE peaks (with `deseq2: enable`) |
| `te/correlation/{contrast}/venn.pdf`                   | GAT enrichment and Venn diagrams of the differential TE peaks and `te: correlation: datasets` |

All files are listed in [Transposable elements](transposable_elements.md#output).

## Integration with RNA-seq (`rnaseq: tables`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `integration/rna-seq/{contrast}/genes.tsv`             | per gene: RNA-seq results and the differential peaks assigned to it |
| `integration/rna-seq/{contrast}/{scatter,distance,enrichment}.pdf` | peak-gene pairs with a significant change in both, enrichment of differential peaks at up- and downregulated genes |
| `te/integration/rna-seq/{contrast}/genes/`             | the same for the peaks of the TE module (with `te: enable`, also with `te: only`) |
| `te/integration/rna-seq/{contrast}/subfamilies/`       | TE subfamily expression against subfamily accessibility (with `te: enable`) |

All files are listed in [Integration with RNA-seq](rnaseq_integration.md#output).

## Integration with H3K9me3 (`h3k9me3: contrasts`)

| file                                                   | content                                                         |
| ------------------------------------------------------ | --------------------------------------------------------------- |
| `integration/h3k9me3/{contrast}/{enrichment,distribution}.pdf` | differential peaks in regions that lose or gain H3K9me3, log2 fold changes by H3K9me3 change |
| `integration/h3k9me3/{contrast}/scatter.pdf`           | peak against H3K9me3 log2 fold changes (DiffBind tables)        |
| `integration/h3k9me3/{contrast}/signal_profile.pdf`    | H3K9me3 signal around the differential peaks (with `bigwigs`)   |
| `te/integration/h3k9me3/{contrast}/peaks/`             | the same for the peaks of the TE module (with `te: enable`, also with `te: only`) |
| `te/integration/h3k9me3/{contrast}/atac_profile.pdf`   | Tn5 insertions over the lost, gained and unchanged H3K9me3 regions |
| `te/integration/h3k9me3/{contrast}/subfamilies/`       | repeat subfamilies and full-length copies in the H3K9me3 regions against their accessibility |

All files are listed in [Integration with H3K9me3](h3k9me3_integration.md#output).

## Integration with DNA methylation (`methylation: contrasts`)

The same files as for H3K9me3, in `integration/methylation/{contrast}/` and
`te/integration/methylation/{contrast}/`, with the DMRs as regions (hypo- and hypermethylated).
`scatter.pdf` needs a methylKit table, and `signal_profile.pdf` needs methylation BigWigs.
All files are listed in [Integration with DNA methylation](methylation_integration.md#output).

## Integration of all data sets (`multiomics: enable`)

| Path (in `results/`) | Content |
| -------------------- | ------- |
| `te/integration/multiomics/{contrast}/subfamilies/heatmap.pdf` | change in accessibility, expression, H3K9me3 and methylation of the significant repeat subfamilies |
| `te/integration/multiomics/{contrast}/subfamilies/profile_heatmap.pdf` | the same for the `te: profile: subfamilies` subfamilies |
| `te/integration/multiomics/{contrast}/peaks/upset.pdf` | UpSet plot of the TE peaks by accessibility, mark changes and nearest gene |
| `te/integration/multiomics/{contrast}/peaks/enrichment.pdf` | accessibility and expression by combination of the marks |
| `te/integration/multiomics/{contrast}/copies/{heatmap,profile}.pdf` | signal of all data sets over the full-length copies (with `te: profile: subfamilies`) |

All files are listed in [Integration of all data sets](multiomics_integration.md#output).

## Methods

`methods/methods.md` (in the analysis directory, not `results/`) is a methods section of the
steps that were run, written after all other outputs. It gives the configuration values and
cites each tool where it is first mentioned, with its version and its reference (PubMed ID,
else DOI, else URL). The cited references follow, then the configuration of the run
(including command line overrides). The versions are those installed in the conda
environments of the run. They are listed in `methods/software_versions.tsv` (columns
environment, package, version, source, cited) for the packages named in each environment
file and the Python and R interpreters; `cited` marks the packages cited in the text.
`methods/references.bib` has the cited references as BibTeX (e.g. for Zotero).

The methods are generated with [snakemake-methods](https://github.com/niekwit/snakemake-methods);
the text is in `workflow/methods/methods.md.j2`.
