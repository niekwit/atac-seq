# Comparing samples

The ENCODE pipeline analyses each condition on its own. To compare all samples, the workflow
also uses [DESeq2](https://doi.org/10.1186/s13059-014-0550-8) (`deseq2: enable: true`).

## Steps

1. **Consensus peaks:** the optimal peak sets of all conditions (`deseq2: peak_method`, overlap
   or IDR) are merged with `bedtools merge`.
2. **Counting:** Tn5 insertions (the 5' end of each read of `{sample}.nodup.bam`) in the
   consensus peaks are counted with `featureCounts`.
3. **PCA and sample correlation:** variance stabilising transformation (blind to the
   conditions) of the peaks with at least 10 insertions.
4. **Differential accessibility** between conditions, for each `[treatment, control]` pair in
   `deseq2: contrasts` or, if that is empty, each condition with at least 2 replicates against
   the reference condition (`reference` column of the
   [sample sheet](configuration.md#sample-sheet)): DESeq2 model `~ condition` fitted to all
   samples, Wald test, and log2 fold changes shrunken with
   [ashr](https://doi.org/10.1093/biostatistics/kxw041).

## Output

In `results/deseq2/`:

| file                                                  | content                                                     |
| ----------------------------------------------------- | ----------------------------------------------------------- |
| `consensus_peaks.bed`                                 | consensus peaks                                             |
| `consensus_peak_counts.txt`                           | Tn5 insertions per peak and sample                          |
| `vst_counts.tsv`                                      | variance stabilised counts                                  |
| `pca.png`, `pca.tsv`                                  | PCA of the most variable peaks (`deseq2: pca_ntop`), samples coloured by condition |
| `sample_correlation.png`, `sample_correlation.tsv`    | Pearson correlation of the samples over all peaks           |
| `consensus_peaks_annotated.txt`                       | consensus peaks annotated with `ChIPseeker`                 |
| `dds.rds`                                             | fitted DESeq2 object                                        |
| `{treatment}_vs_{control}/differential_peaks.tsv`     | peak coordinates, annotation and test results, by p-value   |
| `{treatment}_vs_{control}/ma_plot.png`, `volcano_plot.png` | plots of the shrunken log2 fold changes                |

In `differential_peaks.tsv`, `change` is `more accessible` or `less accessible` in the treatment
if `padj` is below `deseq2: fdr`. The p-values come from the unshrunken `log2FoldChange` (against
`deseq2: lfc_threshold`); use `log2FoldChange_shrunken` for ranking and plots.

## Normalisation

DESeq2 scales the samples with median-of-ratios size factors, calculated over the peaks.

:::{warning}
This assumes that most peaks do not change, so a global change in accessibility, or a large
difference in FRiP between samples, is not measured reliably. Check `FRiP` in the
[QC summary](qc.md#qc-summary) and the MA plot, in which the cloud of peaks should be centred
on 0.
:::
