# Transposable elements

Young transposable elements (TEs), such as L1Md_T/A/Gf LINE-1s and IAPs in mouse or L1HS and
L1PA2/3 in human, have many nearly identical copies. Their reads align to many places, so the
ENCODE filter (at most 4 alignments, one of which is kept at random) removes them, and they
hardly show up in the peaks and signal tracks. Older, diverged copies map uniquely and are
analysed like the rest of the genome.

The optional TE module (`te: enable: true` in `config/config.yaml`) measures the accessibility
of repeats per subfamily, family and class instead of per copy. It does not change the ENCODE
outputs. See the [parameters](configuration.md#parameters) under `te` for its settings.

## TE module only

With `te: only: true`, only the TE module runs (with the genome, bowtie2 index and read trimming
it needs), not the ENCODE peaks, signal tracks, QC and DESeq2 of the peaks.

QC of the TE alignments is always made (`results/te/qc/`, `results/te/multiqc/`): FastQC,
cutadapt, bowtie2 and `samtools stats` in a MultiQC report, and TSS enrichment, fragment length
distribution (nucleosome-free fraction) and GC bias as in the ENCODE QC. These alignments
include multimappers, so the values are close to, but not the same as, the ENCODE QC.

## Steps

### Alignment and counting

1. RepeatMasker annotation from UCSC (`rmsk`), with Ensembl chromosome names.
2. A second alignment of the trimmed reads that keeps multimappers (`bowtie2 --seed 42 -X2000
   --mm`): bowtie2's default mode reports one of the best alignments of each read pair, picked
   at random (reproducibly, with the seed) among equally good ones.
3. Properly paired read pairs are kept; PCR duplicates and reads on `filter_chrs` (MT) are
   removed. Read pairs in which a mate has an equally good alternative alignment (bowtie2
   `XS` ≥ `AS`) that the pairing does not resolve (MAPQ < 2) are marked as multimappers.
4. Tn5 insertions (5' end of each read) in repeats are counted with `featureCounts`, once with
   all read pairs (`-M`) and once with uniquely mapping read pairs only.

### Signal tracks and profiles

5. Signal tracks of Tn5 insertions (5' end of each read, 10 bp bins, per million insertions;
   `deepTools bamCoverage`) per sample, and the mean of the replicates of each condition
   (`deepTools bigwigAverage`). Unlike the ENCODE signal tracks, they include multimappers.
6. With subfamilies in `te: profile: subfamilies`: profile of the signal tracks over the
   full-length copies of these subfamilies (copies that cover at least 90% of the RepeatMasker
   consensus), scaled to the same length from 5' to 3' with 2 kb flanks (`deepTools
   computeMatrix scale-regions` and `plotProfile`), per condition and per sample. The profile
   is the mean over many copies, so it shows, for example, accessibility at the 5' UTR promoter
   of young LINE-1s, even though single copies can't be resolved.

### Differential accessibility of subfamilies

7. With `deseq2: enable`: differential accessibility of subfamilies with `DESeq2`, for the same
   contrasts, `fdr` and `lfc_threshold` as the peaks (Wald test, ashr-shrunken log2 fold changes
   for plots), and a PCA and sample correlation of the subfamily counts.

   Library sizes come from the median-of-ratios over the subfamilies
   (`te: normalisation: deseq2`), which assumes that most subfamilies do not change. If a
   condition changes the accessibility of many repeat families at once, use
   `te: normalisation: insertions` (all kept insertions of each sample) instead.

### TE peaks

8. Peaks from the TE alignments, so that open regions in young repeats are included: Tn5
   shifted tagAlign files and MACS2 as for the main peaks, but with the total chromosome length
   as effective genome size, since multimappers can align anywhere. The peaks of the pooled
   replicates that overlap a peak in at least `peaks: min_replicate_fraction` of the replicates
   by ≥ 50% are kept, and the peaks of all conditions are merged into consensus peaks.

   Each consensus peak gets the repeat with the largest overlap and `frac_unique`, the fraction
   of its insertions from uniquely mapping read pairs. In a peak with a low `frac_unique` the
   signal comes from multimappers, which sit at a random copy: the repeat subfamily is usually
   right, but the copy (locus) is uncertain.

   With `deseq2: enable`: differential accessibility of the consensus peaks with DESeq2, as for
   the main peaks, a PCA and sample correlation, and the peaks more and less accessible in the
   treatment as BED files.

### Comparison with other data sets

9. With data sets in `te: correlation: datasets` (name: BED file, e.g. EM-seq DMRs):
   [GAT](https://doi.org/10.1093/bioinformatics/btt343) tests per contrast whether the more and
   the less accessible TE peaks overlap the regions of each data set more than expected.

   GAT places the differential peaks at random within the TE consensus peaks (`gat_samples`
   times), so they are compared with accessible regions in general, not with the whole genome
   (open chromatin overlaps most regulatory data sets anyway). q-values are over all tests of a
   contrast.

   The result is shown as area-proportional Venn diagrams
   ([eulerr](https://CRAN.R-project.org/package=eulerr)) with the GAT fold enrichment and
   q-value. The areas count merged regions, as in
   [ChIPpeakAnno](https://doi.org/10.1186/1471-2105-11-237): overlapping peaks and regions are
   merged, so a peak that overlaps several regions counts once.

## Output

In `results/te/`:

| file                     | content                                                                 |
| ------------------------ | ----------------------------------------------------------------------- |
| `te_subfamily.tsv`       | per subfamily (repName), with its family and class                      |
| `te_family.tsv`          | per family (repFamily, e.g. L1, ERVK, Alu)                              |
| `te_class.tsv`           | per class (repClass, e.g. LINE, LTR, SINE)                              |
| `te_library_stats.tsv`   | read pairs kept, duplicates, multimappers and fraction of insertions in repeats |
| `qc/te_qc_summary.tsv`   | per sample: read pairs, duplicates, multimappers, insertions in repeats, TSS enrichment, fraction of reads in nucleosome-free regions |
| `multiqc/multiqc_report.html` | MultiQC report of the TE alignments, with the QC summary in the general statistics |
| `bigwig/{sample}.cpm.bigwig` | Tn5 insertions per million insertions, including multimappers (see the caveats below) |
| `bigwig/mean/{condition}.cpm.bigwig` | mean of the replicates of a condition |
| `deseq2/pca.png`, `deseq2/sample_correlation.png` | PCA (samples coloured by condition) and Pearson correlation of all samples from the variance-stabilised subfamily counts, as for the ENCODE peaks; data in `pca.tsv` and `sample_correlation.tsv` |
| `deseq2/{contrast}/differential_subfamilies.tsv` | per subfamily: family, class, copies, `baseMean`, `log2FoldChange` (and shrunken), `pvalue`, `padj`, `change` |
| `deseq2/{contrast}/{ma,volcano}_plot.png` | MA and volcano plots, most significant subfamilies labelled |
| `deseq2/{contrast}/heatmap_classes.png` | heatmap faceted by class and condition: row z-scores of variance-stabilised insertions in square cells. By default the significant subfamilies of LINE, SINE, LTR and DNA (at most 25 per class); set with `te: heatmap` |
| `profile/{conditions,samples}.png` | profile over the full-length copies, one panel per subfamily (with its number of copies), one line per condition or sample |
| `peaks/overlap/{condition}/replicated_peak.narrowPeak` | TE peaks of a condition, replicated in `min_replicate_fraction` of its replicates |
| `peaks/consensus_peaks_repeats.tsv` | per consensus peak: insertions, `frac_unique`, overlapping repeat (subfamily, family, class, bp) |
| `peaks/consensus_peaks_frac_unique.bed` | consensus peaks as BED6, score: `frac_unique` × 1000 |
| `peaks/deseq2/pca.png`, `peaks/deseq2/sample_correlation.png` | PCA and sample correlation from the TE consensus peak counts |
| `peaks/deseq2/{contrast}/differential_peaks.tsv` | differential accessibility of the TE peaks, with gene annotation, repeat and `frac_unique` |
| `peaks/deseq2/{contrast}/{more,less}_accessible.bed` | differential TE peaks as BED6 (fold changes in `differential_peaks.tsv`) |
| `correlation/{contrast}/gat.tsv` | per direction and data set: numbers of peaks and regions and of those overlapping, GAT observed and expected overlap (bp), fold, p- and q-value |
| `correlation/{contrast}/venn.png` | Venn diagrams of the differential peaks and each data set, counted as merged regions (overlapping peaks and regions merged, as in ChIPpeakAnno), with the numbers of peaks and regions and the GAT result |

For each sample the `te_{subfamily,family,class}.tsv` tables have:

| column                       | meaning                                                                    |
| ---------------------------- | -------------------------------------------------------------------------- |
| `{sample}_count`             | insertions of all read pairs, multimappers at one of their best alignments |
| `{sample}_count_unique`      | insertions of uniquely mapping read pairs                                  |
| `{sample}_frac_unique`       | `count_unique / count`: low for young repeats, whose signal comes from multimappers and is missed by the ENCODE outputs |
| `{sample}_CPM`               | insertions per million insertions of the sample                            |
| `{sample}_log2_enrichment`   | log2 of the fraction of insertions in the repeat over the fraction of the genome it covers; > 0: more accessible than the genome average |

## Caveats

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
  [TE alignment mode](te_alignment_comparison/README.md), which reproduces the increased
  accessibility of young LINE-1s in Morc1 knockout gonocytes (Uneme et al., 2024).
- Longer reads resolve more copies; with short reads (< 50 bp) more of the signal is in
  multimappers.
- The TE alignment takes about as long as the ENCODE alignment.
