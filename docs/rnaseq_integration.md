# Integration with RNA-seq

The differential accessibility of a contrast can be compared with RNA-seq differential expression
of the same contrast: are upregulated genes more accessible, and downregulated genes less? This
needs `deseq2: enable: true` and an RNA-seq DESeq2 results table per contrast in
`rnaseq: tables`. It runs for the ENCODE consensus peaks and, with `te: enable: true`, separately
for the peaks of the TE module and the repeat subfamilies (also with `te: only: true`).

The RNA-seq rows are used in two separate groups:

- **Genes (non-TE):** compared with the peaks near their TSS. Without `te: only`, with the
  ENCODE consensus peaks (uniquely mapped reads); with the TE module, also with the peaks of the
  TE module. These are genome-wide peaks of the alignment that keeps multimapping reads, not
  only peaks on TEs (most are at promoters or outside repeats), so they cover genes as well.
- **TE subfamilies:** compared with the accessibility of the subfamily (all copies) and of the
  peaks overlapping it.

## Configuration

```yaml
rnaseq:
  tables:
    # DESeq2 contrast of this workflow ({treatment}_vs_{control}): RNA-seq table
    ko_vs_wt: /path/to/rnaseq/ko_vs_wt_deseq2.tsv
    # or several tables, e.g. genes and TEs in separate files (combined)
    # ko_vs_wt: [/path/to/ko_vs_wt_genes.csv, /path/to/ko_vs_wt_te.csv]
  columns:
    gene_id: ""              # "": the first column (row names written by write.csv)
    log2fc: log2FoldChange
    padj: padj
    base_mean: baseMean
  id_type: auto              # ensembl, symbol or auto
  te_separator: ":"
  fdr: 0.05
  lfc_threshold: 0
  distance_bins: [1000, 10000, 50000, 100000]
  pair_lfc_threshold: 1
  tss_profile_flank: 3000    # TE module: TSS profile, bp on each side
```

- **Tables:** tab-separated, or comma-separated if the name ends in `.csv` (or `.csv.gz`), e.g.
  `write.csv(as.data.frame(res), "ko_vs_wt_deseq2.csv")`. The RNA-seq contrast must have the same
  direction as the ATAC-seq contrast (treatment over control).
- **Genes:** rows are matched to the genes of the workflow's genome annotation (Ensembl GTF), by
  gene ID (a version suffix such as `.5` is removed) or by gene name (names that belong to more
  than one gene are not used). With `id_type: auto`, whichever matches more rows is used. The
  RNA-seq data should use the same genome and preferably the same Ensembl release.
- **TEs** (TE module only): rows whose ID, up to the first `te_separator`, is a RepeatMasker
  subfamily, e.g. `L1MdTf_I:L1:LINE` as written by
  [TEtranscripts](https://doi.org/10.1093/bioinformatics/btv422) and TEcount, or `L1MdTf_I`. A
  table can have both genes and TEs, or they can be in separate tables (a list of paths, with the
  same columns). Gene rows are only used for the peak figures, TE rows only for the subfamily
  figures. Locus-level tables (one row per copy) are not supported.
- Rows that match no gene or subfamily are left out and listed in
  `logs/integration/rnaseq_prepare_{contrast}.log` (`logs/te/integration/` with `te: only`), with
  the numbers of matched rows. Check this log first if the results look empty. A gene or subfamily that is in the table more than once
  keeps the row with the lowest `padj`.
- Genes and subfamilies with `padj < fdr` and `|log2FoldChange| > lfc_threshold` are up- or
  downregulated. Rows without `padj` (not tested, e.g. removed by the independent filtering of
  DESeq2) are left out of all tests, so the background is the genes that could have been called
  differential.

## Genes and peaks

Each consensus peak (or TE peak) is assigned to the gene of its nearest TSS
([ChIPseeker](https://doi.org/10.1093/bioinformatics/btv145) annotation in
`differential_peaks.tsv`), if that TSS is at most the largest of `distance_bins` away. Per tested
gene, the peaks are summarised as:

| `atac_status`     | the gene has                                    |
| ----------------- | ----------------------------------------------- |
| `more accessible` | at least one more accessible peak, no less accessible peak |
| `less accessible` | at least one less accessible peak, no more accessible peak |
| `both`            | more and less accessible peaks                  |
| `unchanged`       | only peaks that do not change                   |
| `no peak`         | no assigned peak                                |

The analyses:

1. **Peak-gene pairs with a significant change in both** (`pairs.tsv`): pairs of a peak and its
   gene where the peak is differential (padj < `deseq2: fdr`) and the gene up- or downregulated
   (padj < `rnaseq: fdr`), both with |log2 fold change| > `pair_lfc_threshold` (default 1;
   shrunken log2 fold change for the peaks). Each pair is *concordant up* (more accessible, up),
   *concordant down* (less accessible, down), *more accessible, down* or *less accessible, up*.
   The larger of the two adjusted p-values (`padj_both`) is a p-value for a change in both
   (intersection-union test, [Berger 1982](https://doi.org/10.1080/00401706.1982.10487790)); the
   product of the two, used in some tutorials, is not a p-value. Two plots show the pairs:
   - **Scatter** (`scatter.pdf`): peak against gene log2 fold change, with Spearman's rho.
   - **Distance** (`distance.pdf`): the distances of the peaks to the TSS, per category: are
     concordant pairs more common at promoters? This shows where the pairs are, not whether
     they are more common than expected; that is the enrichment test below.
2. **Enrichment** (`enrichment.tsv`, `enrichment.pdf`): Fisher's exact tests of whether
   upregulated genes have more accessible peaks more often than the other tested genes (and up
   with less accessible, down with more and with less accessible), per distance bin and for all
   distances. Concordant changes are usually strongest at promoters and weaker for distal peaks,
   which are often assigned to the wrong gene. p-values are adjusted (Benjamini-Hochberg) over all
   tests of a contrast.
3. **TSS profile** (TE module, `genes/tss_profile.pdf`): Tn5 insertions per million (mean of the
   replicates of each condition, the signal tracks of the TE module, scaled as set by
   `te: bigwig_normalisation`) within
   `tss_profile_flank` of the TSS of the up-, down- and unchanged genes of the RNA-seq data, one
   panel per group, with the two conditions of the contrast as lines (deepTools computeMatrix
   reference-point and plotProfile). If accessibility follows expression, the treatment line is
   above the control line at upregulated genes and below it at downregulated genes; the unchanged
   genes show differences that apply to all genes (e.g. normalisation or composition). The TSS
   BED files of the groups and the deepTools matrix are next to it.
4. **Correlation** (`correlation.tsv`): Spearman correlation of the shrunken log2 fold change of
   each peak and the log2 fold change of its gene, per distance bin, for all peaks and for
   differential peaks only. Genome-wide correlations are typically weak (rho of about 0.1 to
   0.3), as many changes in accessibility do not change expression (e.g. poised enhancers) and
   many changes in expression do not need a change in accessibility.

## Repeat subfamilies (TE module)

The TE rows of the RNA-seq table are compared with two measures of the accessibility of each
subfamily:

- the DESeq2 test of the insertions in all its copies (`te/deseq2/{contrast}/differential_subfamilies.tsv`)
- the differential TE peaks whose largest repeat overlap is a copy of the subfamily
  (`te_peak_status`, as `atac_status` above)

The output has the correlation of the subfamily log2 fold changes (all subfamilies and per class),
the cumulative distributions of the RNA-seq log2 fold changes per accessibility group (Wilcoxon
tests against unchanged subfamilies: a curve to the right of that of the unchanged subfamilies
means higher expression) and Fisher's exact tests. With a few hundred subfamilies,
these tests have less power than those of the genes.

## Output

In `results/integration/rna-seq/{contrast}/` (genes against the ENCODE consensus peaks),
`results/te/integration/rna-seq/{contrast}/genes/` (genes against the peaks of the TE module) and
`results/te/integration/rna-seq/{contrast}/subfamilies/` (TE subfamilies):

| file                    | content                                                         |
| ----------------------- | --------------------------------------------------------------- |
| `genes.tsv`             | per tested gene: RNA-seq results, numbers of assigned, more and less accessible peaks, `atac_status`, the most significant peak (ID, distance to the TSS, shrunken log2 fold change, `padj`) and `concordance` (`concordant`: up with a more accessible peak or down with a less accessible peak; `discordant`: the opposite) |
| `concordant_genes.tsv`  | the concordant genes, by RNA-seq `padj`                         |
| `enrichment.tsv`, `enrichment.pdf` | Fisher's exact tests per distance bin: numbers of genes, fractions of the up- or downregulated and the other genes with a differential peak, odds ratio with 95% confidence interval, p and adjusted p. Tests without up- or downregulated genes or without genes with a differential peak are NA and not plotted; odds ratios of infinity or 0 (no genes in one cell of the table) are drawn as triangles at the top or bottom of the plot |
| `correlation.tsv`       | Spearman correlation per distance bin, of all and of the differential peaks |
| `pairs.tsv`             | peak-gene pairs with a significant change in both: peak and gene, distance to the TSS, both log2 fold changes and adjusted p-values, `category`, `combined_score` (product of the log2 fold changes) and `padj_both` (the larger adjusted p-value) |
| `scatter.pdf`, `distance.pdf` | plots of these pairs (see above)                    |
| `tss_profile.pdf`, `tss_profile.matrix.gz`, `tss_{up,down,unchanged}.bed` | TE module (`genes/`): Tn5 insertions per million around the TSS of the up-, down- and unchanged genes, with the deepTools matrix and the TSS of each group |
| `subfamilies/subfamilies.tsv` | per subfamily tested in both: ATAC-seq (`atac_*`) and RNA-seq (`rna_*`) results, numbers of TE peaks and `concordance` |
| `subfamilies/correlation.tsv`, `scatter.pdf`, `cdf.pdf`, `cdf_tests.tsv`, `enrichment.tsv` | correlation, scatter plot, cumulative distributions of the RNA-seq log2 fold changes per accessibility group with Wilcoxon tests, and Fisher's exact tests, for the subfamilies |

`rnaseq_table.tsv` in `results/integration/rna-seq/{contrast}/` (in
`results/te/integration/rna-seq/{contrast}/` with `te: only`, so that all output is in
`results/te/`) has the matched genes and subfamilies of the
RNA-seq table, with their up- or downregulation.
`input_tables/` next to it has unchanged copies of the RNA-seq tables of `rnaseq: tables`, so
the results keep the input they were made from.

## Example: Morc1 knockout gonocytes

These figures show the output of the integration on public data of P0 *Morc1* knockout and
heterozygous gonocytes ([Uneme et al., 2024](https://doi.org/10.1073/pnas.2317095121); GEO
[GSE235429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235429)). For the biology, see
the publication. The ATAC-seq data (2 knockout, 2 heterozygous replicates) were run with
`te: only: true` (GRCm39, Ensembl 115; see the
[Transposable elements example](transposable_elements.md#example-morc1-knockout-gonocytes)). The
RNA-seq data (3 + 3 replicates) were aligned with STAR to the same genome and annotation, and genes
and TEs were counted with TEtranscripts and tested with DESeq2. They gave two tables, genes and
TEs, both used for the contrast `KO_vs_het`:

```yaml
rnaseq:
  tables:
    KO_vs_het:
      - /path/to/KO_P0_vs_Het_P0_genes.csv
      - /path/to/KO_P0_vs_Het_P0_te.csv
  columns:
    gene_id: ensembl_gene_id
```

`logs/te/integration/rnaseq_prepare_KO_vs_het.log` lists how many rows matched a gene or a
RepeatMasker subfamily, and the numbers of up- and downregulated genes and subfamilies.

### Genes and TE module peaks

![Morc1: peak-gene pairs with a significant change in both](images/rnaseq_integration/morc1_genes_scatter.png)

**Scatter** (`genes/scatter.pdf`): each point is a peak and the gene of its nearest TSS (within
100 kb), both with a significant change. The x axis is the shrunken ATAC-seq log2 fold change of
the peak, the y axis the RNA-seq log2 fold change of the gene. The colours give the four
combinations (concordant up or down, and the two discordant ones), with their numbers in the
legend, and the caption gives Spearman's rho of the pairs.

![Morc1: distance of the pairs to the TSS](images/rnaseq_integration/morc1_genes_distance.png)

**Distance** (`genes/distance.pdf`): the distance of the peak to the TSS for the pairs of each
combination. The nearest-gene assignment is most reliable close to the TSS.

![Morc1: enrichment of differential peaks at differential genes](images/rnaseq_integration/morc1_genes_enrichment.png)

**Enrichment** (`genes/enrichment.pdf`): per distance bin, the odds ratio (log2, with 95%
confidence interval) that an up- or downregulated gene has a more or less accessible peak,
compared with the other tested genes (Fisher's exact test). Triangles are odds ratios of 0, from
bins without any gene in that combination.

![Morc1: Tn5 insertions around the TSS of up-, down- and unchanged genes](images/rnaseq_integration/morc1_genes_tss_profile.png)

**TSS profile** (`genes/tss_profile.pdf`): Tn5 insertions per million around the TSS of the up-,
down- and unchanged genes, as the mean of the replicates of each condition. The height of the two
lines depends on how the signal tracks are scaled, so compare each gene group with the unchanged
genes rather than the two lines with each other, and use the peak tests above for significance.
The same profile with tracks scaled by the DESeq2 size factors of the TE peaks
(`te: bigwig_normalisation: deseq2`):

![Morc1: TSS profile with tracks scaled by DESeq2 size factors](images/rnaseq_integration/morc1_genes_tss_profile_deseq2.png)

See [Scaling of the signal tracks](transposable_elements.md#scaling-of-the-signal-tracks) for when
to use which scaling.

### Repeat subfamilies

![Morc1: subfamily log2 fold changes in ATAC-seq and RNA-seq](images/rnaseq_integration/morc1_subfamilies_scatter.png)

**Scatter** (`subfamilies/scatter.pdf`): each point is a subfamily, with the shrunken log2 fold
change of the insertions in all its copies against its RNA-seq log2 fold change, coloured by class.
Subfamilies that are significant in both are solid, and those with the largest changes are
labelled. Spearman's rho per class is in `subfamilies/correlation.tsv`.

![Morc1: RNA-seq log2 fold changes of subfamilies per accessibility group](images/rnaseq_integration/morc1_subfamilies_cdf.png)

**Cumulative distributions** (`subfamilies/cdf.pdf`): the RNA-seq log2 fold changes of the
subfamilies, grouped by accessibility: on the left by the test of all copies, on the right by the
differential TE peaks in the subfamily. A curve to the right of the unchanged group means larger
RNA-seq log2 fold changes (knockout over heterozygous) than in the unchanged group. It does not by
itself mean that the group is upregulated, as both groups can have negative fold changes. The tests are in `subfamilies/cdf_tests.tsv` (Wilcoxon) and
`subfamilies/enrichment.tsv` (Fisher's exact test).

## Caveats

- **Nearest gene:** assigning a peak to the nearest TSS is right for most promoter peaks but
  often wrong for distal peaks, whose target genes can be hundreds of kb away. Use the distance
  bins to see where the relationship holds.
- **Fold changes:** the scatter plots and correlations use the shrunken ATAC-seq log2 fold
  changes, the RNA-seq log2 fold changes as given in the table (use shrunken ones there too, if
  available).
- **Same samples:** the analysis compares fold changes, so the RNA-seq and ATAC-seq data do not
  have to come from the same samples, but should be from the same conditions.
- **TEs:** the RNA-seq TE counts depend on how multimapping reads were assigned (e.g. the EM
  algorithm of TEtranscripts), as do the ATAC-seq subfamily counts (see
  [Transposable elements](transposable_elements.md#caveats)).
