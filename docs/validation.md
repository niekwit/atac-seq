# Validation against ENCODE

The workflow was run on ENCODE experiments and its output was compared with the files ENCODE
processed with its own ATAC-seq pipeline.

## Human (MCF-7, GRCh38)

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

## Mouse (dendritic cells, mm10 and GRCm39)

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
