# Workflow steps

The workflow follows the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
(v2.2.3) and checks the results against the
[ENCODE ATAC-seq data standards](https://www.encodeproject.org/atac-seq/). Samples of the same
condition are treated as biological replicates of one ENCODE experiment.

1. Download genome sequence and annotation from Ensembl, and the ENCODE exclusion list
   (blacklist), TSS file and annotated regions (GRCh38, mm10; mm39: ENCODE mm10 exclusion list
   and TSS file lifted to mm39; other genomes: blacklist from AnnotationHub)
2. Quality control of reads (`FastQC`)
3. Adapter trimming (`cutadapt -e 0.1 -m 5`, Nextera adapter)
4. Alignment (`bowtie2 -k 5 -X2000 --mm`)
5. Removal of unmapped, unpaired and multimapping (> 4 alignments) reads (`samtools`, `pysam`)
6. Duplicate marking (`Picard MarkDuplicates`) and removal of duplicates and mitochondrial reads
7. Conversion to Tn5 shifted (+4/-5 bp) tagAlign files of both mates
8. Pseudoreplication of each replicate, and pooling of the replicates and of the
   pseudoreplicates per condition
9. Peak calling (`MACS2 callpeak -p 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all
   --call-summits`), keeping the top 300,000 peaks
10. Removal of peaks in blacklisted regions and on non-standard chromosomes
11. Reproducibility: IDR (threshold 0.05) and naive overlap of all pairs of replicates, of the
    pseudoreplicates of each replicate and of the pooled pseudoreplicates, giving the optimal and
    conservative peak set of each condition and, not part of the ENCODE pipeline, the
    [peaks reproducible in several replicates](qc.md#peaks-reproducible-in-several-replicates)
12. Fold enrichment and p-value signal tracks (`MACS2 bdgcmp`) of each replicate and condition
13. QC: library complexity (NRF, PBC1, PBC2), fraction of mitochondrial reads, FRiP, TSS
    enrichment (calculated as by ENCODE), fragment length distribution and nucleosomal pattern
    (`Picard CollectInsertSizeMetrics`), GC bias (`Picard CollectGcBiasMetrics`), fingerprint and
    Jensen-Shannon distance (`deepTools plotFingerprint`), fraction of reads in DHS, promoters,
    enhancers and blacklisted regions, rescue and self-consistency ratios, `ataqv`, `MultiQC`
    (see [QC and reproducibility](qc.md))
14. Annotation of the optimal peak set (`ChIPseeker`)
15. Comparison of all samples (not part of the ENCODE pipeline): PCA, sample correlation and
    differential accessibility between conditions from Tn5 insertions in the merged optimal
    peaks (`DESeq2`; see [Comparing samples](comparing_samples.md))
16. Optional (not part of the ENCODE pipeline): accessibility of transposable elements per
    RepeatMasker subfamily, family and class, keeping multimapping reads (see
    [Transposable elements](transposable_elements.md))

## Differences from the ENCODE pipeline

- Not included: its optional (off by default) preseq, cross-correlation and Roadmap comparison,
  its HTML/JSON QC report and the bigBed/starch/hammock peak formats.
- The fraction of reads in DHS, promoters and enhancers is only calculated for GRCh38 and mm10
  (`mm38`), for which ENCODE provides these regions.
- TSS enrichment uses ENCODE's TSS files (one TSS per protein coding gene; GRCh38, mm10, and for
  mm39 the mm10 file lifted to mm39), so that scores can be compared with the ENCODE thresholds.
  For other genomes it uses the 5' end of each Ensembl protein coding gene.
