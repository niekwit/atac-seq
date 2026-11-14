# QC and reproducibility

## QC summary

`results/qc/encode_qc_summary.tsv` has one row per sample with the QC metrics of the ENCODE
pipeline. Most are graded as `ideal`, `acceptable` or `concerning` against the
[ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/#standards):

| Column                 | Meaning                                                                                                   | Ideal        | Acceptable |
| ---------------------- | --------------------------------------------------------------------------------------------------------- | ------------ | ---------- |
| `alignment_rate`       | Fraction of read pairs aligned by bowtie2                                                                 | > 0.95       | > 0.80     |
| `frac_mito_reads`      | Fraction of aligned reads that align to the mitochondrial genome, as by ENCODE (not graded; lower is better) |              |            |
| `nodup_non_mito_reads` | Reads left after filtering and removal of duplicates and mitochondrial reads; these are used for peak calling | ≥ 50 million |            |
| `NRF`                  | Non-redundant fraction: distinct fragments / all fragments. Low values mean many PCR duplicates           | > 0.9        | > 0.8      |
| `PBC1`                 | PCR bottlenecking coefficient 1: positions with exactly one fragment / positions with at least one        | > 0.9        | > 0.8      |
| `PCR_bottlenecking`    | Bottlenecking by PBC1 as in the ENCODE QC report: none (> 0.9), mild (0.8-0.9), moderate (0.5-0.8), severe (< 0.5) |  |  |
| `PBC2`                 | PCR bottlenecking coefficient 2: positions with one fragment / positions with two                         | > 3          | > 1        |
| `num_peaks`            | Peaks of the replicate (MACS2 p < 0.01, top 300,000, blacklist filtered)                                  |              |            |
| `FRiP`                 | Fraction of Tn5 cutting sites in the replicate's peaks                                                    | > 0.3        | > 0.2      |
| `TSS_enrichment`       | Signal at transcription start sites relative to the background 2 kb away, calculated as by ENCODE        | genome dependent¹ |       |
| `TSS_enrichment_ataqv` | TSS enrichment calculated by ataqv (not graded, see below)                                                |              |            |
| `frac_fragments_in_NFR`, `NFR_mono_nuc_ratio` | Fraction of fragments < 150 bp (nucleosome free region, NFR), and NFR / mononucleosome (150-300 bp) fragments (ENCODE flags < 0.4 and < 2.5) | | |
| `NFR_peak`, `mono_nuc_peak`, `di_nuc_peak` | Peaks in the fragment length distribution at 20-90, 120-250 and 300-500 bp            |              |            |
| `nucleosome_pattern_standard` | ENCODE standard: a NFR and a mononucleosome peak are present                                       | both present |            |
| `synthetic_JS_distance` | Jensen-Shannon distance from uniform coverage (deepTools fingerprint); higher means more enriched signal |              |            |
| `frac_reads_in_DHS`, `_promoters`, `_enhancers`, `_blacklist` | Fraction of Tn5 cutting sites in ENCODE's DNase I hypersensitive sites, promoters, enhancers and exclusion list (DHS, promoters and enhancers for GRCh38 and mm10 only) | | |
| `AT_dropout`, `GC_dropout` | Picard GC bias: reads missing from AT- or GC-rich regions (%)                                         |              |            |

¹ hg38: > 7 ideal, > 5 acceptable; hg19: > 10 ideal, > 6 acceptable;
mm10/mm39: > 15 ideal, > 10 acceptable.

### TSS enrichment

`TSS_enrichment` follows the ENCODE pipeline (`encode_task_tss_enrich.py`): reads are centred
on their 5' end, their coverage is averaged over ±2 kb windows around the TSSs of protein
coding genes, and the maximum is divided by the coverage at the window edges.

The TSSs are those of ENCODE's TSS files (GRCh38, mm10; lifted over for mm39), as the score
depends on the annotation: in recent Ensembl releases the 5' end of many genes lies upstream of
the main TSS, which lowers the score (by ~15% for GRCm39, Ensembl 116), whereas the TSSs of
canonical transcripts raise it (by ~50% for GRCh38) above the scale of the ENCODE thresholds.
As in ENCODE, reads near the window edges are partly excluded, which raises the score by about a
third.

The ENCODE thresholds only apply to scores calculated this way: ataqv calculates TSS enrichment
on a different scale (about 5 where ENCODE reports about 29), so `TSS_enrichment_ataqv` is not
graded.

## Reproducibility

`results/qc/encode_reproducibility_summary.tsv` has one row per condition and reproducibility
method (IDR and overlap): the optimal and conservative peak sets, their number of peaks
(ENCODE: > 150,000 ideal, > 100,000 acceptable) and FRiP, and the rescue and self-consistency
ratios. A ratio above 2 means poor agreement between replicates (rescue ratio) or within a
replicate (self-consistency ratio); with one ratio above 2 reproducibility is `borderline`,
with both `fail`.

## Peaks reproducible in several replicates

The optimal and conservative peak sets come from a single pair of replicates (ENCODE picks the
pair with the most peaks), so with three or more replicates a peak only has to be reproducible
in two of them. `results/{idr,overlap}/{condition}/replicated_peak.narrowPeak` (not part of the
ENCODE pipeline) has the peaks of the pooled sample that are reproducible in at least *k* of the
*n* replicates of a condition, with *k* = `peaks: min_replicate_fraction` × *n*, rounded up and
at least 2 (default 1: all replicates):

- **Naive overlap:** the peak overlaps a peak of at least *k* replicates by at least 50% (of
  either peak), the criterion of the ENCODE naive overlap (`bedtools intersect -f 0.5 -F 0.5 -e`).
- **IDR:** IDR only compares pairs, so there must be *k* replicates of which every pair passed
  the IDR threshold for the peak; its score (column 5) is the lowest scaled IDR of these pairs.

With two replicates, or with *k* = *n*, the peaks are in the peak sets of all pairs of
replicates; with two replicates this is the conservative peak set.

Requiring all replicates makes the set smaller with every replicate added, and a single weaker
replicate limits it: the log (`logs/replicated_peaks/`) lists the peaks of each replicate
(overlap) or pair (IDR).

:::{note}
With different numbers of replicates per condition the sets are not equally strict, so compare
conditions with the [DESeq2 analysis](comparing_samples.md), not with these sets.
:::

As in the other peak sets, a region can be listed more than once, with different summits
(`--call-summits`).
