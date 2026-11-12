# bowtie2 alignment mode for the TE module

The TE module keeps multimapping reads, which make up most of the reads in young
transposable elements. This comparison tests whether bowtie2's default mode, which reports
one best alignment per read pair and picks one at random among equally good ones (reproducible
with `--seed 42`), gives the same per-subfamily counts as reporting up to 100 alignments per
read pair (`-k 100`) and weighting each by 1/NH. Default mode is several times faster.

## Data

ATAC-seq of FACS-sorted P0 gonocytes from Morc1 heterozygous (het, control) and knockout (KO)
mice ([Uneme et al., 2024, PNAS 121:e2317095121](https://doi.org/10.1073/pnas.2317095121);
GEO [GSE235429](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235429), BioProject
PRJNA986047). In Morc1 KO gonocytes chromatin is more accessible at young LINE1 elements
(Fig. 2J of the paper), which the comparison should reproduce in every mode.

| sample | SRA run | GEO sample | read pairs (2 × 151 bp) |
|---|---|---|---|
| P0_het_1 | SRR24979836 | GSM7502306 | 141.7M |
| P0_het_2 | SRR24979837 | GSM7502307 | 42.9M |
| P0_KO_1 | SRR24979838 | GSM7502308 | 98.5M |
| P0_KO_2 | SRR24979839 | GSM7502309 | 52.6M |

The first 10M read pairs of each run are used.

## Methods

- Reads: first 10M spots of each run (`fastq-dump -X`), adapters trimmed as in the workflow
  (`cutadapt -e 0.1 -m 5 -a CTGTCTCTTATA -A CTGTCTCTTATA`).
- Genome and repeats: mm39 (Ensembl, release in `config/config.yaml`) and the UCSC
  RepeatMasker annotation, built by the workflow.
- Alignment: bowtie2 2.5.4 (Langmead & Salzberg, 2012,
  [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)) with `-X2000 --mm` in five
  modes: default `--seed 42`, `--very-sensitive --seed 42`, `-k 10`, `-k 20` and `-k 100`.
- Filtering: `te_filter_multi.py`, the module's filter from before it switched to default
  mode (best properly paired alignments, duplicates removed, MT reads removed, NH tag = number
  of best alignments).
- Counting: featureCounts 2.1.1 (Liao et al., 2014,
  [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)),
  `-p -F SAF --read2pos 5 --largestOverlap -s 0`, with `-M --fraction` (all reads) and
  without (uniquely aligned reads only), per subfamily.
- Comparison with `-k 100` (`compare_modes.py`): correlation of per-subfamily counts,
  library-size-normalised log2 ratios, with young L1 and IAP subfamilies reported
  separately.
- Differential accessibility KO vs het per mode (`deseq2_modes.R`): DESeq2 (Love et al., 2014,
  [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)) on rounded
  per-subfamily counts.
- Run time and peak memory per mode with `/usr/bin/time`, 12 threads.

## Running it

```bash
# Tools (create the environment outside any synced folder)
mamba env create -f env.yaml -p /path/to/envs/te_alignment_comparison
conda activate /path/to/envs/te_alignment_comparison

# Analysis directory with the configuration of this folder
OUT=/path/to/analysis
mkdir -p $OUT && cp -r config $OUT/

./download_reads.sh $OUT 10000000
SNAKEMAKE="snakemake --sdm conda" THREADS=12 ./run_comparison.sh $OUT
mkdir -p $OUT/te_alignment_comparison/comparison
python compare_modes.py $OUT > $OUT/te_alignment_comparison/comparison/report.txt
Rscript deseq2_modes.R $OUT > $OUT/te_alignment_comparison/comparison/deseq2_report.txt
```

`deseq2_modes.R` needs DESeq2, e.g. from `workflow/envs/deseq2.yaml`. `run_comparison.sh` skips
finished steps, so it can be restarted, and `MODES="default k100"` limits it to some of the
modes.

## Results

Run on 2026-11-11/12 (mm39, Ensembl release 115; 10M read pairs per sample, 12 threads).

### Run time and peak memory per sample

| mode | minutes per sample | peak memory |
|---|---|---|
| default `--seed 42` | 15–19 | 3.4 GB |
| `--very-sensitive --seed 42` | 17–23 | 3.4 GB |
| `-k 10` | 18–23 | 3.5 GB |
| `-k 20` | 25–36 | 3.7 GB |
| `-k 100` | 88–137 | 5.0 GB |

### Counts per subfamily compared with `-k 100`

For the 760–816 subfamilies with at least 50 counts, the correlation of log2 counts with
`-k 100` is at least 0.9998 in every mode and sample, and no subfamily differs by more than
0.5 log2.

Library-size-normalised log2(count / `-k 100` count), averaged over the four samples:

| subfamily | default | very_sensitive | k10 | k20 |
|---|---|---|---|---|
| L1MdTf_I | 0.03 | 0.02 | 0.02 | 0.01 |
| L1MdTf_II | 0.02 | 0.03 | 0.01 | 0.01 |
| L1MdTf_III | 0.05 | 0.05 | 0.03 | 0.03 |
| L1MdA_I | 0.03 | 0.03 | 0.01 | 0.01 |
| L1MdA_II | 0.09 | 0.09 | 0.05 | 0.04 |
| L1MdA_III | 0.06 | 0.07 | 0.04 | 0.03 |
| L1MdGf_I | 0.04 | 0.04 | 0.02 | 0.02 |
| L1MdGf_II | 0.07 | 0.08 | 0.05 | 0.04 |
| L1MdN_I | 0.15 | 0.14 | 0.08 | 0.05 |
| L1MdF_I | 0.05 | 0.05 | 0.02 | 0.01 |
| L1MdF_II | 0.11 | 0.10 | 0.06 | 0.04 |
| IAPEz-int | 0.08 | 0.08 | 0.06 | 0.05 |
| IAPLTR1_Mm | 0.03 | 0.03 | 0.00 | -0.01 |
| IAPLTR1a_Mm | 0.07 | 0.07 | 0.06 | 0.06 |

The largest change among subfamilies with at least 200 counts in every sample is GSAT_MM
(satellite), −0.34 log2 in default mode.

Only 15–35% of the `-k 100` counts of L1MdTf_I–III and L1MdA_I–II come from uniquely
aligned reads, so counting only those would lose most of the signal in young L1s.

### Differential accessibility, Morc1 KO vs het

DESeq2 log2 fold changes:

| subfamily | default | very_sensitive | k10 | k20 | k100 |
|---|---|---|---|---|---|
| L1MdTf_I | 0.69 | 0.69 | 0.67 | 0.67 | 0.68 |
| L1MdTf_II | 0.68 | 0.68 | 0.69 | 0.69 | 0.69 |
| L1MdTf_III | 0.80 | 0.80 | 0.80 | 0.80 | 0.80 |
| L1MdA_I | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 |
| L1MdA_II | 0.85 | 0.85 | 0.85 | 0.85 | 0.86 |
| L1MdA_III | 0.58 | 0.59 | 0.58 | 0.59 | 0.59 |
| L1MdGf_I | 0.30 | 0.30 | 0.31 | 0.31 | 0.32 |
| L1MdGf_II | 0.51 | 0.51 | 0.52 | 0.52 | 0.53 |
| L1MdN_I | 0.24 | 0.24 | 0.25 | 0.27 | 0.28 |
| L1MdF_I | 0.17 | 0.18 | 0.18 | 0.18 | 0.18 |
| L1MdF_II | 0.32 | 0.33 | 0.33 | 0.34 | 0.36 |
| IAPEz-int | 0.09 | 0.09 | 0.09 | 0.10 | 0.11 |
| IAPLTR1_Mm | 1.77 | 1.77 | 1.77 | 1.76 | 1.79 |
| IAPLTR1a_Mm | 0.99 | 0.99 | 0.99 | 1.00 | 0.98 |

All of these are significant (adjusted p < 1e-5) in every mode except IAPEz-int, which is
adjusted p 0.02 with `-k 100` and 0.06–0.07 in the other modes. Accessibility rises most at
the youngest L1 subfamilies (L1MdA_I, L1MdTf, L1MdA_II), in line with Fig. 2J of the paper.

Agreement with `-k 100` over the 786 subfamilies with a mean normalised count of at least 50:

| mode | Pearson r of log2 fold changes | adjusted p < 0.05 | shared with `-k 100` (121) |
|---|---|---|---|
| default | 0.9979 | 125 | 120 |
| very_sensitive | 0.9980 | 122 | 118 |
| k10 | 0.9994 | 123 | 120 |
| k20 | 0.9997 | 123 | 120 |

### Conclusion

For per-subfamily counts, default mode with `--seed 42` gives the same result as `-k 100` at
about 1/6 of the alignment time and less memory. The TE module therefore uses default mode
with `--seed 42`. In this mode a multimapping read pair is assigned to one of its best
alignments at random, so per-copy (locus) counts are less reliable than subfamily totals.
