#!/usr/bin/env bash
# Compares bowtie2 alignment modes for the TE module: default mode --seed 42,
# --very-sensitive --seed 42, -k 10, -k 20 and -k 100 (reference). Each mode
# uses the same reads, trimming, filtering (te_filter.py) and featureCounts
# settings as the TE module.
#
# Usage: run_comparison.sh <analysis directory>
#   The analysis directory holds config/ (copied from this folder) and
#   reads/{sample}_R{1,2}_001.fastq.gz (download_reads.sh).
# Environment variables:
#   THREADS     threads for alignment and counting (default 12)
#   MODES       space-separated subset of the modes (default: all)
#   SNAKEMAKE   snakemake command used to build the genome, bowtie2 index and
#               RepeatMasker annotation (default: snakemake --sdm conda)
# Needs the tools of env.yaml on PATH.
set -euo pipefail

OUT=$(readlink -f "$1")
HERE=$(dirname "$(readlink -f "$0")")
REPO=$(readlink -f $HERE/../..)
THREADS=${THREADS:-12}
MODES=${MODES:-"default very_sensitive k10 k20 k100"}
SNAKEMAKE=${SNAKEMAKE:-"snakemake --sdm conda"}
EXTRA="-X2000 --mm"
C=$OUT/te_alignment_comparison
SAMPLES=$(tail -n +2 $OUT/config/samples.csv | cut -d, -f1)

# bowtie2 options and te_filter max_alignments (1 in default mode) per mode
declare -A OPTS=(
    [default]="--seed 42"
    [very_sensitive]="--very-sensitive --seed 42"
    [k10]="-k 10"
    [k20]="-k 20"
    [k100]="-k 100"
)
declare -A KMAX=([default]=1 [very_sensitive]=1 [k10]=10 [k20]=20 [k100]=100)

mkdir -p $C/reads $C/bam $C/counts $C/logs
cd $OUT

echo "## $(date '+%F %T') genome, bowtie2 index and RepeatMasker annotation"
$SNAKEMAKE -s $REPO/workflow/Snakefile --cores $THREADS \
    resources/bowtie2_index/index.1.bt2 resources/te/repeats.saf resources/te/subfamilies.tsv \
    > $C/logs/snakemake_resources.log 2>&1

echo "## $(date '+%F %T') adapter trimming as in the workflow"
for s in $SAMPLES; do
    [ -s $C/reads/${s}_R2.trim.fastq.gz ] && continue
    # Outputs are renamed only when complete, so that an interrupted run is redone
    cutadapt -j $THREADS -e 0.1 -m 5 -a CTGTCTCTTATA -A CTGTCTCTTATA \
        -o $C/reads/${s}_R1.tmp.fastq.gz -p $C/reads/${s}_R2.tmp.fastq.gz \
        reads/${s}_R1_001.fastq.gz reads/${s}_R2_001.fastq.gz \
        > $C/logs/cutadapt_${s}.log 2>&1
    mv $C/reads/${s}_R1.tmp.fastq.gz $C/reads/${s}_R1.trim.fastq.gz
    mv $C/reads/${s}_R2.tmp.fastq.gz $C/reads/${s}_R2.trim.fastq.gz
done

for mode in $MODES; do
    for s in $SAMPLES; do
        out=$C/bam/${s}.${mode}.filtered.bam
        [ -s $out ] && continue
        echo "## $(date '+%F %T') align $s $mode (${OPTS[$mode]})"
        /usr/bin/time -f "%e s wall, %M kB max RSS" -o $C/logs/time_${s}_${mode}.txt \
            bash -c "bowtie2 ${OPTS[$mode]} $EXTRA --threads $THREADS \
                -x resources/bowtie2_index/index \
                -1 $C/reads/${s}_R1.trim.fastq.gz -2 $C/reads/${s}_R2.trim.fastq.gz \
                2> $C/logs/bowtie2_${s}_${mode}.log \
                | samtools view -b -@ 2 -o $C/bam/${s}.${mode}.bam -"
        echo "## $(date '+%F %T') filter $s $mode"
        python $HERE/run_te_filter.py $C/bam/${s}.${mode}.bam $C/bam/${s}.${mode}.tmp.bam \
            $C/counts/${s}.${mode}.filter_stats.tsv ${KMAX[$mode]} 4 \
            $C/logs/te_filter_${s}_${mode}.log
        mv $C/bam/${s}.${mode}.tmp.bam $out
        rm $C/bam/${s}.${mode}.bam
    done
    echo "## $(date '+%F %T') count $mode"
    bams=$(for s in $SAMPLES; do echo -n "$C/bam/${s}.${mode}.filtered.bam "; done)
    opt="-p -F SAF --read2pos 5 --largestOverlap -s 0 -T $THREADS"
    featureCounts $opt -M --fraction -a resources/te/repeats.saf \
        -o $C/counts/${mode}.multi.txt $bams > $C/logs/featurecounts_${mode}.log 2>&1
    featureCounts $opt -a resources/te/repeats.saf \
        -o $C/counts/${mode}.unique.txt $bams >> $C/logs/featurecounts_${mode}.log 2>&1
done
echo "## $(date '+%F %T') done"
