#!/usr/bin/env bash
# Downloads the first PAIRS read pairs of each SRA run in config/sra_runs.csv to
# reads/{sample}_R{1,2}_001.fastq.gz in the analysis directory.
# Usage: download_reads.sh <analysis directory> [pairs, default 10000000]
# Needs fastq-dump (sra-tools) and pigz on PATH (env.yaml).
set -euo pipefail

OUT=$(readlink -f "$1")
PAIRS=${2:-10000000}

mkdir -p $OUT/reads $OUT/logs
tail -n +2 $OUT/config/sra_runs.csv | while IFS=, read -r sample run geo; do
    [ -s $OUT/reads/${sample}_R2_001.fastq.gz ] && continue
    echo "## $(date '+%F %T') $sample ($run, $geo): first $PAIRS read pairs"
    tmp=$OUT/reads/tmp_$run
    mkdir -p $tmp
    fastq-dump -X $PAIRS --split-files -O $tmp $run > $OUT/logs/fastq-dump_$sample.log 2>&1
    # Renamed only when complete, so that an interrupted run is redone
    for r in 1 2; do
        pigz -p 4 -c $tmp/${run}_$r.fastq > $tmp/R$r.fastq.gz
    done
    mv $tmp/R1.fastq.gz $OUT/reads/${sample}_R1_001.fastq.gz
    mv $tmp/R2.fastq.gz $OUT/reads/${sample}_R2_001.fastq.gz
    rm -r $tmp
done
echo "## $(date '+%F %T') done"
