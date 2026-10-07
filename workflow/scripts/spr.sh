#!/usr/bin/env bash
# Splits a paired-end tagAlign file into two pseudoreplicates by randomly
# assigning read pairs, as in the ENCODE ATAC-seq pipeline (encode_task_spr.py)

set -euo pipefail

LOG=${snakemake_log[0]}
TA=${snakemake_input[0]}
PR1=${snakemake_output[pr1]}
PR2=${snakemake_output[pr2]}
SEED=${snakemake_params[seed]}

exec 2> "$LOG"

PREFIX="$PR1.split"
trap 'rm -f "$PREFIX".*' EXIT

NLINES=$(( ($(zcat -f "$TA" | wc -l) / 2 + 1) / 2 ))

if [ "$SEED" -eq 0 ]; then
    SEED=$(zcat -f "$TA" | wc -c)
fi
echo "Random seed for pseudoreplication: $SEED" >&2

zcat -f "$TA" | sed 'N;s/\n/\t/' |
    shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null) |
    split -d -l "$NLINES" - "$PREFIX."

for pair in "00 $PR1" "01 $PR2"; do
    set -- $pair
    awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$PREFIX.$1" |
        gzip -nc > "$2"
done
