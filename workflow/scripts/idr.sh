#!/usr/bin/env bash
# Runs IDR on two peak sets as in the ENCODE ATAC-seq pipeline
# (encode_task_idr.py), keeping peaks with IDR <= threshold

set -euo pipefail

LOG=${snakemake_log[0]}
PEAK1=${snakemake_input[peak1]}
PEAK2=${snakemake_input[peak2]}
POOLED=${snakemake_input[pooled]}
CHRSZ=${snakemake_input[chrom_sizes]}
NPEAK=${snakemake_output[narrowpeak]}
UNTHRESHOLDED=${snakemake_output[unthresholded]}
PLOT=${snakemake_output[plot]}
IDR_LOG=${snakemake_output[idr_log]}
THRESH=${snakemake_params[threshold]}

exec 2> "$LOG"

IDR_OUT="$NPEAK.unthresholded-peaks.txt"
trap 'rm -f "$IDR_OUT" "$IDR_OUT".*' EXIT

idr --samples "$PEAK1" "$PEAK2" --peak-list "$POOLED" \
    --input-file-type narrowPeak --output-file "$IDR_OUT" \
    --rank p.value --soft-idr-threshold "$THRESH" \
    --plot --use-best-multisummit-IDR --log-output-file "$IDR_LOG"
mv "$IDR_OUT.png" "$PLOT"

# Clip peaks to 0-chromSize
bedClip "$IDR_OUT" "$CHRSZ" "$IDR_OUT.clipped" -truncate -verbose=2

# Column 12 is -log10(global IDR); sort by p-value (column 8)
NEG_LOG10_THRESH=$(awk -v t="$THRESH" 'BEGIN{print -log(t)/log(10)}')
awk -v t="$NEG_LOG10_THRESH" 'BEGIN{OFS="\t"} $12>=t {if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$IDR_OUT.clipped" |
    sort | uniq | sort -grk8,8 | cut -f 1-10 > "$NPEAK"

gzip -nc "$IDR_OUT.clipped" > "$UNTHRESHOLDED"

if [ ! -s "$NPEAK" ]; then
    echo "WARNING: no IDR peaks found. The IDR threshold might be too stringent or replicates have very poor concordance." >&2
fi
