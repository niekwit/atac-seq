#!/usr/bin/env bash
# Calls peaks on a tagAlign file with MACS2 as in the ENCODE ATAC-seq pipeline
# (encode_task_macs2_atac.py): peaks are sorted by p-value, renamed, capped at
# cap_num_peak and clipped to the chromosome sizes. The pileup and lambda
# bedGraphs are kept for the signal tracks.

set -euo pipefail

LOG=${snakemake_log[0]}
TA=${snakemake_input[ta]}
CHRSZ=${snakemake_input[chrom_sizes]}
NPEAK=${snakemake_output[narrowpeak]}
TREAT=${snakemake_output[treat]}
CONTROL=${snakemake_output[control]}
GSIZE=${snakemake_params[gsize]}
PVAL=${snakemake_params[pval_thresh]}
SMOOTH_WIN=${snakemake_params[smooth_win]}
CAP=${snakemake_params[cap_num_peak]}

exec 2> "$LOG"

OUTDIR=$(dirname "$NPEAK")
NAME=$(basename "$NPEAK" .narrowPeak)
PREFIX="$OUTDIR/$NAME"
SHIFT=$(awk -v w="$SMOOTH_WIN" 'BEGIN{x = w / 2; printf "%d", -(x == int(x) ? x : int(x + 0.5))}')
trap 'rm -f "${PREFIX}_peaks.xls" "${PREFIX}_summits.bed" "${PREFIX}_peaks.narrowPeak" "$NPEAK.tmp"' EXIT

macs2 callpeak \
    -t "$TA" -f BED -n "$NAME" --outdir "$OUTDIR" -g "$GSIZE" -p "$PVAL" \
    --shift "$SHIFT" --extsize "$SMOOTH_WIN" \
    --nomodel -B --SPMR --keep-dup all --call-summits

# Keep the top CAP peaks in awk, not with head, which would end the pipe
# early (SIGPIPE) when there are more peaks than CAP
LC_COLLATE=C sort -k 8gr,8gr "${PREFIX}_peaks.narrowPeak" |
    awk -v cap="$CAP" 'BEGIN{OFS="\t"} NR<=cap {$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) $10=$2+int(($3-$2+1)/2.0); print $0}' > "$NPEAK.tmp"

# Clip peaks to 0-chromSize
bedClip "$NPEAK.tmp" "$CHRSZ" "$NPEAK" -truncate -verbose=2

if [ "${PREFIX}_treat_pileup.bdg" != "$TREAT" ]; then
    mv "${PREFIX}_treat_pileup.bdg" "$TREAT"
    mv "${PREFIX}_control_lambda.bdg" "$CONTROL"
fi
