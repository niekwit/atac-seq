#!/usr/bin/env bash
# Removes unmapped, low quality, orphan and multimapping reads as in
# the ENCODE ATAC-seq pipeline (encode_task_filter.py)

set -euo pipefail

LOG=${snakemake_log[0]}
BAM=${snakemake_input[0]}
OUTPUT=${snakemake_output[0]}
MULTIMAPPING=${snakemake_params[multimapping]}
MAPQ=${snakemake_params[mapq_thresh]}
ASSIGN_MULTIMAPPERS=${snakemake_params[assign_multimappers]}
THREADS=${snakemake[threads]}

exec 2> "$LOG"

TMP_FILT="$OUTPUT.tmp_filt.bam"
FIXMATE="$OUTPUT.fixmate.bam"
trap 'rm -f "$TMP_FILT" "$FIXMATE"' EXIT

if [ "$MULTIMAPPING" -gt 0 ]; then
    # Keep reads with at most MULTIMAPPING alignments (no MAPQ filter)
    samtools view -F 524 -f 2 -u "$BAM" |
        samtools sort -@ "$THREADS" -n -T "$OUTPUT.sort1" -o "$TMP_FILT" -
    samtools view -h "$TMP_FILT" |
        python3 "$ASSIGN_MULTIMAPPERS" -k "$MULTIMAPPING" --paired-end |
        samtools fixmate -r - "$FIXMATE"
else
    samtools view -F 1804 -f 2 -q "$MAPQ" -u "$BAM" |
        samtools sort -@ "$THREADS" -n -T "$OUTPUT.sort1" -o "$TMP_FILT" -
    samtools fixmate -r "$TMP_FILT" "$FIXMATE"
fi

# Keep only properly paired primary alignments
samtools view -F 1804 -f 2 -u "$FIXMATE" |
    samtools sort -@ "$THREADS" -T "$OUTPUT.sort2" -o "$OUTPUT" -

if [ "$(samtools view -c "$OUTPUT")" -eq 0 ]; then
    echo "No reads left after filtering $BAM" >&2
    exit 1
fi
