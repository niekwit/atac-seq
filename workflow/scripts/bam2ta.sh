#!/usr/bin/env bash
# Converts a paired-end BAM file to a Tn5 shifted tagAlign file
# (both mates as separate reads) as in the ENCODE ATAC-seq pipeline
# (encode_task_bam2ta.py)

set -euo pipefail

LOG=${snakemake_log[0]}
BAM=${snakemake_input[bam]}
OUTPUT=${snakemake_output[0]}
SUBSAMPLE=${snakemake_params[subsample]}
THREADS=${snakemake[threads]}

exec 2> "$LOG"

TA="$OUTPUT.tmp.tagAlign.gz"
trap 'rm -f "$TA" "$OUTPUT".nmsrt*' EXIT

# BAM -> BEDPE -> tagAlign
samtools sort -@ "$THREADS" -n -T "$OUTPUT.nmsrt" -o - "$BAM" |
    LC_COLLATE=C bedtools bamtobed -bedpe -mate1 -i stdin |
    awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' |
    gzip -nc > "$TA"

# Subsample read pairs (seeded with the tagAlign size, as ENCODE)
if [ "$SUBSAMPLE" -gt 0 ]; then
    SUBSAMPLED="$OUTPUT.tmp.subsampled.tagAlign.gz"
    zcat -f "$TA" | sed 'N;s/\n/\t/' |
        shuf -n $((SUBSAMPLE / 2)) --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f "$TA" | wc -c) -nosalt </dev/zero 2>/dev/null) |
        awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' |
        gzip -nc > "$SUBSAMPLED"
    mv "$SUBSAMPLED" "$TA"
fi

# Tn5 shift: +4 bp on the plus strand, -5 bp on the minus strand
zcat -f "$TA" |
    awk 'BEGIN {OFS = "\t"} {if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} if ($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1} } print $0}' |
    gzip -nc > "$OUTPUT"

if [ "$(zcat -f "$OUTPUT" | head -n 1 | wc -l)" -eq 0 ]; then
    echo "Empty tagAlign: $OUTPUT" >&2
    exit 1
fi
