#!/usr/bin/env bash
# Generates fold enrichment and p-value bigWigs from MACS2 bedGraphs as in the
# ENCODE ATAC-seq pipeline (encode_task_macs2_signal_track_atac.py)

set -euo pipefail

LOG=${snakemake_log[0]}
TREAT=${snakemake_input[treat]}
CONTROL=${snakemake_input[control]}
TA=${snakemake_input[ta]}
CHRSZ=${snakemake_input[chrom_sizes]}
FC_BW=${snakemake_output[fc]}
PVAL_BW=${snakemake_output[pval]}

exec 2> "$LOG"

PREFIX="${FC_BW%.fc.signal.bigwig}"
trap 'rm -f "${PREFIX}"_FE.bdg "${PREFIX}"_ppois.bdg "${PREFIX}".*.bedgraph' EXIT

to_bigwig() {
    local bdg=$1 bedgraph=$2 bigwig=$3
    bedtools slop -i "$bdg" -g "$CHRSZ" -b 0 | bedClip stdin "$CHRSZ" "$bedgraph"
    # Sort and remove overlapping regions
    LC_COLLATE=C sort -k1,1 -k2,2n "$bedgraph" |
        awk 'BEGIN{OFS="\t"}{if (NR==1 || NR>1 && (prev_chr!=$1 || prev_chr==$1 && prev_chr_e<=$2)) {print $0}; prev_chr=$1; prev_chr_e=$3;}' > "$bedgraph.srt"
    bedGraphToBigWig "$bedgraph.srt" "$CHRSZ" "$bigwig"
    rm -f "$bedgraph" "$bedgraph.srt"
}

# Fold enrichment
macs2 bdgcmp -t "$TREAT" -c "$CONTROL" --o-prefix "$PREFIX" -m FE
to_bigwig "${PREFIX}_FE.bdg" "${PREFIX}.fc.bedgraph" "$FC_BW"

# -log10(p-value); sval is the number of tags per million
SVAL=$(zcat -f "$TA" | wc -l | awk '{print $1 / 1000000.0}')
macs2 bdgcmp -t "$TREAT" -c "$CONTROL" --o-prefix "$PREFIX" -m ppois -S "$SVAL"
to_bigwig "${PREFIX}_ppois.bdg" "${PREFIX}.pval.bedgraph" "$PVAL_BW"
