# Fold enrichment and p-value signal tracks from the MACS2 pileups of
# each replicate and each pooled condition (ENCODE: macs2_signal_track)
# -----------------------------------------------------
rule macs2_signal_track:
    input:
        treat="results/macs2/{prefix}_treat_pileup.bdg",
        control="results/macs2/{prefix}_control_lambda.bdg",
        ta="results/tagalign/{prefix}.tagAlign.gz",
        chrom_sizes="resources/chrom_sizes.txt",
    output:
        fc="results/bigwig/{prefix}.fc.signal.bigwig",
        pval="results/bigwig/{prefix}.pval.signal.bigwig",
    threads: 1
    log:
        "logs/macs2_signal_track/{prefix}.log",
    conda:
        "../envs/encode.yaml"
    script:
        "../scripts/macs2_signal_track.sh"
