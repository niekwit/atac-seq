"""Runs te_filter_multi.py (the TE module filter for bowtie2 -k and default
mode alignments, before the module switched to default mode) outside Snakemake"""

import os
import runpy
import sys
import types

bam, bam_out, stats, max_alignments, threads, log = sys.argv[1:]
snakemake = types.SimpleNamespace(
    input=[bam],
    output={"bam": bam_out, "stats": stats},
    params={"filter_chrs": ["MT"], "max_alignments": int(max_alignments)},
    threads=int(threads),
    log=[log],
)
runpy.run_path(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "te_filter_multi.py",
    ),
    init_globals={"snakemake": snakemake},
)
