"""
Removes unmapped, low quality, orphan and multimapping reads as in the
ENCODE ATAC-seq pipeline (encode_task_filter.py):

1. Without multimapping: keep properly paired reads with MAPQ >= mapq_thresh.
   With multimapping: keep properly paired reads with at most multimapping
   alignments per mate (no MAPQ filter), as assign_multimappers.py of the
   ENCODE pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline,
   MIT License, Copyright (c) 2017 ENCODE DCC)
2. samtools fixmate -r on the query name sorted reads
3. Keep only properly paired primary alignments, sorted by coordinate

samtools view/sort/fixmate are run through pysam.
"""

import itertools
import logging
import os

import pysam

# Load Snakemake variables
bam = snakemake.input[0]
output = snakemake.output[0]
multimapping = int(snakemake.params["multimapping"])
mapq_thresh = int(snakemake.params["mapq_thresh"])
threads = str(snakemake.threads)
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

flagged = f"{output}.tmp_flagged.bam"
name_sorted = f"{output}.tmp_filt.bam"
multimap_filtered = f"{output}.tmp_multimap.bam"
fixmate = f"{output}.fixmate.bam"


def samtools(command, *args):
    """Runs a samtools command through pysam"""
    logging.info(f"Running: samtools {command} {' '.join(args)}")
    # catch_stdout=False, as output goes to a file (-o)
    getattr(pysam, command)(*args, catch_stdout=False)


def remove_multimappers(name_sorted_bam, out_bam, cutoff):
    """
    Writes the reads of a query name sorted BAM file with at most cutoff
    alignments; reads with more alignments are removed altogether.
    """
    kept = removed = 0
    with pysam.AlignmentFile(name_sorted_bam, "rb") as fh_in, pysam.AlignmentFile(
        out_bam, "wbu", template=fh_in
    ) as fh_out:
        for _, alignments in itertools.groupby(fh_in, key=lambda a: a.query_name):
            alignments = list(alignments)
            if len(alignments) <= cutoff:
                for alignment in alignments:
                    fh_out.write(alignment)
                kept += 1
            else:
                removed += 1
    logging.info(f"Kept {kept} reads, removed {removed} multimapping reads")


try:
    if multimapping > 0:
        # -F 524: remove unmapped, mate unmapped, QC fail; -f 2: proper pairs
        samtools("view", "-F", "524", "-f", "2", "-u", "-o", flagged, bam)
        samtools(
            "sort", "-@", threads, "-n", "-T", f"{output}.sort1",
            "-o", name_sorted, flagged,
        )  # fmt: skip
        os.remove(flagged)
        # Both mates of a pair have at most multimapping alignments
        remove_multimappers(name_sorted, multimap_filtered, multimapping * 2)
        os.remove(name_sorted)
        samtools("fixmate", "-r", multimap_filtered, fixmate)
        os.remove(multimap_filtered)
    else:
        # -F 1804: also remove secondary alignments and duplicates
        samtools(
            "view", "-F", "1804", "-f", "2", "-q", str(mapq_thresh),
            "-u", "-o", flagged, bam,
        )  # fmt: skip
        samtools(
            "sort", "-@", threads, "-n", "-T", f"{output}.sort1",
            "-o", name_sorted, flagged,
        )  # fmt: skip
        os.remove(flagged)
        samtools("fixmate", "-r", name_sorted, fixmate)
        os.remove(name_sorted)

    # Keep only properly paired primary alignments
    samtools("view", "-F", "1804", "-f", "2", "-u", "-o", flagged, fixmate)
    os.remove(fixmate)
    samtools("sort", "-@", threads, "-T", f"{output}.sort2", "-o", output, flagged)
    os.remove(flagged)

    n_reads = int(pysam.view("-c", output))
    if n_reads == 0:
        raise ValueError(f"No reads left after filtering {bam}")
    logging.info(f"{n_reads} alignments left after filtering")
except Exception:
    logging.exception("Filtering failed")
    for tmp in [flagged, name_sorted, multimap_filtered, fixmate]:
        if os.path.exists(tmp):
            os.remove(tmp)
    raise

logging.info("Done!")
