"""
Converts a paired-end BAM file to a Tn5 shifted tagAlign file, as in the
ENCODE ATAC-seq pipeline (encode_task_bam2ta.py).

Why
---
Downstream steps (pseudoreplication, pooling, MACS2, FRiP) work on reads
as plain intervals, so the alignments are converted to tagAlign: a BED6
file with one line per read (chrom, start, end, N, 1000, strand). Both
mates are written as separate reads, on consecutive lines, so that the
pseudoreplicates can keep the mates of a fragment together.

The Tn5 transposase inserts the sequencing adapters as a dimer, with the
two cut sites 9 bp apart; the ends of the reads are therefore shifted
(+4 bp on the plus strand, -5 bp on the minus strand) to mark the centre
of the insertion, which is what the peak caller should see.

What
----
1. The BAM file is sorted by read name (samtools sort -n) so that the
   mates of each pair are next to each other. The sorted reads are streamed
   uncompressed into pysam, rather than written to a temporary file
2. Each pair is converted to two tagAlign lines, read 1 first (as
   `bedtools bamtobed -bedpe -mate1`); reads without a mate are skipped
3. If subsample_reads > 0, that many reads (half as many pairs) are
   randomly selected. The random seed is the size of the uncompressed
   tagAlign file, as in ENCODE
4. Tn5 shift: start + 4 on the plus strand, end - 5 on the minus strand;
   reads that would become empty are reduced to 1 bp at the read end
   (plus strand) or start (minus strand)

Differences from ENCODE:
- ENCODE subsamples with `shuf` and an openssl keystream as random source;
  this script uses Python's random number generator instead, so a different
  (equally random) subset is selected. Subsampled pairs are written in the
  original order, rather than in random order.
- Without subsampling, the Tn5 shift is applied while converting, rather
  than in a second pass over the tagAlign file.

Speed: on 20 million reads this takes about as long as the ENCODE shell
pipeline. Sorting dominates; output is compressed by pigz in a separate
process, as Python's gzip module would double the conversion time.
"""

import gzip
import logging
import os
import random
import subprocess

import pysam

# Load Snakemake variables
bam = snakemake.input["bam"]
output = snakemake.output[0]
subsample = int(snakemake.params["subsample"])
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

unshifted = f"{output}.tmp.tagAlign.gz"


class PigzWriter:
    """
    Context manager that gzips everything written to it into a file with
    pigz (-n: no file name and time stamp in the header, so that reruns
    give identical files)
    """

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.out_fh = open(self.path, "wb")
        self.log_fh = open(log, "a")
        self.process = subprocess.Popen(
            ["pigz", "-p", threads, "-n", "-c"],
            stdin=subprocess.PIPE,
            stdout=self.out_fh,
            stderr=self.log_fh,
        )
        return self.process.stdin

    def __exit__(self, *exc):
        self.process.stdin.close()
        returncode = self.process.wait()
        self.out_fh.close()
        self.log_fh.close()
        if returncode != 0 and exc[0] is None:
            raise subprocess.CalledProcessError(returncode, "pigz")


def read_pairs(bam_file):
    """
    Sorts a BAM file by read name and yields (read 1, read 2) of each pair.
    Reads whose mate is not next to them are skipped.
    """
    command = [
        "samtools", "sort", "-@", threads, "-n", "-u",
        "-T", f"{output}.nmsrt", "-o", "-", bam_file,
    ]  # fmt: skip
    logging.info(f"Running: {' '.join(command)}")
    unpaired = 0
    with open(log, "a") as log_fh, subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=log_fh
    ) as sort, pysam.AlignmentFile(sort.stdout, "rb") as fh:
        previous = None
        for read in fh:
            if previous is not None and read.query_name == previous.query_name:
                yield (previous, read) if previous.is_read1 else (read, previous)
                previous = None
            else:
                if previous is not None:
                    unpaired += 1
                previous = read
        if previous is not None:
            unpaired += 1
    if sort.returncode != 0:
        raise subprocess.CalledProcessError(sort.returncode, command)
    if unpaired:
        logging.warning(f"Skipped {unpaired} reads without a mate next to them")


def tn5_shift(start, end, strand):
    """Shifts a read for the Tn5 insertion (+4 bp plus, -5 bp minus strand)"""
    if strand == "+":
        start += 4
        if start >= end:
            start = end - 1
    else:
        end -= 5
        if start >= end:
            end = start + 1
    return start, end


def tagalign_lines(pairs, shift):
    """Yields the two tagAlign lines (bytes) of each read pair"""
    for pair in pairs:
        lines = []
        for read in pair:
            start, end = read.reference_start, read.reference_end
            strand = "-" if read.is_reverse else "+"
            if shift:
                start, end = tn5_shift(start, end, strand)
            lines.append(f"{read.reference_name}\t{start}\t{end}\tN\t1000\t{strand}\n")
        yield "".join(lines).encode()


try:
    if subsample == 0:
        # Convert to tagAlign and shift in one pass
        # -----------------------------------------------------
        n_pairs = 0
        with PigzWriter(output) as fh:
            for pair_lines in tagalign_lines(read_pairs(bam), shift=True):
                fh.write(pair_lines)
                n_pairs += 1
        logging.info(f"Wrote {n_pairs} read pairs")
    else:
        # Convert to tagAlign, subsample, then shift
        # -----------------------------------------------------
        n_pairs = 0
        n_bytes = 0
        with PigzWriter(unshifted) as fh:
            for pair_lines in tagalign_lines(read_pairs(bam), shift=False):
                fh.write(pair_lines)
                n_pairs += 1
                n_bytes += len(pair_lines)
        logging.info(f"Converted {n_pairs} read pairs")

        n_keep = min(subsample // 2, n_pairs)
        logging.info(f"Subsampling {n_keep} read pairs (random seed: {n_bytes})")
        keep = bytearray([1]) * n_keep + bytearray(n_pairs - n_keep)
        random.Random(n_bytes).shuffle(keep)

        with gzip.open(unshifted, "rb") as fh_in, PigzWriter(output) as fh_out:
            for selected in keep:
                mates = [next(fh_in), next(fh_in)]
                if not selected:
                    continue
                for mate in mates:
                    chrom, start, end, name, score, strand = mate.rstrip(b"\n").split(
                        b"\t"
                    )
                    start, end = tn5_shift(int(start), int(end), strand.decode())
                    fh_out.write(
                        b"\t".join(
                            [chrom, b"%d" % start, b"%d" % end, name, score, strand]
                        )
                        + b"\n"
                    )
        os.remove(unshifted)

    if n_pairs == 0:
        raise ValueError(f"Empty tagAlign: {output}")
except Exception:
    logging.exception("Conversion to tagAlign failed")
    for tmp in [unshifted, output]:
        if os.path.exists(tmp):
            os.remove(tmp)
    raise

logging.info("Done!")
