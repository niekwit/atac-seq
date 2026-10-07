"""
Removes reads with more than k alignments from a query name sorted SAM
stream (stdin to stdout). The remaining reads keep all their alignments;
secondary alignments are removed afterwards with samtools.

Based on assign_multimappers.py of the ENCODE ATAC-seq pipeline
(https://github.com/ENCODE-DCC/atac-seq-pipeline, MIT License,
Copyright (c) 2017 ENCODE DCC). Unlike the original, the alignments of
the last read in the stream are also written.
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Keeps reads with at most k alignments and discards all others"
    )
    parser.add_argument("-k", type=int, required=True, help="Alignment number cutoff")
    parser.add_argument("--paired-end", action="store_true", help="Data is paired-end")
    return parser.parse_args()


def write_reads(reads, cutoff):
    if 0 < len(reads) <= cutoff:
        sys.stdout.writelines(reads)


def main():
    args = parse_args()
    cutoff = args.k * 2 if args.paired_end else args.k

    current_reads = []
    current_qname = None

    for line in sys.stdin:
        if line.startswith("@"):
            sys.stdout.write(line)
            continue

        qname = line.split("\t", 1)[0]
        if qname != current_qname:
            write_reads(current_reads, cutoff)
            current_reads = []
            current_qname = qname
        current_reads.append(line)

    write_reads(current_reads, cutoff)


if __name__ == "__main__":
    main()
