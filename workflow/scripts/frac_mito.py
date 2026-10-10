"""
Fraction of mitochondrial reads, as in the ENCODE ATAC-seq pipeline
(encode_task_frac_mito.py):

    frac_mito_reads = Rm / (Rn + Rm)

ENCODE counts mapped reads with SAMstats (each mate separately):
- Rm: reads that align to the mitochondrial genome, from a separate alignment
  of all reads to the mitochondrial chromosome only
- Rn: reads with an alignment outside the mitochondrial chromosome, from the
  raw alignment (bowtie2 -k) after removing the mitochondrial alignments

Here both are counted on the raw alignment (bowtie2 -k), which contains all
alignments of reads with up to multimapping + 1 alignments:
- Rm: reads with any (primary or secondary) alignment on the mitochondrial
  chromosome. This includes reads from nuclear copies of mitochondrial DNA
  (NUMTs) that also align to the mitochondrial genome, as ENCODE's
  mitochondria-only alignment does.
- Rn: mapped reads minus the reads whose alignments are all mitochondrial

Read names are collected with samtools and sort, which keeps memory use low.
"""

import logging
import os
import subprocess
import tempfile

# Load Snakemake variables
bam = snakemake.input["bam"]
output = snakemake.output[0]
mito = snakemake.params["mito"]
threads = snakemake.threads
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def run(command):
    """Runs a shell command and returns its output"""
    logging.info(f"Running: {command}")
    with open(log, "a") as fh:
        return subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=fh,
            text=True,
        ).stdout


chroms = run(f"samtools idxstats {bam} | cut -f 1").split()
# Mapped reads: each mapped read has one primary alignment
mapped = int(run(f"samtools view -c -@ {threads} -F 0x904 {bam}"))

mito_reads = only_mito_reads = 0
if mito in chroms:
    tmp_dir = tempfile.mkdtemp(dir=os.path.dirname(output))
    sort = f"sort -u -S 1G -T {tmp_dir}"
    # Read 1 and read 2 are counted separately
    for flag in ["0x40", "0x80"]:
        names = os.path.join(tmp_dir, f"mito_{flag}.txt")
        run(
            f"samtools view -F 0x4 -f {flag} {bam} {mito} | cut -f 1 | {sort} > {names}"
        )
        n_mito = int(run(f"wc -l < {names}"))
        # Of these reads, those with an alignment outside the mitochondrial
        # chromosome
        n_also_nuclear = int(
            run(
                f"samtools view -@ {threads} -F 0x4 -f {flag} -N {names} {bam} | "
                f"awk -v m={mito} '$3 != m' | cut -f 1 | {sort} | wc -l"
            )
        )
        mito_reads += n_mito
        only_mito_reads += n_mito - n_also_nuclear
        os.remove(names)
    os.rmdir(tmp_dir)

non_mito_reads = mapped - only_mito_reads
total = non_mito_reads + mito_reads
frac = mito_reads / total if total > 0 else 0
logging.info(
    f"mapped reads: {mapped}, mito reads: {mito_reads}, "
    f"of which only mitochondrial: {only_mito_reads}, non-mito reads: {non_mito_reads}"
)

with open(output, "w") as fh:
    fh.write("non_mito_reads\tmito_reads\tfrac_mito_reads\n")
    fh.write(f"{non_mito_reads}\t{mito_reads}\t{frac:.6f}\n")
logging.info("Done!")
