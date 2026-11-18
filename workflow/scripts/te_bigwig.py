"""
Signal track of the Tn5 insertions (5' end of each read, --Offset 1) of a
sample of the TE module, with bamCoverage: per million kept insertions
(--normalizeUsing CPM) or, with the size factors of the TE consensus peaks,
scaled with --scaleFactor 1e6 / (size factor x geometric mean of the library
sizes) (te: bigwig_normalisation: deseq2).
"""

import csv
import subprocess

sample = snakemake.wildcards["sample"]
log = snakemake.log[0]

size_factors = snakemake.input.get("size_factors", [])
if size_factors:
    with open(size_factors[0]) as fh:
        rows = {row["sample"]: row for row in csv.DictReader(fh, delimiter="\t")}
    scale = float(rows[sample]["scale_factor"])
    normalisation = ["--scaleFactor", f"{scale:.10g}"]
else:
    normalisation = ["--normalizeUsing", "CPM"]

command = [
    "bamCoverage",
    "--binSize",
    str(snakemake.params["bin_size"]),
    *normalisation,
    "--Offset",
    "1",
    *snakemake.params["extra"].split(),
    "-p",
    str(snakemake.threads),
    "-b",
    snakemake.input["bam"],
    "-o",
    snakemake.output[0],
]
with open(log, "w") as fh:
    fh.write(" ".join(command) + "\n")
    fh.flush()
    subprocess.run(command, check=True, stdout=fh, stderr=fh)
