"""
deepTools computeMatrix scale-regions matrix of the TE signal tracks over the
full-length copies of the te: profile: subfamilies entries
(resources/te/full_length_copies.bed). The copies of an entry are followed by
"#<entry> (n=...)", which deepTools 4 no longer reads as a group label, so the
copies are split into one BED file per entry, named after the label (deepTools
group label: the file name).
"""

import logging
import os
import shlex
import subprocess
import tempfile

bed_file = snakemake.input["bed"]
bigwigs = list(snakemake.input["bigwigs"])
out_matrix = snakemake.output[0]
log = snakemake.log[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(log, mode="w")],
    force=True,
)

with tempfile.TemporaryDirectory(dir=os.path.dirname(out_matrix)) as tmp:
    beds, pending = [], []
    for line in open(bed_file):
        if line.startswith("#"):
            label = line[1:].split("\t")[0].strip()
            path = os.path.join(tmp, f"{label}.bed")
            with open(path, "w") as fh:
                fh.writelines(pending)
            beds.append(path)
            pending = []
        elif line.strip():
            pending.append(line)
    if pending:
        raise ValueError(f"{bed_file}: copies after the last group label")

    command = [
        "computeMatrix",
        "scale-regions",
        "-R",
        *beds,
        "-S",
        *bigwigs,
        "--samplesLabel",
        *snakemake.params["labels"],
        "--regionBodyLength",
        str(snakemake.params["body_length"]),
        "--beforeRegionStartLength",
        str(snakemake.params["flank"]),
        "--afterRegionStartLength",
        str(snakemake.params["flank"]),
        "--binSize",
        "50",
        "--missingDataAsZero",
        *shlex.split(snakemake.params["extra"]),
        "-p",
        str(snakemake.threads),
        "-o",
        out_matrix,
    ]
    logging.info(shlex.join(command))
    with open(log, "a") as fh:
        subprocess.run(command, check=True, stdout=fh, stderr=fh)
