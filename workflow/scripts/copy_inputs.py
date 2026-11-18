"""
Copies external input files of the config (e.g. te: correlation: datasets),
unchanged, into a directory of the results, so that the results keep the
input they were made from. Files keep their names (numbered if names repeat
or a file is named datasets.tsv); datasets.tsv lists the name in the config,
the source and the copy.
"""

import logging
import os
import shutil

files = snakemake.params["files"]
out_dir = snakemake.output[0]

logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    handlers=[logging.FileHandler(snakemake.log[0])],
    force=True,
)

os.makedirs(out_dir, exist_ok=True)
names = [os.path.basename(path) for path in files.values()]
rows = ["name\tsource\tcopy"]
for i, ((name, path), base) in enumerate(zip(files.items(), names), start=1):
    # Numbered if names repeat or the name is that of the manifest
    copy = f"{i}_{base}" if names.count(base) > 1 or base == "datasets.tsv" else base
    shutil.copy2(path, os.path.join(out_dir, copy))
    rows.append(f"{name}\t{os.path.abspath(path)}\t{copy}")
    logging.info(f"Copied {path} ({name}) to {out_dir}/{copy}")
with open(os.path.join(out_dir, "datasets.tsv"), "w") as fh:
    fh.write("\n".join(rows) + "\n")
