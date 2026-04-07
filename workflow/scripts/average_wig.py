import logging
from snakemake.shell import shell
import pandas as pd

# Load Snakemake variables
log = snakemake.log[0]
log_fmt = snakemake.log_fmt_shell(stdout=True, stderr=True)

all_bw = snakemake.input
condition = snakemake.wildcards["condition"]
wig = snakemake.output["wig"]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)
logging.info(f"Condition processed: {condition}")
logging.info(f"All bigwig files: {all_bw}")

# Load sample info and get list of samples for the current condition
csv = pd.read_csv("config/samples.csv")
samples = csv[csv["condition"] == condition]["sample"].tolist()

logging.info(f"Samples in condition '{condition}': {samples}")

# Subset bigwigs to only those from the samples in the current condition
bw = [x for x in all_bw if any(sample in x for sample in samples)]
logging.info(f"Bigwig files for condition '{condition}': {bw}")

# Create average bigwig file
command = "wiggletools write {wig} mean {bw} {log_fmt}"
logging.debug(command)
shell(command)
logging.info("Done!")
