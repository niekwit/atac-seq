import logging
import pybedtools as pb
import pandas as pd


# Load Snakemake variables
peaks = snakemake.input["peaks"]
chrom_sizes_file = snakemake.input["chrom_sizes"]
method = snakemake.params["method"]
keep = snakemake.params["keep"]
condition = snakemake.wildcards["condition"]
log = snakemake.log[0]
bed = snakemake.output["bed"]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)

# Load sample info and get list of samples for the current condition
csv = pd.read_csv("config/samples.csv")
samples = csv[csv["condition"] == condition]["sample"].tolist()

# Load chrom_sizes into dictionary
# This is needed to make sure that the extended regions
# do not extend past the beginning/end of the chromosome
chrom_sizes = {}
with open(chrom_sizes_file) as f:
    for line in f:
        (key, value) = line.split()
        chrom_sizes[key] = (0, int(value))

# Subset peak files to only those from the samples in the current condition
peak_files = [x for x in peaks if any(sample in x for sample in samples)]
logging.info(f"Peak files for condition '{condition}': {peak_files}")

# BedTool objects and sort them immediately
logging.info(f"Loading and sorting {len(peak_files)} peak files")
bedtools_objs = [pb.BedTool(f).sort(g=chrom_sizes) for f in peak_files]


# If method is "merge", use bedtools merge to create consensus peaks
if method == "merge":
    logging.info("Using bedtools merge to create consensus peaks")

    # Cat, sort, then merge
    all_peaks = (
        bedtools_objs[0].cat(*bedtools_objs[1:], postmerge=False).sort(g=chrom_sizes)
    )

    # Merge sorted peaks
    consensus_peaks = all_peaks.merge()
else:  # Otherwise, use bedtools multiinter to find consensus peaks
    logging.info("Using bedtools multiinter to create consensus peaks")

    # multi_intersect requires sorted inputs
    raw_intersect = pb.BedTool().multi_intersect(i=bedtools_objs)

    # Filter by the 'keep' parameter (column index 3 is the count)
    consensus_peaks = raw_intersect.filter(lambda x: int(x[3]) >= int(keep))

# Final cleanup to ensure no peaks exceed chromosome boundaries
consensus_peaks = consensus_peaks.truncate_to_chrom(chrom_sizes)

# Consensus peaks only have 3 columns
# Add a 4th column with unique peak IDs
consensus_peaks = consensus_peaks.each(
    lambda x: x.append(f"peak_{x.chrom}_{x.start}_{x.end}") or x
)

logging.info(f"Saving consensus peaks to {bed}")
consensus_peaks.saveas(bed)
logging.info("Done!")
