"""
Collects the ENCODE QC metrics of all samples and conditions and compares
them with the ENCODE ATAC-seq data standards
(https://www.encodeproject.org/atac-seq/).
"""

import gzip
import json
import logging
import re

# Load Snakemake variables
samples = snakemake.params["samples"]
conditions = snakemake.params["conditions"]
replicates = snakemake.params["replicates"]
methods = snakemake.params["methods"]
samples_out = snakemake.output["samples"]
conditions_out = snakemake.output["conditions"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def read_table(f):
    """Returns the first data row of a TSV file with header as a dict"""
    with open(f) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        values = fh.readline().rstrip("\n").split("\t")
    return dict(zip(header, values))


def read_value(f):
    with open(f) as fh:
        return float(fh.read().strip())


def num_lines(f):
    with open(f) as fh:
        return sum(1 for _ in fh)


def alignment_rate(f):
    with open(f) as fh:
        m = re.search(r"([\d.]+)% overall alignment rate", fh.read())
    return float(m.group(1)) / 100 if m else float("nan")


def mapped_reads(f):
    """Reads mapped from samtools stats output"""
    with open(f) as fh:
        for line in fh:
            if line.startswith("SN\treads mapped:"):
                return int(line.split("\t")[2])
    return 0


def find_key(obj, key):
    """Returns the first value of key in a nested JSON object"""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        obj = list(obj.values())
    if isinstance(obj, list):
        for item in obj:
            value = find_key(item, key)
            if value is not None:
                return value
    return None


def ataqv_tss_enrichment(f):
    with gzip.open(f, "rt") as fh:
        value = find_key(json.load(fh), "tss_enrichment")
    return float(value) if value is not None else float("nan")


def grade(value, ideal, acceptable=None):
    """ENCODE standard: ideal / acceptable / concerning"""
    if value != value:  # NaN
        return "NA"
    if value >= ideal:
        return "ideal"
    if acceptable is not None and value >= acceptable:
        return "acceptable"
    return "concerning"


def bottlenecking(pbc1):
    """PCR bottlenecking category of PBC1, as in the ENCODE QC report"""
    if pbc1 >= 0.9:
        return "none"
    if pbc1 >= 0.8:
        return "mild"
    if pbc1 >= 0.5:
        return "moderate"
    return "severe"


# ENCODE TSS enrichment thresholds (ideal, acceptable) depend on the genome
# (https://www.encodeproject.org/atac-seq/#standards)
TSS_THRESHOLDS = {
    "hg19": (10, 6),
    "hg38": (7, 5),
    "test": (7, 5),
    "mm38": (15, 10),
    "mm39": (15, 10),
}
tss_ideal, tss_acceptable = TSS_THRESHOLDS.get(
    snakemake.params["genome"], (float("nan"), float("nan"))
)


def by_sample(name, sample):
    return [f for f in snakemake.input[name] if re.search(rf"/{sample}\.", f)][0]


# Per sample metrics
# Standards: https://www.encodeproject.org/atac-seq/#standards
sample_rows = []
for sample in samples:
    logging.info(f"Collecting QC metrics of {sample}")
    rate = alignment_rate(by_sample("bowtie2", sample))
    mito = read_table(by_sample("frac_mito", sample))
    pbc = read_table(by_sample("lib_complexity", sample))
    nodup = mapped_reads(by_sample("nodup_stats", sample))
    frip = read_value(by_sample("frip", sample))
    tss = read_value(by_sample("tss_enrich", sample))
    tss_ataqv = ataqv_tss_enrichment(by_sample("ataqv", sample))
    nrf, pbc1, pbc2 = float(pbc["NRF"]), float(pbc["PBC1"]), float(pbc["PBC2"])

    sample_rows.append(
        {
            "sample": sample,
            "alignment_rate": f"{rate:.4f}",
            "alignment_rate_standard": grade(rate, 0.95, 0.80),
            "frac_mito_reads": mito["frac_mito_reads"],
            "nodup_non_mito_reads": nodup,
            "nodup_non_mito_reads_standard": grade(nodup, 50e6),
            # Library complexity: ENCODE prefers NRF > 0.9, PBC1 > 0.9 and
            # PBC2 > 3; its QC report accepts NRF > 0.8 and calls
            # PBC1 0.8-0.9 mild bottlenecking
            "NRF": f"{nrf:.4f}",
            "NRF_standard": grade(nrf, 0.9, 0.8),
            "PBC1": f"{pbc1:.4f}",
            "PBC1_standard": grade(pbc1, 0.9, 0.8),
            "PCR_bottlenecking": bottlenecking(pbc1),
            "PBC2": f"{pbc2:.4f}",
            "PBC2_standard": grade(pbc2, 3, 1),
            "num_peaks": num_lines(by_sample("peaks", sample)),
            "FRiP": f"{frip:.4f}",
            "FRiP_standard": grade(frip, 0.3, 0.2),
            "TSS_enrichment": f"{tss:.2f}",
            "TSS_enrichment_standard": grade(tss, tss_ideal, tss_acceptable),
            # ataqv calculates TSS enrichment on a different scale, so it is
            # not graded with the ENCODE thresholds
            "TSS_enrichment_ataqv": f"{tss_ataqv:.2f}",
        }
    )

with open(samples_out, "w") as fh:
    fh.write("\t".join(sample_rows[0].keys()) + "\n")
    for row in sample_rows:
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

# Per condition reproducibility metrics
condition_rows = []
for method in methods:
    for condition in conditions:
        qc_file = f"results/{method}/{condition}/reproducibility.qc"
        qc = read_table(qc_file)
        frip_file = f"results/{method}/{condition}/{qc['opt_set']}.frip.qc"
        n_opt = int(qc["N_opt"])
        condition_rows.append(
            {
                "method": method,
                "condition": condition,
                "replicates": len(replicates[condition]),
                "replicates_standard": (
                    "ideal" if len(replicates[condition]) >= 2 else "concerning"
                ),
                "optimal_set": qc["opt_set"],
                "num_optimal_peaks": n_opt,
                "num_optimal_peaks_standard": grade(n_opt, 150000, 100000),
                "FRiP_optimal": f"{read_value(frip_file):.4f}",
                "FRiP_optimal_standard": grade(read_value(frip_file), 0.3, 0.2),
                "conservative_set": qc["consv_set"],
                "num_conservative_peaks": qc["N_consv"],
                "rescue_ratio": qc["rescue_ratio"],
                "self_consistency_ratio": qc["self_consistency_ratio"],
                "reproducibility": qc["reproducibility"],
            }
        )

with open(conditions_out, "w") as fh:
    fh.write("\t".join(condition_rows[0].keys()) + "\n")
    for row in condition_rows:
        fh.write("\t".join(str(v) for v in row.values()) + "\n")

logging.info("Done!")
