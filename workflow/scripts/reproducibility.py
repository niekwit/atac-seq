"""
Selects the optimal and conservative peak sets of a condition and calculates
the rescue and self-consistency ratios, as in the ENCODE ATAC-seq pipeline
(encode_task_reproducibility.py).
"""

import logging
import shutil

# Load Snakemake variables
peaks = list(snakemake.input["peaks"])
peaks_pr = list(snakemake.input["peaks_pr"])
peak_ppr = snakemake.input["peak_ppr"]
pairs = snakemake.params["pairs"]
optimal_out = snakemake.output["optimal"]
conservative_out = snakemake.output["conservative"]
qc_out = snakemake.output["qc"]
log = snakemake.log[0]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
    force=True,
)


def num_lines(f):
    with open(f) as fh:
        return sum(1 for _ in fh)


def ratio(a, b):
    """max/min ratio; infinite when one of the sets is empty"""
    lo, hi = min(a, b), max(a, b)
    if lo == 0:
        return 0.0 if hi == 0 else float("inf")
    return hi / lo


# N: number of peaks in the self-pseudoreplicate peak sets
N = [num_lines(p) for p in peaks_pr]
num_rep = len(peaks_pr)
pr_pairs = [f"rep{i + 1}-pr1_vs_rep{i + 1}-pr2" for i in range(num_rep)]

if len(peaks) > 0:
    # Multiple replicates
    true_pairs = [p for p in pairs if p not in pr_pairs and not p.startswith("pooled")]
    num_peaks_tr = [num_lines(p) for p in peaks]
    Nt = max(num_peaks_tr)
    Np = num_lines(peak_ppr[0])
    rescue_ratio = ratio(Np, Nt)
    self_consistency_ratio = ratio(max(N), min(N))

    Nt_idx = num_peaks_tr.index(Nt)
    conservative_set = true_pairs[Nt_idx]
    conservative_peak = peaks[Nt_idx]
    N_conservative = Nt
    if Nt > Np:
        optimal_set = conservative_set
        optimal_peak = conservative_peak
        N_optimal = N_conservative
    else:
        optimal_set = "pooled-pr1_vs_pooled-pr2"
        optimal_peak = peak_ppr[0]
        N_optimal = Np
else:
    # Single replicate
    Nt = 0
    Np = 0
    rescue_ratio = 0.0
    self_consistency_ratio = 1.0

    conservative_set = "rep1-pr1_vs_rep1-pr2"
    conservative_peak = peaks_pr[0]
    N_conservative = N[0]
    optimal_set = conservative_set
    optimal_peak = conservative_peak
    N_optimal = N_conservative

reproducibility = "pass"
if rescue_ratio > 2.0 or self_consistency_ratio > 2.0:
    reproducibility = "borderline"
if rescue_ratio > 2.0 and self_consistency_ratio > 2.0:
    reproducibility = "fail"

logging.info(f"Optimal set: {optimal_set} ({N_optimal} peaks)")
logging.info(f"Conservative set: {conservative_set} ({N_conservative} peaks)")
logging.info(f"Reproducibility: {reproducibility}")
shutil.copyfile(optimal_peak, optimal_out)
shutil.copyfile(conservative_peak, conservative_out)

header = (
    ["Nt"]
    + [f"N{i + 1}" for i in range(num_rep)]
    + [
        "Np",
        "N_opt",
        "N_consv",
        "opt_set",
        "consv_set",
        "rescue_ratio",
        "self_consistency_ratio",
        "reproducibility",
    ]
)
values = (
    [Nt]
    + N
    + [
        Np,
        N_optimal,
        N_conservative,
        optimal_set,
        conservative_set,
        rescue_ratio,
        self_consistency_ratio,
        reproducibility,
    ]
)
with open(qc_out, "w") as fh:
    fh.write("\t".join(header) + "\n")
    fh.write("\t".join(str(v) for v in values) + "\n")
logging.info("Done!")
