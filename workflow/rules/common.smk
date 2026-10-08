# import basic packages
import re
import pandas as pd
from snakemake.utils import validate

# validate config file
validate(config, schema="../schemas/config.schema.yaml")

SAMPLE_SHEET = pd.read_csv("config/samples.csv")
validate(SAMPLE_SHEET, schema="../schemas/samples.schema.yaml")


def samples():
    """
    Checks sample names/files and returns sample wildcard values for Snakemake.
    Paired-end data assumed.
    """
    SAMPLES = SAMPLE_SHEET["sample"].tolist()

    # Check if sample names contain any characters that are not alphanumeric or underscore
    illegal = []
    for sample in SAMPLES:
        if not re.match("^[a-zA-Z0-9_]*$", sample):
            illegal.append(sample)
    if len(illegal) != 0:
        illegal = "\n".join(illegal)
        raise ValueError(f"Following samples contain illegal characters:\n{illegal}")

    # Check if each sample name ends with _[0-9]
    wrong = []
    for sample in SAMPLES:
        if not re.match(".*_[0-9]$", sample):
            wrong.append(sample)
    if len(wrong) != 0:
        wrong = "\n".join(wrong)
        raise ValueError(f"Following samples do not end with _[0-9]:\n{wrong}")

    # Check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        r1 = f"reads/{sample}_R1_001.fastq.gz"
        r2 = f"reads/{sample}_R2_001.fastq.gz"
        if not os.path.isfile(r1):
            not_found.append(r1)
        if not os.path.isfile(r2):
            not_found.append(r2)

    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"Following files not found:\n{not_found}")

    return SAMPLES


def conditions():
    conditions = SAMPLE_SHEET["condition"].unique().tolist()

    illegal = [c for c in conditions if not re.match("^[a-zA-Z0-9_]+$", c)]
    if len(illegal) != 0:
        illegal = "\n".join(illegal)
        raise ValueError(f"Following conditions contain illegal characters:\n{illegal}")

    return conditions


def replicates(condition):
    """
    Returns the samples (biological replicates) of a condition,
    in sample sheet order. ENCODE's rep1, rep2, ... refer to this order.
    """
    return SAMPLE_SHEET[SAMPLE_SHEET["condition"] == condition]["sample"].tolist()


def peak_pairs(condition):
    """
    Returns the ENCODE names of the peak set comparisons of a condition:
    all pairs of true replicates, the self-pseudoreplicates of each
    replicate and, with more than one replicate, the pooled pseudoreplicates.
    """
    n = len(replicates(condition))
    reps = range(1, n + 1)
    true_reps = [f"rep{i}_vs_rep{j}" for i in reps for j in reps if i < j]
    pseudo_reps = [f"rep{i}-pr1_vs_rep{i}-pr2" for i in reps]
    pooled_pseudo_reps = ["pooled-pr1_vs_pooled-pr2"] if n > 1 else []
    return true_reps + pseudo_reps + pooled_pseudo_reps


def pair_prefixes(condition, pair):
    """
    Returns the tagAlign/peak file prefixes (relative to results/tagalign/
    and results/macs2/) of the two peak sets and the pooled peak set of a
    comparison, and the prefix of the tagAlign used for FRiP.
    """
    reps = replicates(condition)
    pooled = f"pooled/{condition}"

    m = re.fullmatch(r"rep(\d+)_vs_rep(\d+)", pair)
    if m:
        rep1, rep2 = reps[int(m.group(1)) - 1], reps[int(m.group(2)) - 1]
        return {"peak1": rep1, "peak2": rep2, "pooled": pooled, "ta": pooled}

    m = re.fullmatch(r"rep(\d+)-pr1_vs_rep\1-pr2", pair)
    if m:
        rep = reps[int(m.group(1)) - 1]
        return {"peak1": f"{rep}.pr1", "peak2": f"{rep}.pr2", "pooled": rep, "ta": rep}

    if pair == "pooled-pr1_vs_pooled-pr2":
        return {
            "peak1": f"{pooled}.pr1",
            "peak2": f"{pooled}.pr2",
            "pooled": pooled,
            "ta": pooled,
        }

    raise ValueError(f"Unknown peak comparison: {pair}")


def pair_peaks(wildcards):
    prefixes = pair_prefixes(wildcards.condition, wildcards.pair)
    return {
        key: f"results/macs2/{prefixes[key]}.narrowPeak"
        for key in ["peak1", "peak2", "pooled"]
    }


def frip_tagalign(wildcards):
    """
    Returns the tagAlign used to calculate the FRiP of a peak file.
    """
    m = re.fullmatch(r"results/macs2/(.+)", wildcards.stem)
    if m:
        return f"results/tagalign/{m.group(1)}.tagAlign.gz"

    m = re.fullmatch(r"results/(idr|overlap)/([^/]+)/([^/]+)", wildcards.stem)
    if m:
        ta = pair_prefixes(m.group(2), m.group(3))["ta"]
        return f"results/tagalign/{ta}.tagAlign.gz"

    raise ValueError(f"No tagAlign for peak file: {wildcards.stem}")


def reproducibility_methods():
    return ["overlap", "idr"] if config["idr"]["enable"] else ["overlap"]


def ataqv_organism():
    genome = config["genome"]["ensembl"]
    if re.match("hg", genome) or genome == "test":
        return "human"
    elif re.match("mm", genome):
        return "mouse"
    else:
        raise ValueError(f"Unsupported genome: {genome}")


def runtime(minutes):
    """
    Runtime limit (minutes) for cluster executors, e.g. SLURM's --time.
    The limit is multiplied by the attempt number, so that jobs that run out
    of time get more time when the workflow is run with --retries.
    """
    return lambda wildcards, attempt: minutes * attempt
