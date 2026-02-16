# import basic packages
import re
import pandas as pd
from snakemake.utils import validate

# validate sample sheet and config file
#validate(samples, schema="../schemas/samples.schema.yaml")
#validate(config, schema="../schemas/config.schema.yaml")


def samples():
    """
    Checks sample names/files and returns sample wildcard values for Snakemake.
    Paired-end data assumed.
    """
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"].tolist()  # this gets just IP samples
    
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


def ataqv_organism():
    genome = config["genome"]["ensembl"]
    if re.match("hg", genome):
        return "human"
    elif re.match("mm", genome):
        return "mouse"
    else:
        raise ValueError(f"Unsupported genome: {genome}")