# import basic packages
import os
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


def true_replicate_pairs(condition):
    """Returns the pairs of true replicates of a condition (rep1_vs_rep2, ...)"""
    return [p for p in peak_pairs(condition) if re.match(r"rep\d+_vs_", p)]


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


# ENCODE annotated regions (universal DHS, promoters, enhancers) used for the
# fraction of reads in annotated regions; only available for GRCh38 and mm10
ANNOTATION_REGIONS = {
    "hg38": {"dnase": "ENCFF304XEX", "prom": "ENCFF140XLU", "enh": "ENCFF212UAV"},
    "test": {"dnase": "ENCFF304XEX", "prom": "ENCFF140XLU", "enh": "ENCFF212UAV"},
    "mm38": {"dnase": "ENCFF015KVI", "prom": "ENCFF206BQS", "enh": "ENCFF580RGZ"},
}


def annotation_regions():
    """Annotated region files for the genome (empty if not available)"""
    regions = ANNOTATION_REGIONS.get(config["genome"]["ensembl"], {})
    return {region: f"resources/annotation_regions/{region}.bed" for region in regions}


# ENCODE exclusion lists (blacklists) used by the ENCODE ATAC-seq pipeline
ENCODE_BLACKLISTS = {
    "hg38": "ENCFF356LFX",
    "test": "ENCFF356LFX",
    "mm38": "ENCFF547MET",
}

# Genomes without an ENCODE exclusion list get the list of the previous
# assembly lifted over with a UCSC chain: (ENCODE accession, UCSC chain)
LIFTED_BLACKLISTS = {
    "mm39": ("ENCFF547MET", "mm10/liftOver/mm10ToMm39.over.chain.gz"),
}

# ENCODE TSS files (one TSS per protein coding gene) used by the ENCODE
# ATAC-seq pipeline for the TSS enrichment: (ENCODE accession, UCSC chain to
# lift it to the genome, or None)
ENCODE_TSS = {
    "hg38": ("ENCFF493CCB", None),
    "test": ("ENCFF493CCB", None),
    "mm38": ("ENCFF498BEJ", None),
    "mm39": ("ENCFF498BEJ", "mm10/liftOver/mm10ToMm39.over.chain.gz"),
}


def deseq2_enabled():
    """The DESeq2 analysis of all samples is optional and needs at least 2 samples"""
    return config.get("deseq2", {}).get("enable", False) and len(SAMPLES) > 1


def deseq2_peak_method():
    """Reproducibility method of the optimal peak sets used as consensus peaks"""
    method = config.get("deseq2", {}).get("peak_method", "overlap")
    if deseq2_enabled() and method not in reproducibility_methods():
        raise ValueError(f"deseq2: peak_method {method} needs idr: enable: true")
    return method


def reference_condition():
    """
    Returns the reference (control) condition, marked with yes in the
    reference column of the sample sheet, or None without that column
    """
    if "reference" not in SAMPLE_SHEET.columns:
        return None
    marked = SAMPLE_SHEET.groupby("condition", sort=False)["reference"].unique()
    mixed = [c for c, values in marked.items() if len(values) > 1]
    if mixed:
        raise ValueError(
            f"All samples of a condition need the same reference value: {mixed}"
        )
    reference = [c for c, values in marked.items() if values[0] == "yes"]
    if len(reference) != 1:
        raise ValueError(
            "Exactly one condition needs reference yes in the sample sheet, "
            f"found: {reference}"
        )
    return reference[0]


def deseq2_levels():
    """
    Returns the conditions in the order of the DESeq2 factor levels: the
    reference condition first, then the sample sheet order. The reference
    column is only checked when DESeq2 is run.
    """
    if not deseq2_enabled():
        return CONDITIONS
    reference = reference_condition()
    if reference is None:
        return CONDITIONS
    return [reference] + [c for c in CONDITIONS if c != reference]


def deseq2_contrasts():
    """
    Returns the comparisons of the differential accessibility analysis as
    {name: (treatment, control)}: deseq2: contrasts, or otherwise each
    condition with at least 2 replicates against the reference condition of
    the sample sheet.
    """
    configured = config.get("deseq2", {}).get("contrasts") or []
    for pair in configured:
        unknown = [c for c in pair if c not in CONDITIONS]
        if unknown:
            raise ValueError(f"Unknown conditions in deseq2: contrasts: {unknown}")
        if pair[0] == pair[1]:
            raise ValueError(f"Contrast {pair} must compare different conditions")
        if min(len(replicates(c)) for c in pair) < 2:
            raise ValueError(f"Contrast {pair} needs at least 2 replicates of each")
    if not configured:
        replicated = [c for c in CONDITIONS if len(replicates(c)) > 1]
        if len(replicated) < 2:
            return {}
        reference = reference_condition()
        if reference is None:
            raise ValueError(
                "Mark the reference (control) condition with a reference column "
                "(yes/no) in the sample sheet, or set deseq2: contrasts"
            )
        if reference not in replicated:
            raise ValueError(
                f"The reference condition {reference} needs at least 2 replicates"
            )
        configured = [(c, reference) for c in replicated if c != reference]
    contrasts = {}
    for t, c in configured:
        name = f"{t}_vs_{c}"
        if contrasts.setdefault(name, (t, c)) != (t, c):
            raise ValueError(
                f"Contrasts {contrasts[name]} and {(t, c)} have the same name {name}"
            )
    return contrasts


class GenomeFiles:
    """Gets URLs and file names of fasta and GTF files for a given genome and build"""

    # create genome directory
    os.makedirs("resources/", exist_ok=True)

    def __init__(self, genome, build):
        self.genome = genome
        self.build = build

        # base URLs
        base_url = f"https://ftp.ensembl.org/pub/release-{build}/"

        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
            elif genome == "hg38":
                name = "GRCh38"

            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            )

        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"

            elif genome == "mm39":
                name = "GRCm39"

            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            )

        elif "dm" in genome:
            if genome == "dm6":
                name = "BDGP6.46"

            self.fasta_url = f"{base_url}fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.{name}.dna.toplevel.fa.gz"
            self.gtf_url = f"{base_url}gtf/drosophila_melanogaster/Drosophila_melanogaster.{name}.{build}.gtf.gz"

        elif "test" in genome:
            self.fasta_url = "https://github.com/niekwit/damid-seq/raw/main/.test_pe/Homo_sapiens.GRCh38.dna.primary_assembly_chr11.fa.gz"
            self.gtf_url = "https://github.com/niekwit/atac-seq/raw/main/.test/Homo_sapiens.GRCh38.115_chr11.gtf.gz"

        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)

    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file"""

        return f"resources/{os.path.basename(url).replace('.gz', '')}"


# UCSC names of the genomes, for the RepeatMasker annotation of the TE module
UCSC_GENOMES = {
    "hg19": "hg19",
    "hg38": "hg38",
    "test": "hg38",
    "mm38": "mm10",
    "mm39": "mm39",
    "dm6": "dm6",
}


def te_enabled():
    """The transposable element (TE) module is optional (off by default)"""
    return config.get("te", {}).get("enable", False)
