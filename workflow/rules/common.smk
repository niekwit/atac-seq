# import basic packages
import math
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


def extra(tool):
    """
    Extra command-line arguments of a tool run directly as a command, from
    extra: in the config (empty: none). They come after the workflow's own
    options, so for most tools they override them.
    """
    return config.get("extra", {}).get(tool) or ""


def ataqv_organism():
    genome = config["genome"]["ensembl"]
    if re.match("hg", genome) or genome == "test":
        return "human"
    elif re.match("mm", genome):
        return "mouse"
    else:
        raise ValueError(f"Unsupported genome: {genome}")


def cpus_for_memory(mem_mb, mem_mb_per_cpu=3370):
    """
    CPUs to request (cpus_per_task, e.g. SLURM's --cpus-per-task) for rules
    with fewer threads than CPUs needed for their memory on clusters that
    assign memory per CPU. The rule itself runs with its threads, so these
    jobs do not block cores in local runs.
    """
    return math.ceil(mem_mb / mem_mb_per_cpu)


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


def subsample_min():
    """
    subsample_reads: min subsamples each replicate to the number of reads of
    the smallest library (after deduplication), separately for the ENCODE and
    the TE alignments
    """
    return config["subsample_reads"] == "min"


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


def te_bigwig_normalisation():
    """
    Scaling of the TE signal tracks: cpm (insertions per million kept
    insertions) or deseq2 (DESeq2 size factors of the TE consensus peaks,
    which needs deseq2: enable)
    """
    normalisation = config.get("te", {}).get("bigwig_normalisation", "cpm")
    if normalisation == "deseq2" and te_enabled() and not deseq2_enabled():
        raise ValueError("te: bigwig_normalisation: deseq2 needs deseq2: enable: true")
    return normalisation


def te_profile_config():
    """
    Returns the settings of the profile plot over full-length copies of
    repeat subfamilies (no plot when the list of subfamilies is empty)
    """
    profile = config.get("te", {}).get("profile", {})
    subfamilies = profile.get("subfamilies", [])
    return {
        "subfamilies": [subfamilies] if isinstance(subfamilies, str) else subfamilies,
        "min_consensus_fraction": profile.get("min_consensus_fraction", 0.9),
        "flank": profile.get("flank", 2000),
        "body_length": profile.get("body_length", 6000),
    }


def te_correlation_config():
    """
    Returns the data sets ({name: BED file}) to compare with the differential
    TE peaks with GAT, and the number of GAT samples
    """
    correlation = config.get("te", {}).get("correlation", {})
    datasets = correlation.get("datasets") or {}
    for name, path in datasets.items():
        if not os.path.isfile(path):
            raise ValueError(f"te: correlation: datasets: {name}: no file {path}")
    return {
        "datasets": datasets,
        "gat_samples": correlation.get("gat_samples", 10000),
    }


def rnaseq_config():
    """
    Returns the settings of the integration with RNA-seq (rnaseq: in the
    config), with defaults
    """
    rnaseq = config.get("rnaseq") or {}
    columns = rnaseq.get("columns") or {}
    return {
        "tables": rnaseq.get("tables") or {},
        "columns": {
            "gene_id": columns.get("gene_id") or "",
            "log2fc": columns.get("log2fc") or "log2FoldChange",
            "padj": columns.get("padj") or "padj",
            "base_mean": columns.get("base_mean") or "baseMean",
        },
        "id_type": rnaseq.get("id_type", "auto"),
        "fdr": rnaseq.get("fdr", 0.05),
        "lfc_threshold": rnaseq.get("lfc_threshold", 0),
        "distance_bins": sorted(
            rnaseq.get("distance_bins", [1000, 10000, 50000, 100000])
        ),
        "te_separator": rnaseq.get("te_separator", ":"),
        "pair_lfc_threshold": rnaseq.get("pair_lfc_threshold", 1),
        "tss_profile_flank": rnaseq.get("tss_profile_flank", 3000),
    }


def rnaseq_contrasts():
    """
    Returns the DESeq2 contrasts with RNA-seq DESeq2 results tables
    ({contrast: [tables]}, rnaseq: tables: a path or a list of paths, e.g.
    genes and TEs in separate files); empty without deseq2: enable
    """
    tables = rnaseq_config()["tables"]
    if not tables or not deseq2_enabled():
        return {}
    contrasts = deseq2_contrasts()
    result = {}
    for name, paths in tables.items():
        if name not in contrasts:
            raise ValueError(
                f"rnaseq: tables: {name} is not a DESeq2 contrast "
                f"(contrasts: {', '.join(contrasts) or 'none'})"
            )
        paths = [paths] if isinstance(paths, str) else list(paths)
        if not paths:
            raise ValueError(f"rnaseq: tables: {name}: no table")
        for path in paths:
            if not os.path.isfile(path):
                raise ValueError(f"rnaseq: tables: {name}: no file {path}")
        result[name] = paths
    return result


# Integrations with differential epigenomic data (not part of the ENCODE
# pipeline): the settings of each section of the config ("mark"), the change
# codes for a decrease and an increase in the treatment, and plot labels
# ({treatment} is replaced by the treatment of the contrast)
EPI_MARKS = {
    "h3k9me3": {
        "name": "H3K9me3",
        # Key of the boolean that marks a results table (instead of a BED file)
        "table_key": "diffbind",
        "table_name": "DiffBind",
        "down": "lost",
        "up": "gained",
        "columns": {
            "chr": "seqnames",
            "start": "start",
            "end": "end",
            "effect": "Fold",
            "fdr": "FDR",
        },
        # Config key of the effect column
        "effect_key": "log2fc",
        "fdr": 0.05,
        "min_effect": 0,
        "effect_name": "log2FoldChange",
        "labels": {
            "down": "H3K9me3 lost in {treatment}",
            "up": "H3K9me3 gained in {treatment}",
            "down_short": "H3K9me3 lost",
            "up_short": "H3K9me3 gained",
            "unchanged": "H3K9me3 unchanged",
            "outside": "outside H3K9me3 regions",
            "regions": "H3K9me3 regions",
            "effect": "H3K9me3 log2 fold change",
            "signal": "H3K9me3 signal",
            "sites": "DiffBind sites",
        },
        # Coverage signal: bins without signal are 0
        "missing_as_zero": True,
    },
    "methylation": {
        "name": "DNA methylation",
        "table_key": "methylkit",
        "table_name": "methylKit",
        "down": "hypo",
        "up": "hyper",
        "columns": {
            "chr": "chr",
            "start": "start",
            "end": "end",
            "effect": "meth.diff",
            "fdr": "qvalue",
        },
        "effect_key": "diff",
        "fdr": 0.01,
        "min_effect": 25,
        "effect_name": "meth_diff",
        "labels": {
            "down": "Hypomethylated in {treatment}",
            "up": "Hypermethylated in {treatment}",
            "down_short": "Hypomethylated",
            "up_short": "Hypermethylated",
            "unchanged": "Methylation unchanged",
            "outside": "outside the methylation regions",
            "regions": "methylation regions",
            "effect": "Methylation difference (percentage points)",
            "signal": "CpG methylation (%)",
            "sites": "tested regions",
        },
        # Methylation per CpG: bins without a CpG are missing, not 0%
        "missing_as_zero": False,
    },
}


def epi_config(mark):
    """
    Returns the settings of an integration with differential epigenomic data
    (h3k9me3: or methylation: in the config), with defaults. The effect
    column is under columns: effect (h3k9me3: columns: log2fc, methylation:
    columns: diff in the config)
    """
    defaults = EPI_MARKS[mark]
    section = config.get(mark) or {}
    columns = dict(section.get("columns") or {})
    if defaults["effect_key"] in columns:
        columns["effect"] = columns.pop(defaults["effect_key"])
    return {
        "contrasts": section.get("contrasts") or {},
        "columns": {k: columns.get(k) or v for k, v in defaults["columns"].items()},
        "flip": section.get("flip_fold", section.get("flip_diff", False)),
        "fdr": section.get("fdr", defaults["fdr"]),
        "min_effect": section.get("min_diff", defaults["min_effect"]),
        "profile_flank": section.get("profile_flank", 5000),
        "profile_body_length": section.get("profile_body_length", 5000),
        "min_copies": section.get("min_copies", 20),
    }


def epi_contrasts(mark):
    """
    Returns the DESeq2 contrasts with data of an epigenomic mark
    ({contrast: settings}: regions, table (a results table, not a BED file),
    direction (the change code of a BED file) and bigwigs, None without);
    empty without deseq2: enable
    """
    defaults = EPI_MARKS[mark]
    settings_all = epi_config(mark)
    contrasts = settings_all["contrasts"]
    if not contrasts or not deseq2_enabled():
        return {}
    deseq2 = deseq2_contrasts()
    # deepTools needs the profile lengths to be multiples of the bin size
    bin_size = computematrix_bin_size()
    for setting in ("profile_flank", "profile_body_length"):
        value = settings_all[setting]
        if value % bin_size:
            raise ValueError(
                f"{mark}: {setting} ({value}) must be a multiple of the "
                f"computeMatrix bin size ({bin_size} bp)"
            )
    codes = (defaults["down"], defaults["up"])
    table_key = defaults["table_key"]
    result = {}
    for name, settings in contrasts.items():
        key = f"{mark}: contrasts: {name}"
        if name not in deseq2:
            raise ValueError(
                f"{key} is not a DESeq2 contrast "
                f"(contrasts: {', '.join(deseq2) or 'none'})"
            )
        settings = settings or {}
        regions = settings.get("regions")
        if not regions or not os.path.isfile(regions):
            raise ValueError(f"{key}: regions: no file {regions}")
        table = settings.get(table_key)
        if not isinstance(table, bool):
            raise ValueError(f"{key}: {table_key} must be true or false")
        direction = settings.get("direction", codes[0])
        if direction not in codes:
            raise ValueError(f"{key}: direction must be {' or '.join(codes)}")
        bigwigs = settings.get("bigwigs") or None
        if bigwigs:
            for condition in ("control", "treatment"):
                path = bigwigs.get(condition)
                if not path or not os.path.isfile(path):
                    raise ValueError(f"{key}: bigwigs: {condition}: no file {path}")
        result[name] = {
            "regions": regions,
            "table": table,
            "direction": direction,
            "bigwigs": bigwigs,
        }
    return result


def epi_marks():
    """The epigenomic marks with at least one contrast"""
    return [mark for mark in EPI_MARKS if epi_contrasts(mark)]


def computematrix_bin_size():
    """
    Bin size (bp) of the computeMatrix profiles of the integrations: 50, or
    --binSize/-bs of extra: computematrix (later arguments override earlier)
    """
    found = re.findall(r"(?:--binSize|-bs)[ =](\d+)", extra("computematrix"))
    return int(found[-1]) if found else 50


def epi_changes(mark, contrast):
    """Change codes of a mark in the data of a contrast (decrease first)"""
    settings = epi_contrasts(mark)[contrast]
    defaults = EPI_MARKS[mark]
    if settings["table"]:
        return [defaults["down"], defaults["up"]]
    return [settings["direction"]]


def epi_labels(mark, contrast):
    """Plot labels of a mark, with the treatment of the contrast"""
    treatment = deseq2_contrasts()[contrast][0]
    labels = EPI_MARKS[mark]["labels"]
    return {k: v.format(treatment=treatment) for k, v in labels.items()}


def epi_targets(mark, te_module):
    """
    Output files of the integration with a mark: of the consensus peaks
    (te_module False) or of the TE module (te_module True)
    """
    contrasts = epi_contrasts(mark)
    tables = [c for c, s in contrasts.items() if s["table"]]
    bigwigs = [c for c, s in contrasts.items() if s["bigwigs"]]
    if te_module:
        base = "results/te/integration/{mark}/{contrast}/"
        files = [
            "peaks/enrichment.pdf",
            "peaks/distribution.pdf",
            "atac_profile.pdf",
            "subfamilies/scatter.pdf",
            "subfamilies/full_length.pdf",
        ]
        table_files, bigwig_files = ["peaks/scatter.pdf"], ["peaks/signal_profile.pdf"]
    else:
        base = "results/integration/{mark}/{contrast}/"
        files = ["enrichment.pdf", "distribution.pdf"]
        table_files, bigwig_files = ["scatter.pdf"], ["signal_profile.pdf"]
    targets = []
    for selected, names in (
        (list(contrasts), files),
        (tables, table_files),
        (bigwigs, bigwig_files),
    ):
        targets += expand([base + f for f in names], mark=mark, contrast=selected)
    return targets


def epi_conditions(contrast):
    """[control, treatment] of a contrast"""
    treatment, control = deseq2_contrasts()[contrast]
    return [control, treatment]


def multiomics_config():
    """
    Returns the settings of the integration of all data sets of a contrast
    (multiomics: in the config), with defaults
    """
    section = config.get("multiomics") or {}
    return {
        "enable": section.get("enable", False),
        "min_copies": section.get("min_copies", 20),
        "max_subfamilies": section.get("max_subfamilies", 60),
        "max_combinations": section.get("max_combinations", 25),
    }


def multiomics_contrasts():
    """
    DESeq2 contrasts with all data sets: RNA-seq (rnaseq: tables), and
    H3K9me3 and DNA methylation with BigWig files (h3k9me3: and methylation:
    contrasts). Empty without multiomics: enable; needs the TE module
    """
    if not multiomics_config()["enable"]:
        return []
    if not te_enabled():
        raise ValueError("multiomics: enable: true needs te: enable: true")
    rnaseq = rnaseq_contrasts()
    marks = {mark: epi_contrasts(mark) for mark in EPI_MARKS}
    contrasts = [
        contrast
        for contrast in rnaseq
        if all(
            contrast in marks[mark] and marks[mark][contrast]["bigwigs"]
            for mark in EPI_MARKS
        )
    ]
    if not contrasts:
        raise ValueError(
            "multiomics: enable: true, but no contrast has an RNA-seq table "
            "(rnaseq: tables) and H3K9me3 and DNA methylation data with "
            "bigwigs (h3k9me3: and methylation: contrasts)"
        )
    # deepTools needs the profile lengths to be multiples of the bin size
    bin_size = computematrix_bin_size()
    for setting in ("flank", "body_length"):
        value = te_profile_config()[setting]
        if value % bin_size:
            raise ValueError(
                f"te: profile: {setting} ({value}) must be a multiple of the "
                f"computeMatrix bin size ({bin_size} bp) for multiomics:"
            )
    return contrasts


def multiomics_marks(contrast):
    """
    Per mark (h3k9me3, methylation): its change codes in the data of a
    contrast (decrease first), the code of a decrease, and the plot labels
    """
    return {
        mark: {
            "changes": epi_changes(mark, contrast),
            "down": EPI_MARKS[mark]["down"],
            "labels": epi_labels(mark, contrast),
        }
        for mark in EPI_MARKS
    }


def multiomics_targets():
    """Output files of the integration of all data sets"""
    files = [
        "subfamilies/heatmap.pdf",
        "peaks/upset.pdf",
        "peaks/enrichment.pdf",
    ]
    # Full-length copies only with te: profile: subfamilies
    if te_profile_config()["subfamilies"]:
        files += ["copies/heatmap.pdf", "copies/profile.pdf"]
    return expand(
        "results/te/integration/multiomics/{contrast}/{file}",
        contrast=multiomics_contrasts(),
        file=files,
    )


def te_only():
    """
    With te: only, the workflow runs the TE module only (and the read trimming,
    genome and index it needs), not the ENCODE peaks, signal tracks and QC
    """
    return te_enabled() and config.get("te", {}).get("only", False)


def te_heatmap_config():
    """
    Returns the selection for the TE heatmap of each DESeq2 contrast: lists of
    classes and subfamilies ("all" selects everything; a single name is
    accepted too) and the maximum number of subfamilies per class
    """
    heatmap = config.get("te", {}).get("heatmap", {})

    def as_list(value):
        return [value] if isinstance(value, str) else list(value)

    return {
        "classes": as_list(heatmap.get("classes", ["LINE", "SINE", "LTR", "DNA"])),
        "subfamilies": as_list(heatmap.get("subfamilies", ["all"])),
        "max_per_class": heatmap.get("max_per_class", 25),
    }


def methods_envs():
    """
    Returns the conda environments (workflow/envs/{env}.yaml) used by the
    jobs of this configuration, whose package versions go into the methods
    """
    genome = config["genome"]["ensembl"]
    # Genome, index, trimming, FastQC/MultiQC, TSS file and deepTools/Picard
    envs = ["atac", "encode", "reports", "qc", "tss"]
    if not te_only() and config["idr"]["enable"]:
        envs.append("idr")
    # ChIPseeker annotation (ENCODE peaks, or the TE peaks with DESeq2) and
    # the AnnotationHub blacklist of genomes without an ENCODE exclusion list
    blacklist_from_r = (
        genome not in LIFTED_BLACKLISTS and genome not in ENCODE_BLACKLISTS
    )
    if not te_only() or deseq2_enabled() or blacklist_from_r:
        envs.append("r")
    if genome in LIFTED_BLACKLISTS or ENCODE_TSS.get(genome, (None, None))[1]:
        envs.append("liftover")
    if deseq2_enabled() or te_enabled():
        envs.append("deseq2")
    if te_enabled():
        envs += ["te", "gat"]
        if (
            deseq2_enabled()
            and deseq2_contrasts()
            and te_correlation_config()["datasets"]
        ):
            envs.append("eulerr")
    # Plots of the integration with RNA-seq
    if rnaseq_contrasts():
        envs.append("integration")
    # Integration with H3K9me3 and DNA methylation: overlaps (bedtools), plots
    if epi_marks():
        envs += [e for e in ("integration", "gat") if e not in envs]
    # Integration of all data sets: ComplexHeatmap
    if multiomics_contrasts():
        envs.append("multiomics")
    return envs


def methods_wrappers():
    """{name: (wrapper, [packages])} of the Snakemake wrappers used (none)"""
    return {}


def methods_git():
    """{name: clone directory} of tools cloned at run time (none)"""
    return {}


def methods_string_db():
    """Output files of STRING queries (none)"""
    return []


def methods_targets():
    """All other targets: the methods describe the finished analysis"""
    return ENCODE_TARGETS + TE_TARGETS


def reference_files():
    """
    Returns the sources of the genome, blacklist and TSS files of the genome,
    for the methods
    """
    genome = config["genome"]["ensembl"]
    if genome in LIFTED_BLACKLISTS:
        accession, chain = LIFTED_BLACKLISTS[genome]
        blacklist = {
            "source": "encode",
            "accession": accession,
            "lifted_from": chain.split("/")[0],
        }
    elif genome in ENCODE_BLACKLISTS:
        blacklist = {"source": "encode", "accession": ENCODE_BLACKLISTS[genome]}
    else:
        blacklist = {"source": "annotationhub"}
    if genome in ENCODE_TSS:
        accession, chain = ENCODE_TSS[genome]
        tss = {"source": "encode", "accession": accession}
        if chain:
            tss["lifted_from"] = chain.split("/")[0]
    else:
        tss = {"source": "gtf"}
    return {
        "fasta": os.path.basename(resources.fasta),
        "gtf": os.path.basename(resources.gtf),
        "release": config["genome"]["release"],
        "ucsc_genome": UCSC_GENOMES.get(genome),
        "blacklist": blacklist,
        "tss": tss,
        "annotation_regions": bool(annotation_regions()),
    }


def env_packages(env):
    """
    Package names in the dependency lines of workflow/envs/{env}.yaml
    (bioconda::macs2=2.2.9.1=py310h7c593f9_1 -> macs2). Read while the
    workflow is parsed: workflow.source_path() gives a local copy of the file,
    also when the workflow is used as a module (e.g. deployed with
    snakedeploy), where workflow.basedir is the deployment directory
    """
    names = []
    with open(workflow.source_path(f"../envs/{env}.yaml")) as fh:
        for line in fh:
            match = re.match(r"\s*-\s+(?:[^:\s]+::)?([^=<>!\s#]+)", line)
            if match:
                names.append(match.group(1))
    return names
