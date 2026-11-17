"""
Writes a methods section (Markdown) of the parts of the workflow that were
run: the configuration values, the versions of the software installed in the
conda environments of the run, and references for all tools.

Only the steps enabled in the configuration are described (ENCODE peaks and
QC unless te: only, IDR, DESeq2, TE module, TE profile and TE correlation).
Tools are cited (author, year) where they are first mentioned, with their
version; the reference list contains the cited references only. The
configuration of the run (including command line overrides) is appended.
"""

import logging
import re
import os
from datetime import date

import yaml

log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
)

config = snakemake.config
params = snakemake.params
reference = params["reference"]
output = snakemake.output[0]

# Installed package versions per conda environment
# -----------------------------------------------------
versions = {}
for path in snakemake.input["versions"]:
    env = os.path.basename(path).removesuffix(".tsv")
    versions[env] = {}
    with open(path) as fh:
        for line in fh:
            if "\t" not in line:
                continue
            package, v = line.rstrip("\n").split("\t", 1)
            # conda writes the "-" of R package versions as "_" (ashr 2.2-63)
            if package.startswith("r-") and package != "r-base":
                v = v.replace("_", "-")
            versions[env][package] = v
    logging.info(f"{env}: {len(versions[env])} packages")

# References: key -> (in-text citation, full reference)
# -----------------------------------------------------
REFERENCES = {
    "amemiya2019": (
        "Amemiya et al., 2019",
        "Amemiya HM, Kundaje A, Boyle AP (2019). The ENCODE Blacklist: "
        "identification of problematic regions of the genome. Sci Rep 9:9354. "
        "https://doi.org/10.1038/s41598-019-45839-z",
    ),
    "andrews2010": (
        "Andrews, 2010",
        "Andrews S (2010). FastQC: a quality control tool for high throughput "
        "sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
    ),
    "benjamini1995": (
        "Benjamini & Hochberg, 1995",
        "Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a "
        "practical and powerful approach to multiple testing. J R Stat Soc B "
        "57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    ),
    "broad2019": (
        "Broad Institute, 2019",
        "Broad Institute (2019). Picard toolkit. "
        "https://broadinstitute.github.io/picard/",
    ),
    "danecek2021": (
        "Danecek et al., 2021",
        "Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, "
        "Whitwham A, Keane T, McCarthy SA, Davies RM, Li H (2021). Twelve years "
        "of SAMtools and BCFtools. GigaScience 10:giab008. "
        "https://doi.org/10.1093/gigascience/giab008",
    ),
    "dyer2025": (
        "Dyer et al., 2025",
        "Dyer SC, Austine-Orimoloye O, Azov AG, et al. (2025). Ensembl 2025. "
        "Nucleic Acids Res 53:D948–D957. https://doi.org/10.1093/nar/gkae1071",
    ),
    "encode2020": (
        "ENCODE Project Consortium, 2020",
        "ENCODE Project Consortium, Moore JE, Purcaro MJ, Pratt HE, et al. "
        "(2020). Expanded encyclopaedias of DNA elements in the human and mouse "
        "genomes. Nature 583:699–710. https://doi.org/10.1038/s41586-020-2493-4; "
        "ENCODE ATAC-seq pipeline v2.2.3: "
        "https://github.com/ENCODE-DCC/atac-seq-pipeline",
    ),
    "ewels2016": (
        "Ewels et al., 2016",
        "Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize "
        "analysis results for multiple tools and samples in a single report. "
        "Bioinformatics 32:3047–3048. https://doi.org/10.1093/bioinformatics/btw354",
    ),
    "harris2020": (
        "Harris et al., 2020",
        "Harris CR, Millman KJ, van der Walt SJ, et al. (2020). Array "
        "programming with NumPy. Nature 585:357–362. "
        "https://doi.org/10.1038/s41586-020-2649-2",
    ),
    "heger2013": (
        "Heger et al., 2013",
        "Heger A, Webber C, Goodson M, Ponting CP, Lunter G (2013). GAT: a "
        "simulation framework for testing the association of genomic intervals. "
        "Bioinformatics 29:2046–2048. https://doi.org/10.1093/bioinformatics/btt343",
    ),
    "hinrichs2006": (
        "Hinrichs et al., 2006",
        "Hinrichs AS, Karolchik D, Baertsch R, et al. (2006). The UCSC Genome "
        "Browser Database: update 2006. Nucleic Acids Res 34:D590–D598. "
        "https://doi.org/10.1093/nar/gkj144",
    ),
    "hunter2007": (
        "Hunter, 2007",
        "Hunter JD (2007). Matplotlib: a 2D graphics environment. Comput Sci Eng "
        "9:90–95. https://doi.org/10.1109/MCSE.2007.55",
    ),
    "kent2002": (
        "Kent et al., 2002",
        "Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, "
        "Haussler D (2002). The human genome browser at UCSC. Genome Res "
        "12:996–1006. https://doi.org/10.1101/gr.229102",
    ),
    "kent2010": (
        "Kent et al., 2010",
        "Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D (2010). BigWig and "
        "BigBed: enabling browsing of large distributed datasets. Bioinformatics "
        "26:2204–2207. https://doi.org/10.1093/bioinformatics/btq351",
    ),
    "kolde": (
        "Kolde",
        "Kolde R. pheatmap: Pretty Heatmaps. R package. "
        "https://CRAN.R-project.org/package=pheatmap",
    ),
    "langmead2012": (
        "Langmead & Salzberg, 2012",
        "Langmead B, Salzberg SL (2012). Fast gapped-read alignment with Bowtie "
        "2. Nat Methods 9:357–359. https://doi.org/10.1038/nmeth.1923",
    ),
    "larsson": (
        "Larsson",
        "Larsson J. eulerr: Area-Proportional Euler and Venn Diagrams with "
        "Ellipses. R package. https://CRAN.R-project.org/package=eulerr",
    ),
    "lawrence2009": (
        "Lawrence et al., 2009",
        "Lawrence M, Gentleman R, Carey V (2009). rtracklayer: an R package for "
        "interfacing with genome browsers. Bioinformatics 25:1841–1842. "
        "https://doi.org/10.1093/bioinformatics/btp328",
    ),
    "lawrence2013": (
        "Lawrence et al., 2013",
        "Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan "
        "MT, Carey VJ (2013). Software for computing and annotating genomic "
        "ranges. PLoS Comput Biol 9:e1003118. "
        "https://doi.org/10.1371/journal.pcbi.1003118",
    ),
    "li2011": (
        "Li et al., 2011",
        "Li Q, Brown JB, Huang H, Bickel PJ (2011). Measuring reproducibility of "
        "high-throughput experiments. Ann Appl Stat 5:1752–1779. "
        "https://doi.org/10.1214/11-AOAS466",
    ),
    "liao2014": (
        "Liao et al., 2014",
        "Liao Y, Smyth GK, Shi W (2014). featureCounts: an efficient general "
        "purpose program for assigning sequence reads to genomic features. "
        "Bioinformatics 30:923–930. https://doi.org/10.1093/bioinformatics/btt656",
    ),
    "love2014": (
        "Love et al., 2014",
        "Love MI, Huber W, Anders S (2014). Moderated estimation of fold change "
        "and dispersion for RNA-seq data with DESeq2. Genome Biol 15:550. "
        "https://doi.org/10.1186/s13059-014-0550-8",
    ),
    "martin2011": (
        "Martin, 2011",
        "Martin M (2011). Cutadapt removes adapter sequences from high-throughput "
        "sequencing reads. EMBnet.journal 17:10–12. "
        "https://doi.org/10.14806/ej.17.1.200",
    ),
    "mckinney2010": (
        "McKinney, 2010",
        "McKinney W (2010). Data structures for statistical computing in Python. "
        "Proc 9th Python in Science Conference, 56–61. "
        "https://doi.org/10.25080/Majora-92bf1922-00a",
    ),
    "molder2021": (
        "Mölder et al., 2021",
        "Mölder F, Jablonski KP, Letcher B, Hall MB, Tomkins-Tinch CH, Sochat V, "
        "Forster J, Lee S, Twardziok SO, Kanitz A, Wilm A, Holtgrewe M, Rahmann "
        "S, Nahnsen S, Köster J (2021). Sustainable data analysis with Snakemake. "
        "F1000Research 10:33. https://doi.org/10.12688/f1000research.29032.2",
    ),
    "morgan": (
        "Morgan & Shepherd",
        "Morgan M, Shepherd L. AnnotationHub: Client to access AnnotationHub "
        "resources. Bioconductor package. "
        "https://doi.org/10.18129/B9.bioc.AnnotationHub",
    ),
    "orchard2020": (
        "Orchard et al., 2020",
        "Orchard P, Kyono Y, Hensley J, Kitzman JO, Parker SCJ (2020). "
        "Quantification, dynamic visualization, and validation of bias in "
        "ATAC-seq data with ataqv. Cell Syst 10:298–306. "
        "https://doi.org/10.1016/j.cels.2020.02.009",
    ),
    "pysam": (
        "pysam developers",
        "pysam developers. pysam: a Python module for reading and manipulating "
        "SAM/BAM files. https://github.com/pysam-developers/pysam",
    ),
    "python": (
        "Python Software Foundation",
        "Python Software Foundation. Python language reference. "
        "https://www.python.org",
    ),
    "r": (
        "R Core Team",
        "R Core Team. R: A language and environment for statistical computing. "
        "R Foundation for Statistical Computing, Vienna, Austria. "
        "https://www.R-project.org",
    ),
    "ramirez2016": (
        "Ramírez et al., 2016",
        "Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne "
        "S, Dündar F, Manke T (2016). deepTools2: a next generation web server "
        "for deep-sequencing data analysis. Nucleic Acids Res 44:W160–W165. "
        "https://doi.org/10.1093/nar/gkw257",
    ),
    "quinlan2010": (
        "Quinlan & Hall, 2010",
        "Quinlan AR, Hall IM (2010). BEDTools: a flexible suite of utilities for "
        "comparing genomic features. Bioinformatics 26:841–842. "
        "https://doi.org/10.1093/bioinformatics/btq033",
    ),
    "slowikowski": (
        "Slowikowski",
        "Slowikowski K. ggrepel: Automatically Position Non-Overlapping Text "
        "Labels with 'ggplot2'. R package. "
        "https://CRAN.R-project.org/package=ggrepel",
    ),
    "smit2015": (
        "Smit et al., 2013–2015",
        "Smit AFA, Hubley R, Green P (2013–2015). RepeatMasker Open-4.0. "
        "http://www.repeatmasker.org",
    ),
    "stephens2017": (
        "Stephens, 2017",
        "Stephens M (2017). False discovery rates: a new deal. Biostatistics "
        "18:275–294. https://doi.org/10.1093/biostatistics/kxw041",
    ),
    "virtanen2020": (
        "Virtanen et al., 2020",
        "Virtanen P, Gommers R, Oliphant TE, et al. (2020). SciPy 1.0: "
        "fundamental algorithms for scientific computing in Python. Nat Methods "
        "17:261–272. https://doi.org/10.1038/s41592-019-0686-2",
    ),
    "wickham2016": (
        "Wickham, 2016",
        "Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. "
        "Springer, New York. https://ggplot2.tidyverse.org",
    ),
    "wickham2019": (
        "Wickham et al., 2019",
        "Wickham H, Averick M, Bryan J, et al. (2019). Welcome to the tidyverse. "
        "J Open Source Softw 4:1686. https://doi.org/10.21105/joss.01686",
    ),
    "yu2015": (
        "Yu et al., 2015",
        "Yu G, Wang LG, He QY (2015). ChIPseeker: an R/Bioconductor package for "
        "ChIP peak annotation, comparison and visualization. Bioinformatics "
        "31:2382–2383. https://doi.org/10.1093/bioinformatics/btv145",
    ),
    "zhang2008": (
        "Zhang et al., 2008",
        "Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum "
        "C, Myers RM, Brown M, Li W, Liu XS (2008). Model-based analysis of "
        "ChIP-Seq (MACS). Genome Biol 9:R137. "
        "https://doi.org/10.1186/gb-2008-9-9-r137",
    ),
}

# Tools: name -> (references, [(conda environment, conda package)])
# -----------------------------------------------------
TOOLS = {
    "Snakemake": (["molder2021"], []),
    "FastQC": (["andrews2010"], [("reports", "fastqc")]),
    "cutadapt": (["martin2011"], [("encode", "cutadapt")]),
    "Bowtie 2": (["langmead2012"], [("encode", "bowtie2"), ("te", "bowtie2")]),
    "SAMtools": (["danecek2021"], [("encode", "samtools"), ("te", "samtools")]),
    "pysam": (["pysam"], [("encode", "pysam"), ("tss", "pysam"), ("te", "pysam")]),
    "Picard": (["broad2019"], [("qc", "picard")]),
    "BEDTools": (["quinlan2010"], [("encode", "bedtools"), ("gat", "bedtools")]),
    "MACS2": (["zhang2008"], [("encode", "macs2")]),
    "IDR": (["li2011"], [("idr", "idr")]),
    "UCSC bedClip and bedGraphToBigWig": (
        ["kent2010"],
        [("encode", "ucsc-bedgraphtobigwig")],
    ),
    "liftOver": (["hinrichs2006"], [("liftover", "ucsc-liftover")]),
    "ataqv": (["orchard2020"], [("atac", "ataqv")]),
    "deepTools": (["ramirez2016"], [("qc", "deeptools")]),
    "MultiQC": (["ewels2016"], [("reports", "multiqc")]),
    "rtracklayer": (["lawrence2009"], [("r", "bioconductor-rtracklayer")]),
    "txdbmaker": (["lawrence2013"], [("r", "bioconductor-txdbmaker")]),
    "ChIPseeker": (["yu2015"], [("r", "bioconductor-chipseeker")]),
    "AnnotationHub": (["morgan"], [("r", "bioconductor-annotationhub")]),
    "featureCounts": (
        ["liao2014"],
        [("deseq2", "subread"), ("te", "subread")],
    ),
    "DESeq2": (["love2014"], [("deseq2", "bioconductor-deseq2")]),
    "ashr": (["stephens2017"], [("deseq2", "r-ashr")]),
    "GAT": (["heger2013"], [("gat", "gat")]),
    "eulerr": (["larsson"], [("eulerr", "r-eulerr")]),
    "Python": (
        ["python"],
        [(env, "python") for env in ["encode", "qc", "tss", "te", "gat"]],
    ),
    "NumPy": (["harris2020"], [("qc", "numpy"), ("tss", "numpy"), ("te", "numpy")]),
    "pandas": (
        ["mckinney2010"],
        [("qc", "pandas"), ("te", "pandas"), ("gat", "pandas")],
    ),
    "SciPy": (["virtanen2020"], [("qc", "scipy")]),
    "Matplotlib": (
        ["hunter2007"],
        [("qc", "matplotlib-base"), ("tss", "matplotlib-base")],
    ),
    "R": (["r"], [(env, "r-base") for env in ["r", "deseq2", "eulerr"]]),
    "tidyverse": (["wickham2019"], [("r", "r-tidyverse")]),
    "ggplot2": (["wickham2016"], [("deseq2", "r-ggplot2")]),
    "ggrepel": (["slowikowski"], [("deseq2", "r-ggrepel")]),
    "pheatmap": (["kolde"], [("deseq2", "r-pheatmap")]),
}

cited = []
mentioned = set()


def citation(*keys):
    """In-text citation (without parentheses); adds it to the reference list"""
    for key in keys:
        if key not in cited:
            cited.append(key)
    return "; ".join(REFERENCES[key][0] for key in keys)


def cite(*keys):
    return f"({citation(*keys)})"


def version_key(v):
    """Sort key for versions (3.10.14 after 3.9.2)"""
    return [
        (0, int(x), "") if x.isdigit() else (1, 0, x) for x in re.split(r"[.\-_]", v)
    ]


def version(name):
    """Versions of a tool in the conda environments of the run, sorted"""
    if name == "Snakemake":
        return [params["snakemake_version"]]
    found = {versions.get(env, {}).get(package) for env, package in TOOLS[name][1]}
    return sorted((v for v in found if v), key=version_key)


def tool(name, detail=None):
    """
    Tool with its version(s) and citation at the first mention, the name only
    after that; detail (e.g. options) goes in the same parentheses
    """
    if name in mentioned:
        return f"{name} ({detail})" if detail else name
    mentioned.add(name)
    found = version(name)
    if not found:
        logging.warning(f"No version found for {name}")
    text = f"{name} v{'/'.join(found)}" if found else name
    inside = "; ".join(x for x in [detail, citation(*TOOLS[name][0])] if x)
    return f"{text} ({inside})"


def code(value):
    return f"`{value}`"


def number(value):
    """Thousands separators for integers"""
    return f"{value:,}" if isinstance(value, int) else str(value)


def join(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


# Settings of the run
# -----------------------------------------------------
te_only = params["te_only"]
te = params["te"]
deseq2 = params["deseq2"]
idr = config["idr"]["enable"] and not te_only
replicates = params["replicates"]
replicated = [c for c, n in replicates.items() if n > 1]
contrasts = params["contrasts"]
filter_chrs = config["filter"]["filter_chrs"]
min_fraction = f"{config['peaks'].get('min_replicate_fraction', 1):.0%}"
fdr = config.get("deseq2", {}).get("fdr", 0.05)
lfc_threshold = config.get("deseq2", {}).get("lfc_threshold", 0)
count_options = "-p -F SAF --read2pos 5 --largestOverlap -s 0"
blacklist = reference["blacklist"]
tss = reference["tss"]

sections = []

# Overview
# -----------------------------------------------------
version_text = f" {params['workflow_version']}" if params["workflow_version"] else ""
samples_text = join(f"{c} (n = {n})" for c, n in replicates.items())
encode_text = (
    ""
    if te_only
    else (
        ", following the ENCODE ATAC-seq pipeline v2.2.3 "
        f"{cite('encode2020')} unless stated otherwise"
    )
)
sections.append(
    "Paired-end ATAC-seq data of "
    f"{sum(replicates.values())} samples ({samples_text}) were analysed with "
    f"the atac-seq workflow{version_text} "
    "(https://github.com/niekwit/atac-seq) in "
    f"{tool('Snakemake')}{encode_text}. The genome "
    f"({reference['fasta']}) and gene annotation ({reference['gtf']}) were "
    f"from Ensembl release {reference['release']} {cite('dyer2025')}; "
    "non-nucleotide characters in the genome sequence were replaced by N."
)

# Read processing and alignment
# -----------------------------------------------------
cutadapt = config["cutadapt"]
text = (
    f"## Read processing{'' if te_only else ' and alignment'}\n\n"
    f"Read quality was assessed with {tool('FastQC')}. Adapters "
    f"(read 1: {code(cutadapt['adapter_r1'])}, read 2: "
    f"{code(cutadapt['adapter_r2'])}) were trimmed with "
    f"{tool('cutadapt', code(cutadapt['extra']))}."
)
if not te_only:
    multimapping = config["bowtie2"]["multimapping"]
    k = f"-k {multimapping + 1} " if multimapping else ""
    text += (
        " Read pairs were aligned with "
        f"{tool('Bowtie 2', code(k + config['bowtie2']['extra']))}."
    )
    if multimapping:
        text += (
            f" Reads with more than {multimapping} alignments were removed and "
            "the others were assigned one of their alignments at random "
            "(ENCODE assign_multimappers);"
        )
    else:
        text += (
            " Reads with a mapping quality below "
            f"{config['filter']['mapq_thresh']} were removed;"
        )
    text += (
        " only properly paired primary alignments were kept, with "
        f"{tool('SAMtools')} via {tool('pysam')}. Duplicates were marked "
        f"with {tool('Picard')} MarkDuplicates and removed, together with "
        f"alignments on {join(filter_chrs)} "
        f"({code('samtools view -F 1804 -f 2')})."
    )
sections.append(text)

# Peak calling
# -----------------------------------------------------
macs2 = config["macs2"]
smooth_win = macs2["smooth_win"]
macs2_options = (
    f"callpeak -f BED --nomodel --shift {-int(round(smooth_win / 2))} "
    f"--extsize {smooth_win} -p {macs2['pval_thresh']} --keep-dup all "
    "--call-summits -B --SPMR"
)
tagalign_text = (
    "Tn5-shifted tagAlign files (+4 bp on the plus strand, −5 bp on the minus "
    "strand)"
)
keep_chr_text = code("^(" + config["peaks"]["keep_chr_regex"] + ")$")


def blacklist_text():
    """Source of the blacklist (cites the tools only when it is used)"""
    if blacklist["source"] == "encode":
        text = f"the ENCODE blacklist {blacklist['accession']} {cite('amemiya2019')}"
        if "lifted_from" in blacklist:
            text += (
                f", lifted over from {blacklist['lifted_from']} with "
                f"{tool('liftOver', code('-minMatch=0.5'))},"
            )
        return text
    return f"the blacklist {cite('amemiya2019')} from {tool('AnnotationHub')}"


def chipseeker_text():
    return (
        f"{tool('ChIPseeker', 'TSS region ± 3 kb')}, using a TxDb built from the "
        f"Ensembl GTF with {tool('txdbmaker')} and {tool('rtracklayer')}"
    )


def macs2_text(gsize):
    options = f"{code(macs2_options)}, effective genome size {gsize}"
    return (
        f"Peaks were called with {tool('MACS2', options)}, keeping the "
        f"{number(macs2['cap_num_peak'])} peaks with the lowest p-values. Peaks "
        f"overlapping {blacklist_text()} or on chromosomes not matching "
        f"{keep_chr_text} were removed with {tool('BEDTools')}."
    )


if not te_only:
    seed = config["pseudoreplication_random_seed"]
    seed_text = (
        "seeded with the tagAlign file size" if seed == 0 else f"random seed {seed}"
    )
    text = f"## Peak calling\n\nAlignments were converted to {tagalign_text}"
    if config["subsample_reads"]:
        text += f" and subsampled to {number(config['subsample_reads'])} reads"
    text += (
        f". Each replicate was split into two pseudoreplicates ({seed_text}), "
        "and the replicates and pseudoreplicates of each condition were pooled. "
        + macs2_text(code(macs2["gsize"]))
        + " Reproducible peaks of each condition were selected by "
        "naive overlap (pooled peaks overlapping a peak of both replicates or "
        "pseudoreplicates by ≥ 50%)"
    )
    if idr:
        text += " and by " + tool(
            "IDR", f"IDR < {config['idr']['threshold']}, peaks ranked by p-value"
        )
    text += (
        ", and optimal and conservative peak sets were chosen from the "
        "comparisons of true replicates and pseudoreplicates as in ENCODE."
    )
    if replicated:
        text += (
            " In addition, pooled peaks reproducible in at least "
            f"{min_fraction} of the replicates of a condition (rounded up, at "
            "least 2) were reported (overlap ≥ 50%"
            + (", or IDR below the threshold for all pairs" if idr else "")
            + ")."
        )
    text += (
        " Fold-enrichment and p-value signal tracks were generated with MACS2 "
        f"bdgcmp and {tool('UCSC bedClip and bedGraphToBigWig')}. Optimal peaks "
        f"were annotated with {chipseeker_text()}."
    )
    sections.append(text)


# Quality control
# -----------------------------------------------------
def tss_text():
    """Source of the TSSs (cites liftOver only when it is used)"""
    if tss["source"] == "encode":
        text = f"TSSs of ENCODE file {tss['accession']}"
        if "lifted_from" in tss:
            text += f", lifted over from {tss['lifted_from']} with {tool('liftOver')}"
        return text
    return "5' ends of protein-coding genes of the Ensembl GTF"


if not te_only:
    regions_text = (
        "ENCODE DNase I hypersensitive sites, promoters, enhancers and the blacklist"
        if reference["annotation_regions"]
        else "the blacklist"
    )
    sections.append(
        "## Quality control\n\n"
        "QC metrics were calculated as in ENCODE: library complexity (NRF, PBC1 "
        "and PBC2, without mitochondrial reads), fraction of mitochondrial "
        f"reads ({config['mito_chr_name']}), TSS enrichment ({tss_text()}), "
        "fragment length distribution and GC bias (Picard "
        "CollectInsertSizeMetrics and CollectGcBiasMetrics), fraction of reads "
        f"in {regions_text}, fraction of reads in peaks, and fingerprints and "
        "Jensen–Shannon distances of the replicates with "
        f"{tool('deepTools', 'plotFingerprint, MAPQ ≥ 30')}. Further QC was done "
        f"with {tool('ataqv')}, "
        f"and reports were combined with {tool('MultiQC')}."
    )

# Differential accessibility
# -----------------------------------------------------
pca_ntop = config.get("deseq2", {}).get("pca_ntop", 500)
lfc_text = (
    f"|log2 fold change| > {lfc_threshold}" if lfc_threshold else "log2 fold change ≠ 0"
)
contrast_text = join(f"{t} vs {c}" for t, c in contrasts.values())


def deseq2_sentences(features):
    """PCA, correlation and differential accessibility of count tables"""
    text = (
        f"With {tool('DESeq2')}, sample PCA (the {number(pca_ntop)} most "
        f"variable {features}) and Pearson correlations used "
        "variance-stabilised counts (blind to the conditions)."
    )
    if contrasts:
        text += (
            f" {features[0].upper() + features[1:]} with < 10 insertions in all "
            "samples together were removed, and differential accessibility "
            f"({contrast_text}) was tested (design {code('~ condition')}; "
            f"Wald test of {lfc_text}; Benjamini–Hochberg adjusted p < {fdr}, "
            f"{citation('benjamini1995')}); log2 fold changes were shrunken "
            f"with {tool('ashr')} "
            "for ranking and plots."
        )
    return text


if deseq2 and not te_only:
    sections.append(
        "## Differential accessibility\n\n"
        f"The optimal peak sets ({params['peak_method']}) of all conditions were "
        "merged into consensus peaks (BEDTools merge). Tn5 insertions (5' read "
        f"ends) per peak were counted with "
        f"{tool('featureCounts', code(count_options))}. "
        + deseq2_sentences("peaks")
        + " Consensus peaks were annotated with ChIPseeker as above."
    )

# Transposable elements
# -----------------------------------------------------
if te:
    te_config = config.get("te", {})
    classes = te_config.get("classes") or []
    classes_text = (
        f"repeats of the classes {join(classes)}" if classes else "all repeats"
    )
    text = (
        "## Transposable elements\n\n"
        "Trimmed read pairs were aligned again with "
        f"{tool('Bowtie 2', code('--seed 42 ' + config['bowtie2']['extra']))} in "
        "default mode, which reports one of the best alignments of multimapping "
        "reads at random. Properly paired alignments were kept, and duplicates "
        "(same fragment coordinates and strand) and alignments on "
        f"{join(filter_chrs)} were removed with {tool('pysam')} and "
        f"{tool('SAMtools')}. Read pairs were classed as "
        "multimapping when a mate had an equally good alternative alignment "
        "(XS ≥ AS) and MAPQ < 2. "
        f"The RepeatMasker annotation {cite('smit2015')} of "
        f"{reference['ucsc_genome']} was downloaded from the UCSC Genome "
        f"Browser {cite('kent2002')}, and {classes_text} were used. Tn5 "
        "insertions in repeats were counted with "
        + tool(
            "featureCounts",
            f"{code(count_options)}; with {code('-M')} for all read pairs, "
            "without for uniquely mapping read pairs",
        )
        + " and summed per subfamily, family and class, with counts per million "
        "insertions and log2 enrichment (fraction of insertions over fraction of "
        "the genome covered, pseudocount 1). CPM-normalised insertion tracks "
        "were made with "
        + tool(
            "deepTools",
            f"bamCoverage {code('--binSize 10 --normalizeUsing CPM --Offset 1')}",
        )
        + " and averaged per condition (bigwigAverage)."
    )
    profile = params["profile"]
    if profile["subfamilies"]:
        text += (
            " Insertion profiles over full-length copies (covering ≥ "
            f"{profile['min_consensus_fraction']:.0%} of the consensus) of "
            f"{join(code(s) for s in profile['subfamilies'])} were calculated "
            "with computeMatrix scale-regions (body "
            f"{number(profile['body_length'])} bp, flanks "
            f"{number(profile['flank'])} bp, 50 bp bins) and plotted with "
            "plotProfile."
        )
    sections.append(text)

    # Subfamilies with DESeq2
    if deseq2:
        normalisation = te_config.get("normalisation", "deseq2")
        normalisation_text = (
            "median-of-ratios size factors over the subfamilies"
            if normalisation == "deseq2"
            else "size factors proportional to all kept insertions of a sample"
        )
        heatmap = params["heatmap"]
        if te_only:
            text = (
                "Subfamily counts (all read pairs) were normalised with "
                f"{normalisation_text}. " + deseq2_sentences("subfamilies")
            )
        else:
            text = (
                f"Subfamily counts (all read pairs; {normalisation_text}) were "
                "analysed with DESeq2 as the consensus peaks."
            )
        if contrasts:
            selection = (
                "the significant subfamilies (at most "
                f"{heatmap['max_per_class']} per class)"
                if "all" in heatmap["subfamilies"]
                else join(code(s) for s in heatmap["subfamilies"])
            )
            heatmap_classes = (
                "all classes"
                if "all" in heatmap["classes"]
                else f"the classes {join(heatmap['classes'])}"
            )
            text += (
                " Heatmaps show row z-scores of variance-stabilised counts of "
                f"{selection} of {heatmap_classes}."
            )
        sections.append(text)

    # TE peaks
    gsize_text = "the total chromosome length"
    if te_only:
        text = (
            f"For peak calling, these alignments were converted to {tagalign_text} "
            "and pooled per condition. " + macs2_text(gsize_text)
        )
    else:
        text = (
            "Peaks were called from these alignments as above (MACS2 with "
            f"{gsize_text} as effective genome size, blacklist filtering)."
        )
    # replicated_peaks.py: k = min(n, max(2, ceil(fraction * n)))
    single = (
        ", or the only replicate of a condition with one"
        if len(replicated) < len(replicates)
        else ""
    )
    text += (
        f" Pooled peaks overlapping a peak of at least {min_fraction} of the "
        f"replicates (rounded up, at least 2{single}) by ≥ 50% were merged into "
        "consensus TE peaks (BEDTools merge). Insertions per peak were counted "
        f"with featureCounts ({code('-M')} for all read pairs, and uniquely "
        "mapping read pairs only), and each peak was annotated with the repeat "
        "of largest overlap and the fraction of insertions from uniquely mapping "
        "read pairs."
    )
    if deseq2:
        # The TE peaks always use DESeq2's median-of-ratios size factors
        # (deseq2_fit.R), whatever te: normalisation is
        text += (
            " TE peaks were analysed with DESeq2 as the "
            + (
                "subfamilies, but with median-of-ratios size factors over the " "peaks,"
                if te_only
                else "consensus peaks"
            )
            + " and annotated with "
            + (chipseeker_text() if te_only else "ChIPseeker")
            + "."
        )
    sections.append(text)

    # Correlation with other data sets
    datasets = params["datasets"]
    if deseq2 and contrasts and datasets:
        sections.append(
            "Enrichment of the more and less accessible TE peaks of each contrast "
            f"in {join(datasets)} was tested with "
            + tool("GAT", f"{number(params['gat_samples'])} samples, random seed 42")
            + f", with the genome without {join(filter_chrs)} and the blacklist "
            "as workspace and q-values over all tests of a contrast. Overlaps "
            f"were shown as area-proportional Euler diagrams with {tool('eulerr')}."
        )

    if te_only:
        sections.append(
            "QC of the TE alignments: TSS enrichment as in ENCODE "
            f"({tss_text()}), fragment length distribution and GC bias with "
            + tool("Picard", "CollectInsertSizeMetrics and CollectGcBiasMetrics")
            + f" and SAMtools stats, combined with {tool('MultiQC')}."
        )
    else:
        sections.append(
            "QC of the TE alignments (TSS enrichment, fragment length "
            "distribution, GC bias and SAMtools stats) was done as above and "
            "combined with MultiQC."
        )

# Software used by the workflow's own scripts
# -----------------------------------------------------
# R scripts: GTF annotation and ChIPseeker (ENCODE peaks, or TE peaks with
# DESeq2), the AnnotationHub blacklist, DESeq2 and eulerr
annotation = not te_only or deseq2
r_used = annotation or blacklist["source"] == "annotationhub"
r_packages = []
if annotation:
    r_packages.append(tool("tidyverse"))
if deseq2:
    r_packages += [tool("ggplot2"), tool("ggrepel"), tool("pheatmap")]
text = (
    f"Custom scripts used {tool('Python')} with {tool('NumPy')}, "
    f"{tool('pandas')}, {tool('SciPy')} and {tool('Matplotlib')}"
)
if r_used:
    text += f", and {tool('R')}" + (f" with {join(r_packages)}" if r_packages else "")
sections.append(text + ".")

# References and configuration
# -----------------------------------------------------
references = sorted(
    (REFERENCES[key][1] for key in cited), key=lambda r: r.lower().lstrip("(")
)
config_yaml = yaml.safe_dump(
    dict(config), sort_keys=False, default_flow_style=None, width=88
)

with open(output, "w") as fh:
    fh.write("# Methods\n\n")
    fh.write(
        f"*Generated on {date.today().isoformat()} from the configuration of "
        "the run (appended below). Software versions are those installed in the "
        "conda environments of the run.*\n\n"
    )
    fh.write("\n\n".join(sections) + "\n\n")
    fh.write("## References\n\n")
    fh.write("\n".join(f"- {r}" for r in references) + "\n\n")
    fh.write("## Configuration\n\n```yaml\n" + config_yaml + "```\n")

logging.info(f"Cited {len(cited)} references")
