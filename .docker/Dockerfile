FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="c4757afffe41b10b7beb68efd73a082d843ea38c18c007002ab4d16f0c8c1ded"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/atac.yaml
#   prefix: /conda-envs/3fa5a6410b3ad1a25bf4aa1f6ecdbb6c
#   name: atac
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - pigz=2.8
#       - conda-forge::python=3.12
#       - conda-forge::pip=26.0.1
#       - bioconda::ataqv=1.3.1
#       - bioconda::samtools=1.23
#       - pip:
#             - pyranges1==1.3.0
RUN mkdir -p /conda-envs/3fa5a6410b3ad1a25bf4aa1f6ecdbb6c
COPY workflow/envs/atac.yaml /conda-envs/3fa5a6410b3ad1a25bf4aa1f6ecdbb6c/environment.yaml

# Conda environment:
#   source: workflow/envs/deseq2.yaml
#   prefix: /conda-envs/febae1ec7d9f7c1e400e8f29d9008d1d
#   # Counting reads in peaks and DESeq2 analysis of all samples
#   name: deseq2
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bioconda::bioconductor-deseq2=1.50.2
#     - bioconda::subread=2.1.1
#     - conda-forge::r-ashr
#     - conda-forge::r-ggplot2
#     - conda-forge::r-ggrepel
#     - conda-forge::r-pheatmap
RUN mkdir -p /conda-envs/febae1ec7d9f7c1e400e8f29d9008d1d
COPY workflow/envs/deseq2.yaml /conda-envs/febae1ec7d9f7c1e400e8f29d9008d1d/environment.yaml

# Conda environment:
#   source: workflow/envs/encode.yaml
#   prefix: /conda-envs/f74af83bc8a7ce9a857b8b62bbf328a9
#   # Tools used by the ENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)
#   name: encode
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.10
#     - conda-forge::numpy=1.26
#     - conda-forge::pigz=2.8
#     - conda-forge::openssl
#     - bioconda::cutadapt=4.9
#     - bioconda::bowtie2=2.5.4
#     - bioconda::samtools=1.21
#     - bioconda::pysam=0.23.3
#     - bioconda::bedtools=2.31.1
#     # Later builds reference __log_finite, which glibc >= 2.31 no longer exports
#     - bioconda::macs2=2.2.9.1=py310h7c593f9_1
#     - bioconda::ucsc-bedclip=482
#     - bioconda::ucsc-bedgraphtobigwig=482
RUN mkdir -p /conda-envs/f74af83bc8a7ce9a857b8b62bbf328a9
COPY workflow/envs/encode.yaml /conda-envs/f74af83bc8a7ce9a857b8b62bbf328a9/environment.yaml

# Conda environment:
#   source: workflow/envs/eulerr.yaml
#   prefix: /conda-envs/97a4ed4484bedf357d14749c7ea21e4b
#   # Area-proportional Venn (Euler) diagrams
#   name: eulerr
#   channels:
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - conda-forge::r-base=4.4
#     - conda-forge::r-eulerr=7.0.2
#     - conda-forge::r-gridextra=2.3
RUN mkdir -p /conda-envs/97a4ed4484bedf357d14749c7ea21e4b
COPY workflow/envs/eulerr.yaml /conda-envs/97a4ed4484bedf357d14749c7ea21e4b/environment.yaml

# Conda environment:
#   source: workflow/envs/gat.yaml
#   prefix: /conda-envs/c706002875dacea712684a206a04756c
#   # GAT (Genomic Association Tester); gat 1.3.6 uses numpy.float, removed in
#   # numpy 1.24
#   name: gat
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.10
#     - conda-forge::numpy=1.23
#     - conda-forge::pandas=2.0
#     - conda-forge::scipy=1.11
#     - conda-forge::matplotlib-base=3.7
#     - bioconda::gat=1.3.6
#     - bioconda::bedtools=2.31.1
RUN mkdir -p /conda-envs/c706002875dacea712684a206a04756c
COPY workflow/envs/gat.yaml /conda-envs/c706002875dacea712684a206a04756c/environment.yaml

# Conda environment:
#   source: workflow/envs/idr.yaml
#   prefix: /conda-envs/88052609364672bb9c19593ccc476bd9
#   # IDR as run by the ENCODE ATAC-seq pipeline (https://github.com/ENCODE-DCC/atac-seq-pipeline)
#   # Separate from encode.yaml, as its numpy pin is incompatible with the working macs2 build
#   name: idr
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.10
#     # idr 2.0.4.2 uses np.int, which was removed in numpy 1.24
#     - conda-forge::numpy=1.23
#     - conda-forge::matplotlib-base=3.8
#     - bioconda::idr=2.0.4.2
#     - bioconda::ucsc-bedclip=482
RUN mkdir -p /conda-envs/88052609364672bb9c19593ccc476bd9
COPY workflow/envs/idr.yaml /conda-envs/88052609364672bb9c19593ccc476bd9/environment.yaml

# Conda environment:
#   source: workflow/envs/liftover.yaml
#   prefix: /conda-envs/ad0fd93fc15de5b62ea2f16b5370f844
#   # Lifting the ENCODE exclusion list to a newer assembly
#   name: liftover
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::pigz=2.8
#     - bioconda::ucsc-liftover=482
RUN mkdir -p /conda-envs/ad0fd93fc15de5b62ea2f16b5370f844
COPY workflow/envs/liftover.yaml /conda-envs/ad0fd93fc15de5b62ea2f16b5370f844/environment.yaml

# Conda environment:
#   source: workflow/envs/methods.yaml
#   prefix: /conda-envs/06ce3aa372654153d9fb942c443e5803
#   # Methods section of the run (snakemake-methods)
#   name: methods
#   channels:
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.12
#     - conda-forge::jinja2=3.1
#     - conda-forge::pyyaml=6.0
#     - conda-forge::git=2.47.1
#     - conda-forge::pip
#     - pip:
#         - git+https://github.com/niekwit/snakemake-methods@v0.1.1
RUN mkdir -p /conda-envs/06ce3aa372654153d9fb942c443e5803
COPY workflow/envs/methods.yaml /conda-envs/06ce3aa372654153d9fb942c443e5803/environment.yaml

# Conda environment:
#   source: workflow/envs/qc.yaml
#   prefix: /conda-envs/f7d20974e82d452ad6b06c7558708582
#   # ENCODE ATAC-seq QC: fragment length, GC bias, JSD and annotated regions
#   name: qc
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.11
#     - conda-forge::numpy=2.4
#     - conda-forge::scipy=1.17
#     - conda-forge::pandas=3.0
#     - conda-forge::matplotlib-base=3.11
#     # Picard needs R for its histogram and GC charts
#     - conda-forge::r-base=4.6
#     - bioconda::picard=3.4.0
#     - bioconda::deeptools=3.5.6
#     - bioconda::bedtools=2.31.1
#     - bioconda::samtools=1.21
RUN mkdir -p /conda-envs/f7d20974e82d452ad6b06c7558708582
COPY workflow/envs/qc.yaml /conda-envs/f7d20974e82d452ad6b06c7558708582/environment.yaml

# Conda environment:
#   source: workflow/envs/r.yaml
#   prefix: /conda-envs/309f933dc6e56cf74268a569ddd265d4
#   name: R
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bioconda::bioconductor-rtracklayer=1.70.1
#     - bioconda::bioconductor-annotationhub=4.0.0
#     - bioconda::bioconductor-genomeinfodb=1.46.2
#     - bioconda::bioconductor-txdbmaker=1.6.2
#     - r-tidyverse=2.0.0
#     - bioconductor-chipseeker=1.46.1
RUN mkdir -p /conda-envs/309f933dc6e56cf74268a569ddd265d4
COPY workflow/envs/r.yaml /conda-envs/309f933dc6e56cf74268a569ddd265d4/environment.yaml

# Conda environment:
#   source: workflow/envs/reports.yaml
#   prefix: /conda-envs/69c6246d338f3230d62ca85ccdfabd04
#   # Read QC and the MultiQC report
#   name: reports
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bioconda::fastqc=0.12.1
#     - bioconda::multiqc=1.33
RUN mkdir -p /conda-envs/69c6246d338f3230d62ca85ccdfabd04
COPY workflow/envs/reports.yaml /conda-envs/69c6246d338f3230d62ca85ccdfabd04/environment.yaml

# Conda environment:
#   source: workflow/envs/te.yaml
#   prefix: /conda-envs/5beefbf5e42bcd760dfc74810c91cc46
#   # Transposable element (TE) module: multimapper alignment and repeat counting
#   name: te
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.12
#     - conda-forge::numpy=2.2
#     - conda-forge::pandas=2.2
#     - conda-forge::pigz=2.8
#     - bioconda::bowtie2=2.5.4
#     - bioconda::samtools=1.21
#     - bioconda::pysam=0.23.3
#     - bioconda::subread=2.1.1
RUN mkdir -p /conda-envs/5beefbf5e42bcd760dfc74810c91cc46
COPY workflow/envs/te.yaml /conda-envs/5beefbf5e42bcd760dfc74810c91cc46/environment.yaml

# Conda environment:
#   source: workflow/envs/tss.yaml
#   prefix: /conda-envs/12bec9fae685a436f8f92dae48aa3c36
#   # TSS enrichment as in the ENCODE ATAC-seq pipeline
#   name: tss
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - conda-forge::python=3.12
#     - conda-forge::numpy=2.2
#     - conda-forge::matplotlib-base=3.10
#     - bioconda::pysam=0.23.3
RUN mkdir -p /conda-envs/12bec9fae685a436f8f92dae48aa3c36
COPY workflow/envs/tss.yaml /conda-envs/12bec9fae685a436f8f92dae48aa3c36/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/3fa5a6410b3ad1a25bf4aa1f6ecdbb6c --file /conda-envs/3fa5a6410b3ad1a25bf4aa1f6ecdbb6c/environment.yaml && \
    conda env create --prefix /conda-envs/febae1ec7d9f7c1e400e8f29d9008d1d --file /conda-envs/febae1ec7d9f7c1e400e8f29d9008d1d/environment.yaml && \
    conda env create --prefix /conda-envs/f74af83bc8a7ce9a857b8b62bbf328a9 --file /conda-envs/f74af83bc8a7ce9a857b8b62bbf328a9/environment.yaml && \
    conda env create --prefix /conda-envs/97a4ed4484bedf357d14749c7ea21e4b --file /conda-envs/97a4ed4484bedf357d14749c7ea21e4b/environment.yaml && \
    conda env create --prefix /conda-envs/c706002875dacea712684a206a04756c --file /conda-envs/c706002875dacea712684a206a04756c/environment.yaml && \
    conda env create --prefix /conda-envs/88052609364672bb9c19593ccc476bd9 --file /conda-envs/88052609364672bb9c19593ccc476bd9/environment.yaml && \
    conda env create --prefix /conda-envs/ad0fd93fc15de5b62ea2f16b5370f844 --file /conda-envs/ad0fd93fc15de5b62ea2f16b5370f844/environment.yaml && \
    conda env create --prefix /conda-envs/06ce3aa372654153d9fb942c443e5803 --file /conda-envs/06ce3aa372654153d9fb942c443e5803/environment.yaml && \
    conda env create --prefix /conda-envs/f7d20974e82d452ad6b06c7558708582 --file /conda-envs/f7d20974e82d452ad6b06c7558708582/environment.yaml && \
    conda env create --prefix /conda-envs/309f933dc6e56cf74268a569ddd265d4 --file /conda-envs/309f933dc6e56cf74268a569ddd265d4/environment.yaml && \
    conda env create --prefix /conda-envs/69c6246d338f3230d62ca85ccdfabd04 --file /conda-envs/69c6246d338f3230d62ca85ccdfabd04/environment.yaml && \
    conda env create --prefix /conda-envs/5beefbf5e42bcd760dfc74810c91cc46 --file /conda-envs/5beefbf5e42bcd760dfc74810c91cc46/environment.yaml && \
    conda env create --prefix /conda-envs/12bec9fae685a436f8f92dae48aa3c36 --file /conda-envs/12bec9fae685a436f8f92dae48aa3c36/environment.yaml && \
    conda clean --all -y
