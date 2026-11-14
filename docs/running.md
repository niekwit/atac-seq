# Running the workflow

## On a single computer

Check what will be run with a dry run:

```bash
snakemake -n
```

Then run the workflow, letting Snakemake create the conda environments:

```bash
snakemake --sdm conda --cores 16 --resources mem_mb=100000
```

Each rule declares how much memory it needs (`mem_mb`), so `--resources mem_mb=100000` makes
Snakemake run only as many jobs at the same time as fit in 100 GB of RAM. Set it to somewhat
less than the RAM of your computer. The largest jobs need about 10-16 GB (genome index,
alignment, signal tracks, annotation).

Genome files and the blacklist are downloaded on the first run. The main results are the
reproducible peak sets in `results/idr/{condition}/` and the QC summary in
`results/qc/encode_qc_summary.tsv` (see [Output files](outputs.md)).

## On a compute cluster

Use an [executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/) with a profile,
for example:

```bash
snakemake --sdm conda --executor slurm --jobs 50
```

### Runtime

Each rule has a runtime limit (Snakemake's `runtime` resource, in minutes), which cluster
executors pass on as the job's time limit. The limits are about 3-5 times the runtimes for two
replicates of ~50 million read pairs (human), and are multiplied by the attempt number, so jobs
that run out of time are resubmitted with more time when the workflow is run with `--retries`
(e.g. `--retries 2`). For much deeper data, or to change a limit, override it in a profile or on
the command line, e.g. `--set-resources bowtie2_align:runtime=1440`.

### Memory and CPUs

The memory limits (`mem_mb`) are passed on as well, for clusters that schedule by memory.

For clusters that assign memory per CPU, each rule requests enough CPUs to cover its `mem_mb`
at 3,370 MiB per CPU:

- multithreaded rules through their threads;
- single-threaded rules with more memory through the `cpus_per_task` resource (SLURM's
  `--cpus-per-task`), so that locally they take one core.

On clusters with less memory per CPU, raise the CPUs with `--set-threads` (rules without
`cpus_per_task`) or `--set-resources rule:cpus_per_task=N`.
