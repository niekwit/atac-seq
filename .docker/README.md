# Dockerfiles for each version

This directory contains Dockerfiles associated with each release.

The Docker image derived from this file contains all Conda environments for each rule, i.e. the whole workflow is run in one image.

These images are shared via [Docker Hub](https://hub.docker.com/repository/docker/niekwit/atac-seq/general) and are generated as follows (from directory with workflow code):

```shell
snakemake --containerize > Dockerfile
docker build -t niekwit/atac-seq:v0.4.0 .
docker login
docker push niekwit/atac-seq:v0.4.0

```
