# Chromosome 21 assembler containers

This directory contains one Dockerfile per workflow tool:

- `flye2`: Flye 2.9.6
- `goldrush`: GoldRush 1.2.2
- `ntlink`: ntLink 1.3.11
- `verkko2`: Verkko 2.3.2

Each Dockerfile extends the exact BioContainers image corresponding to the
Bioconda build recovered from the tested local Snakemake environment.

The Dockerfiles have not yet been built locally because the Docker daemon is
not available in the current WSL environment.

The Snakemake container URIs are stored separately in:

`benchmark_chr21_real/config/containers.yaml`

This allows the public Quay images to be replaced later by institutional or
GitHub Container Registry images without editing every workflow rule.
