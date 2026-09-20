# lrs_benchmarking

Snakemake workflows for benchmarking long-read and whole-genome variant calling pipelines.

## Current project status

This repository currently contains workflow files for:

- **Long-read read mapping**
  - ONT with `minimap2`
  - ONT with `pbmm2`
  - PacBio HiFi with `minimap2`
  - PacBio HiFi with `pbmm2`

- **SNV/indel calling**
  - ONT with `Clair3`
  - ONT with `DeepVariant`
  - PacBio HiFi with `Clair3`
  - PacBio HiFi with `DeepVariant`
  - Illumina with `Clair3`
  - WGS with `DeepVariant`

- **Structural variant (SV) calling**
  - `Sniffles2`
  - `cuteSV`
  - `pbsv`
  - `Sawfish`

- **Benchmarking and reporting**
  - `hap.py` + `rep.py` for small variant benchmarking
  - `Truvari` for structural variant benchmarking

At the moment, the repository is best understood as a **workflow collection / benchmarking sandbox** rather than a fully packaged, end-user-ready pipeline. The workflows are functional but still tightly coupled to the local compute/storage environment in which they were developed.

## Repository structure

### Core workflow files

#### Read mapping
- `ont.read_mapping.minimap2.smk`
- `ont.read_mapping.pbmm2.smk`
- `pb.read_mapping.minimap2.smk`
- `pb.read_mapping.pbmm2.smk`

#### SNV/indel calling
- `ont.snv_calling.clair3.smk`
- `ont.snv_calling.deepvariant.smk`
- `pb.snv_calling.clair3.smk`
- `pb.snv_calling.deepvariant.smk`
- `ilmn.snv_calling.clair3.smk`
- `wgs.snv_calling.deepvariant.smk`
- `ont.snv_filtering.v3.smk`

#### Structural variant calling
- `sniffles2.hg38.smk`
- `cuteSV.hg38.smk`
- `pbsv.hg38.smk`
- `sawfish.hg38.smk`

#### Benchmarking / evaluation
- `run_happy.smk`

#### Shared configuration
- `header.smk`

### Output / working directories currently present
- `happy_results/`
- `run_metrics/`
- `sawfish/`
- `truvari/`
- `vcf_called/`

### Other files
- `lrs_benchmarking_tools_selection.xlsx`

## Reference build and assumptions

Most workflows are currently written around **GRCh38 / hg38** resources and naming conventions.

Examples visible in the current workflow files include:
- `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta`
- GIAB-style small variant benchmarking inputs for `hap.py`
- GIAB/T2TQ100-style SV benchmarking inputs for `Truvari`

Some configuration logic in `header.smk` also supports `hs1`, but the main benchmarking workflows in this repository currently appear to be centered on **hg38**.

## Environment coupling

The workflows currently assume access to a specific local environment, including:

- fixed storage locations such as `/mnt/storage/...`
- `srun` / Slurm-based execution
- local `mamba` environments
- Docker images hosted locally or externally
- a `VARCAD_PATH` environment variable and associated VarCAD configuration/database layout

Because of this, the repository is **not yet plug-and-play** on a fresh machine or generic cluster without adaptation.

## Tools currently used in the workflows

### Read mapping
- `minimap2`
- `pbmm2`
- `samtools`

### SNV/indel calling
- `DeepVariant`
- `Clair3`
- `bcftools`

### SV calling
- `Sniffles2`
- `cuteSV`
- `pbsv`
- `Sawfish`

### Benchmarking
- `hap.py`
- `rep.py`
- `Truvari`

## What is implemented right now

### Implemented
- Separate Snakemake workflows for multiple mapper/caller combinations
- Small variant benchmarking workflow via `run_happy.smk`
- SV benchmarking steps embedded in the SV caller workflows via `Truvari`
- CPU/GPU branching for DeepVariant workflows
- Support for both FASTQ and unmapped BAM inputs in the read-mapping workflows

### Not yet standardized
- A single top-level entrypoint workflow
- Centralized user-facing configuration
- Reproducible installation instructions
- Portable paths for databases, references, and containers
- Consolidated documentation on expected inputs/outputs
- A documented end-to-end “run this project from scratch” procedure

## Expected input patterns

The current workflows infer samples from files in the working directory using patterns such as:

- `fastq/{dataset}.fastq.gz`
- `bam_unmapped/{dataset}.bam`
- `cram/{dataset}.hg38.{mapper}.cram`

Outputs are written into directories such as:

- `cram/`
- `vcf_called/snv_indel/`
- `sniffles2/`
- `cuteSV/`
- `pbsv/`
- `sawfish/`
- `truvari/`
- `happy_results/`

## Running the workflows

There is currently no single documented command for running the full project end to end.

In practice, workflows appear to be executed by calling Snakemake directly on individual `.smk` files, for example:

```bash
snakemake -s ont.read_mapping.minimap2.smk
snakemake -s ont.snv_calling.deepvariant.smk
snakemake -s sniffles2.hg38.smk
snakemake -s run_happy.smk
```

However, before doing so you will likely need to adapt:

- storage paths
- reference paths
- container/image access
- Slurm settings
- `VARCAD_PATH`
- conda/mamba environments

## Editor / VS Code setup

The repository ships a `.vscode/` config so the workflows and helper scripts get proper editor support out of the box:

- **Terminal shell integration** is enabled for the bash profile, giving command decorations, run-recent-command, and scroll-to-command in the integrated terminal — useful since most `.smk` workflows are launched via `snakemake` from the terminal and several helper scripts live under `alignment_analysis/scripts/`, `run_metrics/`, and `docs/`.
- `.smk` files get Snakemake syntax highlighting via the recommended `Snakemake.snakemake-lang` extension.
- `.sh` files get linting (`timonwong.shellcheck`) and formatting (`foxundermoon.shell-format`) on the recommended extensions list. Install `shellcheck` and `shfmt` locally (e.g. `brew install shellcheck shfmt` or your distro's package manager) for these to work.
- Opening the repo in VS Code will prompt you to install the recommended extensions (`.vscode/extensions.json`); accept the prompt or run "Extensions: Show Recommended Extensions".

## Recommended next steps for this repository

To make this project easier to use and maintain, the next improvements would likely be:

1. Add a top-level workflow or documented execution order
2. Move all hard-coded paths into a config file
3. Document required references, truth sets, and benchmark resources
4. Provide environment setup instructions
5. Clarify which workflows are production-ready vs experimental
6. Add example commands for common benchmark scenarios
7. Standardize output locations and naming conventions

## Summary

`lrs_benchmarking` currently contains a useful set of Snakemake workflows for evaluating long-read and WGS mapping/variant-calling strategies across small variants and structural variants. The code reflects active benchmarking work, but the repository is still environment-specific and would benefit from additional documentation and configuration cleanup before broader reuse.