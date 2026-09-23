# Assembly Quality Assessment

This directory contains workflow-level quality assessment for completed whole-genome assemblies.

The assessment is intentionally separated from `assembly_analysis/`:

- `assemblers/whole_genome_asm/assessment/` **runs assessment tools**
- `assembly_analysis/` **aggregates metrics, performs downstream analysis and creates figures**

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. Directories remain compact; use this module's curated sections for canonical files.

- [`README.md`](README.md)
- [`assembly_quality_busco.smk`](assembly_quality_busco.smk)
- [`assembly_quality_quast.smk`](assembly_quality_quast.smk)

<!-- AUTO_REPOSITORY_TREE_END -->

## Input contract

Assessment workflows discover completed assemblies using the standardized layout:

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

The assessed methods can include:

- Flye
- GoldRush
- Verkko


ntLink is not included in the current assembly-quality comparison.

## Complementary assessment tools

| Tool | Main scientific role |
|---|---|
| QUAST / QUAST-LG | contiguity and reference-based structural agreement |
| BUSCO | conserved-gene completeness |
| Merqury | k-mer completeness and consensus quality (QV) |

No single metric should be used as the sole assembly-quality criterion.

## QUAST

Workflow:

```text
assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk
```

Expected output pattern:

```text
assembly_quality/quast/{assembler}/{dataset}/report.tsv
assembly_quality/quast/{assembler}/{dataset}/quast.log
```

The workflow uses QUAST 5.3.0 with the pinned Biocontainers image:

```text
quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2
```

### Reference genome

The workflow expects the GRCh38 reference:

```text
/data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

Before execution, verify the reference path and the assembly to be assessed:

```bash
ls -lh /data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.*
ls -lh assemblies/{assembler}/{dataset}/assembly.fasta
```


### Server layout

The current production layout is:

Flye assemblies:
```text
  /data/genmedbfx/yu_j/lrs_benchmarking/assemblies/flye/
```

GoldRush assemblies:
```text
  /home/stoiber_l/smbshare/lrs_benchmarking/assemblies/goldrush/
```

Verkko assemblies:
```text
  /home/stoiber_l/smbshare/lrs_benchmarking/assemblies/verkko/
```

Genome reference:
```text
  /data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/
  GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```

QUAST workflow and output repository:
```text
  /data/genmedbfx/yu_j/lrs_benchmarking
```

### Server execution

Run Snakemake from the actual repository root on the machine that owns the Docker daemon:

```bash
cd /data/genmedbfx/yu_j/lrs_benchmarking
```

Verify the working directory:
pwd

Expected:
/data/genmedbfx/yu_j/lrs_benchmarking


Avoid running this workflow from inside an additional Docker container and then launching QUAST through another `docker run`. A path visible inside the outer container may not represent the same path on the Docker host, which can produce invalid bind mounts and missing input files inside the QUAST container.

For example, a container-visible repository path such as `/lrs_benchmarking` may correspond to a different host-side path such as `/data/.../lrs_benchmarking`.

The workflow should be launched from the host environment that owns the Docker
daemon. Avoid launching the workflow inside an additional Docker container and
then invoking another docker run, because host and container paths may not
refer to the same filesystem locations.

Input checks

Before running QUAST, the production assembly locations and genome reference
can be checked with:

```bash
find /data/genmedbfx/yu_j/lrs_benchmarking/assemblies/flye -type f -name 'assembly.fasta' -print

find /home/stoiber_l/smbshare/lrs_benchmarking/assemblies/goldrush -type f -name 'assembly.fasta' -print

find /home/stoiber_l/smbshare/lrs_benchmarking/assemblies/verkko -type f -name 'assembly.fasta' -print

ls -lh /data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
```



### Dry run

From:
```bash
cd /data/genmedbfx/yu_j/lrs_benchmarking
```

run:
```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk --cores 16 --resources mem_gb=128 --config reference=/data/genmedbfx/schilling_m/repos/lrs_benc
```

A `MissingInputException` for the reference means that the configured FASTA path is not accessible from the Snakemake execution environment.

Docker exit code `125` together with missing repository files inside the QUAST container generally indicates a bind-mount or execution-path problem rather than a QUAST installation problem.


### Real Run

The '--dry-run' command is removed to proceed with the real run based on the reference path in genmedbfx. Any other user changes such as user should be adjust in the directory path.

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_quast.smk --cores 16 --resources mem_gb=128 --config reference=/data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta --printshellcmds
```

## BUSCO

Workflow when present:

```text
assembly_quality_busco.smk
```

For human whole-genome assemblies, the F2 workflow is designed around a primate BUSCO lineage. The lineage directory must be supplied explicitly and should contain a valid `dataset.cfg`.

### BUSCO lineage

Human whole-genome assemblies are evaluated using the BUSCO lineage:

```text
primates_odb12.2
```

### Dry run

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_busco.smk --cores 16 --resources mem_gb=64 --config busco_lineage=/data/genmedbfx/ref/busco/lineages/primates_odb12.2 --dry-run --printshellcmds
```

Do not copy the example lineage path literally; use the real extracted lineage directory on the execution server.

### Real Run

The '--dry-run' command is removed to proceed with the real run based on the reference path in genmedbfx.

```bash
snakemake --snakefile assemblers/whole_genome_asm/assessment/assembly_quality_busco.smk --cores 16 --resources mem_gb=64 --config busco_lineage=/data/genmedbfx/ref/busco/lineages/primates_odb12.2 --printshellcmds
```


## Merqury

Workflow when present:

```text
assembly_quality_merqury.smk
```

Merqury should use the intended trusted k-mer source for the scientific comparison. The exact input data, k-mer database construction, container image, and output contract should be documented in the workflow and here once finalized.

Do not invent missing Merqury provenance or substitute a different k-mer source silently.

## Interpretation

### Contiguity

N50 is useful descriptively but is not sufficient to establish assembly correctness. Reference-aware metrics such as NG50 / NGA50, structural discrepancies, genome fraction and duplication should be considered where appropriate.

### Completeness

BUSCO completeness and k-mer completeness measure different properties and should be reported separately.

### Verkko representation

A haplotype-resolved / diploid Verkko output can have a different total representation from a collapsed assembly. Duplicated BUSCOs, total assembly size and related metrics must therefore be interpreted in the context of assembly representation rather than treated automatically as errors.

## Downstream aggregation

Final cross-tool tables and plots belong in:

```text
assembly_analysis/
```

This keeps raw tool outputs reproducible while allowing downstream analysis scripts to evolve independently.
