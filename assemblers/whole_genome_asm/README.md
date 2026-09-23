# Production Whole-Genome Assembly Workflows

These are the production whole-genome assembly workflows for Flye, GoldRush, Verkko and ntLink.

**Run from the repository root.** The workflow files intentionally remain in this directory because they include the shared root-level `header_assembler.smk` and use repository-root-based paths.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. Directories remain compact; use this module's curated sections for canonical files.

- [`README.md`](README.md)
- [`hybrid.assembly.verkko.smk`](hybrid.assembly.verkko.smk)
- [`ont.assembly.flye2.smk`](ont.assembly.flye2.smk)
- [`ont.assembly.goldrush.smk`](ont.assembly.goldrush.smk)
- [`ont.scaffolding.ntlink.smk`](ont.scaffolding.ntlink.smk)
- [`pb.assembly.flye2.smk`](pb.assembly.flye2.smk)
- [`pb.assembly.goldrush.smk`](pb.assembly.goldrush.smk)
- [`pb.scaffolding.ntlink.smk`](pb.scaffolding.ntlink.smk)

- **assessment/** — 3 files; [`guide`](assessment/README.md)

<!-- AUTO_REPOSITORY_TREE_END -->

## Input convention

Production inputs are expected under:

```text
fastq/
```

with names such as:

```text
HG002.ont.30x.fastq.gz
HG002.pb.30x.fastq.gz
HG003.ont.30x.fastq.gz
HG003.pb.30x.fastq.gz
HG004.ont.30x.fastq.gz
HG004.pb.30x.fastq.gz
```

Reduced `1k`, `chr21`, smoke-test and local-test inputs are validation fixtures and should not be mixed into the final whole-genome 30x comparison.

## Workflows

```text
ont.assembly.flye2.smk
pb.assembly.flye2.smk
ont.assembly.goldrush.smk
pb.assembly.goldrush.smk
hybrid.assembly.verkko.smk
ont.scaffolding.ntlink.smk
pb.scaffolding.ntlink.smk
```

## Output convention

```text
assemblies/{assembler}/{dataset}/assembly.fasta
```

Examples:

```text
assemblies/flye/HG002.ont.30x/assembly.fasta
assemblies/goldrush/HG002.pb.30x/assembly.fasta
assemblies/ntlink/HG002.ont.30x/assembly.fasta
assemblies/verkko/HG002/assembly.fasta
```

## Flye

Current workflow version on the repository main branch:

```text
Flye 2.9.6
Docker: nicolasardila1/lrs-flye2:2.9.6
```

Dry-run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk --cores 32 --resources mem_gb=240 --dry-run --printshellcmds
```

Run ONT:

```bash
snakemake --snakefile assemblers/whole_genome_asm/ont.assembly.flye2.smk --cores 32 --resources mem_gb=240 --rerun-incomplete --printshellcmds --show-failed-logs
```

Dry-run / run PacBio HiFi by replacing the Snakefile with:

```text
assemblers/whole_genome_asm/pb.assembly.flye2.smk
```

## GoldRush

Use the technology-specific workflows:

```text
assemblers/whole_genome_asm/ont.assembly.goldrush.smk
assemblers/whole_genome_asm/pb.assembly.goldrush.smk
```

Always dry-run before execution. Resource values should follow the values defined by the active workflow / server configuration rather than being silently duplicated in multiple documentation locations.

## Verkko

Verkko is a hybrid assembler. It requires matching ONT and PacBio HiFi 30x FASTQs for each sample.

Current workflow version on the repository main branch:

```text
Verkko 2.3.2
Docker: nicolasardila1/lrs-verkko2:2.3.2
```

Dry-run:

```bash
snakemake --snakefile assemblers/whole_genome_asm/hybrid.assembly.verkko.smk --configfile assemblers/config/server.yaml --cores 32 --resources mem_mb=200000 --dry-run --printshellcmds
```

Run:

```bash
snakemake --snakefile assemblers/whole_genome_asm/hybrid.assembly.verkko.smk --configfile assemblers/config/server.yaml --cores 32 --resources mem_mb=200000 --rerun-incomplete --printshellcmds --show-failed-logs
```

## ntLink

ntLink is a post-assembly scaffolder. In this benchmark it should be interpreted separately from the three native assemblers.

Run only after the corresponding draft assembly exists.

Technology-specific workflows:

```text
assemblers/whole_genome_asm/ont.scaffolding.ntlink.smk
assemblers/whole_genome_asm/pb.scaffolding.ntlink.smk
```

## Quality assessment

After assemblies are complete, continue with:

[`assessment/README.md`](assessment/README.md)

The assessment layer combines complementary metrics from QUAST, BUSCO and Merqury rather than relying on N50 alone.

## Validation before commit

For every modified workflow:

```bash
git diff --check
snakemake --snakefile path/to/workflow.smk --dry-run --printshellcmds
```

Do not move these workflow files as part of a documentation-only cleanup. Their relative include paths make location a functional part of the implementation.
