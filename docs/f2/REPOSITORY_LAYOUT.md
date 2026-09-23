# F2 Repository Layout

## Design principle

The repository is organized by scientific responsibility rather than by file creation date.

```text
alignment workflows -> alignment analysis
assembly workflows  -> assembly assessment -> assembly analysis
variant workflows   -> variant-calling analysis
```

## Files intentionally kept in place

### Root mapper workflows

The active mapping workflows remain at repository root because the constitution defines their naming and execution model there.

### Whole-genome assembler workflows

`assemblers/whole_genome_asm/*.smk` remain in place because their current relative include path to `header_assembler.smk` is functional and location-sensitive.

### Shared variant-caller workflows

Existing caller workflows from the wider project remain where they are. The F2 reorganization adds a `variant_calling_analysis/` analysis layer without relocating shared workflow code.

## New canonical analysis modules

```text
alignment_analysis/
assembly_analysis/
variant_calling_analysis/
```

Each analysis module follows the same mental model:

```text
README.md
scripts/
tables/
figures/
```

## Workflow outputs versus analysis outputs

Workflow outputs remain in stable production locations:

```text
cram/
assemblies/
assembly_quality/
```

Analysis products remain under their scientific analysis module.

This separation prevents plot and metric scripts from becoming mixed with Snakemake production workflows.

## Documentation hierarchy

```text
README.md
├── docs/FILE_INDEX.md
├── docs/REPOSITORY_TREE.md
├── alignment_analysis/README.md
├── assemblers/README.md
│   └── assemblers/whole_genome_asm/README.md
│       └── assessment/README.md
├── assembly_analysis/README.md
└── variant_calling_analysis/README.md
```

The root README is navigation and scope. `FILE_INDEX.md` locates canonical files, while `REPOSITORY_TREE.md` defines directory responsibilities. Detailed run commands belong only in the most specific relevant README, reducing duplicated documentation that can become inconsistent.

## Final, exploratory, and archived outputs

Selected tables and figures use predictable status directories:

```text
final/         canonical input to reporting
derived/       reproducible intermediate tables
source/        imported source material
exploratory/   diagnostics and alternative visualizations
archive/       superseded versions retained for provenance
```
