# Assembly Analysis

This module contains downstream metric aggregation, statistical analysis and figure generation for the whole-genome assembly benchmark.

It does **not** run the assemblers themselves. Production assembly workflows are under `assemblers/whole_genome_asm/`, and external assessment tools are run under `assemblers/whole_genome_asm/assessment/`.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. Directories remain compact; use this module's curated sections for canonical files.

- [`README.md`](README.md)

- **figures/** — 1 file
- **scripts/** — 3 files; [`guide`](scripts/README.md)
- **tables/** — 1 file

<!-- AUTO_REPOSITORY_TREE_END -->

## Intended structure

```text
assembly_analysis/
├── README.md
├── scripts/
│   ├── README.md
│   ├── metrics/
│   │   ├── calculate_quality_metrics.py
│   │   ├── extract_quast_metrics.py
│   │   ├── extract_busco_metrics.py
│   │   └── extract_merqury_metrics.py
│   └── plots/
├── tables/
└── figures/
```

Only files that actually exist should be documented as implemented. The names above define the preferred location when those scripts are present.

## Raw assembly metric extraction

FASTA-derived metrics can include:

- total assembly length
- number of contigs / sequences
- largest and smallest sequence
- mean and median sequence length
- N50 / L50
- N90 / L90
- nucleotide counts

These descriptive metrics should be kept distinct from reference-aware metrics and k-mer/gene completeness metrics.

## External assessment metrics

Expected sources include:

```text
assembly_quality/quast/
assembly_quality/busco/
assembly_quality/merqury/
```

The final table should retain provenance so that each value can be traced back to a specific source file, assembler, sample, technology and assembly representation.

## Canonical summary table

Recommended final location:

```text
assembly_analysis/tables/assembler_metrics.tsv
```

A long/tidy or well-documented wide table is preferable to multiple disconnected hand-edited spreadsheets.

## Figures

Generated assembly figures belong in:

```text
assembly_analysis/figures/
```

Plotting scripts belong in:

```text
assembly_analysis/scripts/plots/
```

Recommended figure families include:

- contiguity: NG50/N50, contig count, assembly span
- reference agreement: NGA50 / misassembly-related metrics where valid
- sequence accuracy: QV / reference discordance metrics
- completeness: BUSCO and k-mer completeness
- computational resources: runtime and peak RAM
- ntLink before/after comparisons

## Scientific comparison rules

1. Compare the same biological sample and input regime whenever possible.
2. Keep ONT, PacBio HiFi and hybrid inputs explicit.
3. Do not treat ntLink as an independent assembler.
4. Record whether an assembly is collapsed or haplotype-resolved.
5. Never replace unavailable metrics with zero.
6. Preserve exact source-file and tool provenance.
