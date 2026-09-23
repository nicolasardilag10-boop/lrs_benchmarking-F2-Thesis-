# Variant-Calling Benchmark Analysis

This directory is reserved for the variant-calling benchmarking component of the F2 project.

The wider repository already contains variant-calling workflows created or modified by multiple contributors. The repository reorganization **does not move, rename or claim ownership of those existing workflows**.

This module is the location for the F2 comparison layer as it is developed:

```text
variant_calling_analysis/
├── README.md
├── scripts/
│   ├── metrics/
│   └── plots/
├── tables/
└── figures/
```

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. Directories remain compact; use this module's curated sections for canonical files.

- [`README.md`](README.md)

- **figures/** — 1 file
- **scripts/** — 2 files
- **tables/** — 1 file

<!-- AUTO_REPOSITORY_TREE_END -->

## Intended responsibilities

- collect caller outputs for a defined benchmark design
- extract comparable metrics
- preserve caller / aligner / sample / technology provenance
- benchmark against appropriate truth sets and confident regions
- create standardized summary tables
- create figures from the canonical tables

## Existing caller workflows

Existing root-level or project-level caller `.smk` files remain in their current locations unless a dedicated, reviewed refactor is performed later.

This prevents the F2 organizational cleanup from breaking workflows that belong to the wider repository.

## Future README expansion

When the final caller set is fixed, document:

- callers included in the comparison
- exact versions and Docker images
- required aligned inputs
- truth sets / confident regions
- evaluation tools and metric definitions
- canonical output table
- figure scripts
- known caller-specific caveats
