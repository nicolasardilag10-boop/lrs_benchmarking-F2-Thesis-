# Alignment Analysis Scripts

This directory contains scripts supporting the long-read alignment benchmark.

## Canonical production script

```text
30x/pipeline/quality_check_aligners.py
```

This is the canonical metric-extraction script for the 30x alignment benchmark.

## Recommended organization

```text
scripts/
├── 30x/           # pipeline, final plots, exploratory plots, and tests
├── plots/         # figure-generation scripts
├── validation/    # targeted test / validation runners
├── utils/         # FASTQ and small helper utilities
├── legacy/        # retained historical scripts, not canonical
└── archive/       # timestamped development snapshots
```

The 30x scripts use project-root discovery based on the repository structure,
so they remain portable after being grouped into subdirectories. See
[`30x/README.md`](30x/README.md) for the current layout.

Files directly in `scripts/` belong to the earlier 1k analysis and remain in
place for compatibility. They must not be confused with the current 30x
pipeline. Timestamped `before_*` copies have been moved to
`archive/snapshots/` so only runnable legacy entry points remain visible.
