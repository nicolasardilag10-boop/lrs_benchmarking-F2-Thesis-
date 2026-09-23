# Data map

Updated: 20 September 2026

## Data classes

| Class | Location | Policy |
|---|---|---|
| Source reads | [`../fastq/`](../fastq/) | Protected; do not rename or move while workflows use root-relative paths |
| Reference genome/indexes | [`../reference/`](../reference/) | Protected; replace only with provenance and compatible indexes |
| Small test reads | `samples_try/` | Retained for smoke tests and legacy validation |
| Stable data aliases | [`../data/`](../data/README.md) | Symlinks to protected FASTQ and reference locations |
| Alignment statistics | `alignment_analysis/statistics_cram_files/`, `samtools_stats_30x_Christian/` | Source evidence for metric extraction |
| Canonical analysis tables | `*/tables/.../final/` | Preferred inputs for final plots and reporting |
| Derived tables | `*/tables/.../derived/` | Re-creatable analysis intermediates |
| Generated workflow outputs | `results/`, `assemblies/`, `happy_results/`, `truvari/`, `vcf_called/` | Do not edit manually |
| Historical data products | `archive/`, `*/archive/`, `*/archived_runs/` | Provenance only; not current input |

## Alignment data flow

```text
fastq/ + reference/
        │
        ▼
root mapper workflows or integrated aligner workflow
        │
        ▼
CRAM/BAM + samtools statistics + runtime metrics
        │
        ▼
alignment_analysis/scripts/30x/pipeline/
        │
        ▼
alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv
        │
        ▼
alignment_analysis/scripts/30x/plots/final/
        │
        ▼
alignment_analysis/figures/30x/final/
```

## Source versus generated files

A file is a source input only when its provenance is recorded and a workflow consumes it directly. FASTQ-like files inside `results/`, `work/`, `debug/`, or assembler output directories are intermediates unless a module README explicitly says otherwise.

Large sequence files and indexes are intentionally ignored by Git. Their absence from Git status does not mean they are disposable.
