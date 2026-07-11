# LRS Benchmarking Study Guide

> A structured, reproducible guide for the long-read alignment benchmark using HG002, HG003, and HG004 ONT and PacBio datasets.

## Project at a glance

This project compares long-read alignment behaviour across:

- **Three samples:** HG002, HG003, HG004
- **Two sequencing technologies:** Oxford Nanopore and PacBio HiFi
- **Two aligners:** minimap2 and pbmm2
- **Four tested configurations:** `mm2-ont`, `mm2-pb`, `pbmm2-ont`, `pbmm2-pb`

The workflow moves from raw FASTQ files to quality control, alignment, summary tables, figures, statistics, interpretation, and version control.

## Documentation map

| Document | Purpose |
|---|---|
| [Project workflow](PROJECT_WORKFLOW.md) | Complete colored Mermaid flowchart and project stages |
| [File input/output map](FILE_INPUT_OUTPUT_MAP.md) | Shows how each file type becomes the next |
| [Reusable code patterns](REUSABLE_CODE_PATTERNS.md) | Important Python blocks repeated throughout the project |
| [Interpretation guide](ANALYSIS_INTERPRETATION.md) | Scientific wording and cautions for the results |
| [Reproduction checklist](REPRODUCTION_CHECKLIST.md) | Exact commands and validation steps |

## Core project files

```text
alignment_analysis/tables/alignment_summary.tsv
alignment_analysis/scripts/build_alignment_summary.py
alignment_analysis/scripts/plot_general_technology_comparison.py
alignment_analysis/scripts/plot_final_alignment_benchmark.py
alignment_analysis/scripts/plot_mapped_bases_cigar_final.py
```

## Main analytical logic

```text
FASTQ.gz
   ↓
FASTQ quality control
   ↓
SAM/BAM alignments
   ↓
samtools reports
   ↓
alignment_summary.tsv
   ↓
Pandas filtering and validation
   ↓
Matplotlib figures
   ↓
paired statistics
   ↓
scientific interpretation
```

!!! tip "Recommended study order"
    Read the documents in this order: workflow → file map → code patterns → interpretation → reproduction checklist.