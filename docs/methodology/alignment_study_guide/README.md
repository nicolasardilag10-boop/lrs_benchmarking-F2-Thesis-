# Alignment benchmarking study guide

> Conceptual and educational notes for the long-read alignment benchmark. For current executable paths and final outputs, use the [canonical file index](../../FILE_INDEX.md).

## Project at a glance

This project compares long-read alignment behaviour across:

- **Three samples:** HG002, HG003, HG004
- **Two sequencing technologies:** Oxford Nanopore and PacBio HiFi
- **Four aligners:** minimap2, pbmm2, VACMap, and VG Giraffe
- **Eight platform-specific configurations:** ONT and PacBio/HiFi modes for each aligner

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

| Purpose | Current path |
|---|---|
| Final benchmark table | `alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv` |
| Production pipeline | `alignment_analysis/scripts/30x/pipeline/` |
| Final plot scripts | `alignment_analysis/scripts/30x/plots/final/` |
| Final figures | `alignment_analysis/figures/30x/final/` |

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
alignment_benchmark_30x.tsv
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
