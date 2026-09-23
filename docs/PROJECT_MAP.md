# LRS benchmarking project map

Updated: 20 September 2026

## Scientific modules

| Stage | Active workflows | Analysis and selected outputs |
|---|---|---|
| Alignment | Root `*.read_mapping.*.smk` files; `final_report_files/snakemake_aligners_benchmarking/` | [`../alignment_analysis/`](../alignment_analysis/README.md) |
| Assembly | [`../assemblers/`](../assemblers/README.md) | [`../assembly_analysis/`](../assembly_analysis/README.md) |
| Small variants | Root `*.snv_calling.*.smk` files and `run_happy.smk` | [`../variant_calling_analysis/`](../variant_calling_analysis/README.md) |
| Structural variants | Root caller workflows; [`../SV aligners call/`](<../SV aligners call/README.md>) | `truvari/`, `happy_results/`, and `variant_calling_analysis/` |
| Report | Analysis-module figures and tables | [`../final_report_files/`](../final_report_files/README.md) |

## Repository mental model

```text
PROTECTED INPUTS
├── fastq/
└── reference/

MAINTAINED WORKFLOWS
├── root alignment and variant .smk entry points
├── assemblers/
├── final_report_files/snakemake_aligners_benchmarking/
└── SV aligners call/

ANALYSIS MODULES
├── alignment_analysis/
├── assembly_analysis/
└── variant_calling_analysis/

NAVIGATION AND PRESENTATION
├── docs/
├── figures/
├── reports/
└── README.md

GENERATED OR HISTORICAL
├── results/, happy_results/, truvari/, vcf_called/
└── archive/, */archive/, */archived_runs/
```

## Current alignment benchmark

The canonical 30× analysis compares HG002, HG003, and HG004 across ONT and PacBio HiFi using minimap2, pbmm2, VACMap, and VG Giraffe.

- final table: [`../alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv`](../alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv)
- final figures: [`../alignment_analysis/figures/30x/final/`](../alignment_analysis/figures/30x/final/)
- final scripts: [`../alignment_analysis/scripts/30x/plots/final/`](../alignment_analysis/scripts/30x/plots/final/)

## Navigation

Use [`FILE_INDEX.md`](FILE_INDEX.md) when looking for a specific important file and [`REPOSITORY_TREE.md`](REPOSITORY_TREE.md) when deciding where new material belongs. Optional shell shortcuts are defined in [`lrs_shortcuts.sh`](lrs_shortcuts.sh).

## Path-stability rule

Do not move an active workflow or protected data directory solely for visual tidiness. First audit references, update configurable paths, run a Snakemake dry run, perform the move, and run the dry run again.
