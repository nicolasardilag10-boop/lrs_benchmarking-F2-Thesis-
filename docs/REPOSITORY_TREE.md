# Repository tree and responsibilities

Updated: 20 September 2026

This is a curated tree. It shows where maintained code and primary outputs belong while intentionally omitting thousands of generated BAM/CRAM/VCF, validation, cache, and benchmark-detail files.

```text
lrs_benchmarking/
├── README.md                       # project landing page
├── CONSTITUTION.md                 # binding workflow and naming rules
├── *.read_mapping.*.smk            # active mapper entry points (root by policy)
├── *.snv_calling.*.smk             # small-variant workflow entry points
├── {cuteSV,pbsv,sawfish,sniffles2}.hg38.smk
│                                     # structural-variant entry points
│
├── alignment_analysis/             # alignment metrics, tables, and figures
│   ├── scripts/30x/
│   │   ├── pipeline/               # production table-building code
│   │   ├── plots/final/            # scripts for selected figures
│   │   ├── plots/exploratory/      # alternatives and diagnostics
│   │   └── tests/fixtures/         # small test data
│   ├── tables/30x/
│   │   ├── final/                  # canonical alignment table
│   │   ├── derived/                # reproducible secondary tables
│   │   ├── source/                 # imported/source tables
│   │   └── archive/                # superseded versions
│   └── figures/30x/
│       ├── final/                  # numbered publication figures
│       ├── exploratory/            # non-final analyses
│       └── archive/                # superseded drafts
│
├── assemblers/                     # active assembly workflows
│   ├── whole_genome_asm/           # production whole-genome rules
│   ├── benchmark_chr21_real/       # reduced chromosome-21 benchmark
│   ├── containers/                 # reproducible container definitions
│   └── scripts/                    # assembly utilities
├── assembly_analysis/              # assembly tables and figures
├── variant_calling_analysis/       # variant metric/plot analysis layer
├── SV aligners call/               # SV and aligner validation workspace
│
├── final_report_files/
│   ├── final_f2_report.tex         # main manuscript/report source
│   ├── final_f2_report.pdf         # compiled report
│   └── snakemake_aligners_benchmarking/
│       ├── snakemake_aligners.smk  # integrated alignment workflow
│       ├── config.yaml             # active configuration
│       ├── samples.tsv             # active sample sheet
│       ├── scripts/                # current analysis helpers
│       ├── envs/                   # reproducible environments
│       └── archived_runs/          # historical code and runs
│
├── docs/                           # maintained documentation and indexes
│   ├── FILE_INDEX.md               # canonical-file lookup table
│   ├── REPOSITORY_TREE.md          # this document
│   ├── PROJECT_MAP.md              # project modules
│   ├── WORKFLOW_MAP.md             # scientific data flow
│   ├── DATA_MAP.md                 # source/generated data map
│   ├── inventories/                # tool inventories
│   └── methodology/                # methods and study guides
├── scripts/                        # repository maintenance utilities
├── tests/                          # fixtures, smoke tests, validation
├── figures/                        # cross-project figures
├── reports/                        # organized report-facing exports
├── results/                        # generated-output links and groupings
├── data/                           # stable data entry points/symlinks
└── archive/                        # retired code, snapshots, and duplicates
```

## Directory status

| Status | Meaning | Examples |
|---|---|---|
| Canonical | Preferred source or result for current work | `alignment_analysis/tables/30x/final/`, root mapper workflows |
| Active | Maintained workflow or ongoing project area | `assemblers/`, `SV aligners call/` |
| Generated | Re-creatable or result-heavy data | `results/`, `happy_results/`, `truvari/` |
| Exploratory | Scientifically useful but not selected as final | `*/exploratory/` |
| Historical | Retained only for provenance | `archive/`, `archived_runs/` |
| External | Third-party checkout/tool; not project source | `overleaf-toolkit/` |

## Why some workflows remain at root

The root-level `.smk` files are intentional. Their naming and execution model are specified by the project constitution, and several use repository-root-relative paths. Moving them only for visual symmetry would reduce reproducibility. The README and file index provide organization without breaking those stable entry points.
