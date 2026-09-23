# Repository organization changelog

## 20 September 2026

This pass improved discoverability while preserving active workflow paths and scientific results.

### Canonical navigation added

- `docs/FILE_INDEX.md`: lookup table for important workflows, scripts, tables, figures, and reports.
- `docs/REPOSITORY_TREE.md`: curated repository tree and directory responsibilities.
- concise READMEs for documentation, results, reports, tests, scripts, data aliases, external validation areas, and historical outputs.
- generated README explorers changed from exhaustive file dumps to compact scientific entry points.

### Files moved from repository root

| Previous location | New location |
|---|---|
| `Workflow_and_code_implemented/` | `docs/methodology/alignment_study_guide/` |
| `caller_tools_inventory.csv` | `docs/inventories/caller_tools_inventory.csv` |
| `lrs_benchmarking_tools_selection.R.Rmd` | `docs/methodology/tool_selection/lrs_benchmarking_tools_selection.Rmd` |
| `lrs_benchmarking_tools_selection.xlsx` | `docs/methodology/tool_selection/lrs_benchmarking_tools_selection.xlsx` |
| `sv_f1_barplots.png` | `figures/variants/sv_f1_barplots.png` |
| repository-organization shell scripts | `archive/legacy/organizer_scripts/` |
| empty typo placeholders | `archive/unclassified/empty_placeholders/` |
| duplicate test table | `archive/duplicates/test_alignment_metrics.root-copy.tsv` |

### Historical material separated from active code

- Aligner `before_*` and `backup_*` files moved to `final_report_files/snakemake_aligners_benchmarking/archived_runs/code_snapshots/`.
- Old alignment-analysis snapshots moved to `alignment_analysis/scripts/archive/snapshots/` and `alignment_analysis/tables/archive/snapshots/`.
- Legacy alignment documentation moved to `alignment_analysis/archive/documentation/`.
- Duplicate navigation maps moved from `docs_maps/` to `archive/legacy/docs_maps/`.
- Snakemake tutorial projects moved from `final_report_files/` to `archive/tutorials/`.
- Manual alignment test outputs moved to `alignment_analysis/tests/manual_runs/`.

### Metadata cleanup

Forty-six untracked Windows `Zone.Identifier` sidecar files were removed. These contained operating-system download metadata, not scientific data. Python bytecode caches generated during validation were also removed.

### Paths deliberately retained

- root-level mapper and variant-caller `.smk` files, because the constitution and active relative paths require them;
- `fastq/` and `reference/`, because they are protected workflow inputs;
- established generated-result directories such as `happy_results/`, `truvari/`, `vcf_called/`, and `sawfish/`;
- `SV aligners call/`, despite the space in its name, to avoid breaking active validation paths;
- the report and integrated aligner workflow under `final_report_files/`.

### Validation performed

- generated README explorer check;
- local-link validation across maintained Markdown files;
- Python compilation for current 30× analysis and repository utilities;
- shell syntax check for `docs/lrs_shortcuts.sh`;
- execution of all seven final 30× plotting scripts;
- Git whitespace/error check.
