# Integrated aligner benchmark workflow

This is the integrated Snakemake implementation for minimap2, pbmm2, VACMap, and VG Giraffe.

## Active entry points

| File | Purpose |
|---|---|
| [`snakemake_aligners.smk`](snakemake_aligners.smk) | Main workflow |
| [`config.yaml`](config.yaml) | Production configuration |
| [`samples.tsv`](samples.tsv) | Sample sheet |
| [`envs/`](envs/) | Conda environment definitions |
| [`scripts/`](scripts/) | Current metric, summary, and statistical scripts |
| [`run_dryrun_only_test.sh`](run_dryrun_only_test.sh) | Dry-run validation helper |

Run from this directory only if the configuration expects it; otherwise follow the repository-root execution rule documented in [`../../README.md`](../../README.md). Always perform a dry run before production execution.

## Non-production material

- `config.dryrun.yaml` and `snakemake_aligners.dryrun.smk`: dry-run validation inputs.
- `config.vg_test_working.yaml`: VG-focused test configuration.
- `final_validation_snapshot/`: frozen validation snapshot.
- `archived_runs/`: historical runs and code snapshots; do not edit or execute as the current workflow.
- `.snakemake/`, logs, and `results/`: generated runtime state and outputs.

Downstream curated 30× tables and figures are in [`../../alignment_analysis/`](../../alignment_analysis/README.md).
