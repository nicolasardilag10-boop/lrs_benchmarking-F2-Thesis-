# 30x alignment tables

This directory is organized by the role of each table. Start with the file in
`final/` for routine analysis.

## Canonical table

- `final/alignment_benchmark_30x.tsv` — primary 30x benchmarking table. It has
  24 rows and combines alignment statistics with the best available runtime,
  CPU, RAM, and indel-event measurements. In this table, `insertions` and
  `deletions` mean the number of insertion and deletion events calculated from
  samtools `ID` records; they are not inserted/deleted base totals.

View its performance columns in a terminal:

```bash
cut -f1-5,30-34,40-41 final/alignment_benchmark_30x.tsv \
  | column -t -s $'\t' \
  | less -S
```

View the indel event columns:

```bash
cut -f1-5,18-19 final/alignment_benchmark_30x.tsv \
  | column -t -s $'\t' \
  | less -S
```

## Directory guide

- `final/` — approved tables intended for reporting and routine analysis.
- `source/` — source and active intermediate tables consumed by scripts.
- `derived/` — tables generated for plots, diagnostics, and statistical tests.
- `archive/` — superseded tables and backups retained for provenance.

## Important supporting tables

- `source/alignment_summary_30x_Samtools_Christian.tsv` — core samtools
  alignment statistics used by many plotting scripts.
- `source/mm2.run_metrics.tsv` — raw minimap2 runtime records.
- `source/alignment_summary_30x_combined_metrics_without_ram.tsv` — input to
  the runtime/RAM merge script.
- `source/alignment_summary_30x_combined_trusted_metrics.tsv` — trusted metric
  subset used by the MAPQ analysis.
- `derived/alignment_metrics_30x.tsv` — wider alignment-quality table,
  including indel and coverage fields where available.
- `derived/alignment_benchmark_30x_indel_recovery.tsv` — detailed recovery
  table containing event counts, affected-base counts, and normalized rates.
- `derived/alignment_benchmark_30x_indel_recovery_report.tsv` — provenance and
  recovery status for each benchmark row.

## Notes

- Files in `archive/` are not current analysis inputs.
- A filename containing `.backup` or `.bak` is historical and should not be
  used for new analysis.
- Measured peak RAM is not available for every aligner in the canonical table;
  `command_ram_limit_gb` is a configured limit, not measured usage.
- Runtime values are populated for all 24 combinations. The
  `runtime_configuration` and `wallclock_runtime_source_log` columns preserve
  the provenance of each value and are exported with the runtime plot data.
- Indel event counts are available for 17 of 24 configurations. Six ONT
  minimap2/pbmm2 configurations lack matching full stats files, and the HG004
  ONT VG stats file contains no `ID` records.

After regenerating the canonical runtime table, restore its recovered indel
event columns with:

```bash
python alignment_analysis/scripts/30x/pipeline/merge_indels_into_benchmark.py
```
