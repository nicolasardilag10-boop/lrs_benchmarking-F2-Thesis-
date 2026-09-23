# 30x analysis scripts

The scripts are grouped by purpose so the report-facing workflow is easy to
find without mixing it with exploratory work.

## Most relevant plotting scripts

The scripts in `plots/final/` produce the numbered files in
`alignment_analysis/figures/30x/final/`:

1. `plot_alignment_error_rate_30x.py`
2. `plot_aligned_base_yield_30x.py`
3. `plot_read_mapping_yield_30x.py`
4. `plot_mq0_reads_30x.py`
5. `plot_runtime_30x.py`
6. `plot_cigar_composition_30x.py`
7. `plot_alignment_record_types_30x.py`
8. `plot_memory_30x.py`

The final scripts use the canonical table
`tables/30x/final/alignment_benchmark_30x.tsv`. The runtime figure includes all
four aligners and explicitly displays the allocated thread count. All 24
sample, technology, and aligner combinations recorded in the canonical table
are plotted; configuration provenance remains available in the plot-data TSV.
The memory figure uses three vertically stacked panels with horizontal bars:
configured RAM limits, measured ONT peak RSS, and measured PacBio HiFi peak
RSS. Separating the panels prevents allocation limits from being mistaken for
observed usage. All panels retain the GIAB sample colors; hatching identifies
the configured-limit bars in panel A.

Run the new alignment-record figure from the repository root:

```bash
python alignment_analysis/scripts/30x/plots/final/plot_alignment_record_types_30x.py
```

## Other directories

- `pipeline/` — table construction, metric recovery, and merge utilities.
- `plots/exploratory/` — alternative visualizations and diagnostic analyses.
- `tests/fixtures/` — small test inputs.
- `archive/` — superseded or misnamed scripts retained for provenance.

Every active script finds the project root from its own location. No username
or fixed `/home/...` path is required.
