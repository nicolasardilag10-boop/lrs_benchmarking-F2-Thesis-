# 30x alignment figures

## Final figure set

Use the numbered figures in `final/` for reports and presentations:

1. `01_alignment_error_rate_30x` — alignment error rate by technology,
   aligner, and GIAB sample.
2. `02_aligned_base_yield_30x` — aligned and unaligned base yield.
3. `03_read_mapping_yield_30x` — read retention and mapping yield.
4. `04_mq0_reads_30x` — fraction of mapped reads with MAPQ 0.
5. `05_runtime_30x` — observed wall-clock runtime for minimap2, pbmm2, VACMap,
   and VG Giraffe, with ONT and PacBio HiFi panels and thread allocations.
6. `06_cigar_composition_30x` — mapped, mismatch, and unaligned CIGAR-derived
   fractions.
7. `07_alignment_record_types_30x` — supplementary and secondary alignment
   records per 100 mapped reads.
8. `08_memory_30x` — vertically stacked three-panel memory summary with
   horizontal bars: configured limits, measured ONT peak RSS, and measured
   PacBio HiFi peak RSS.

PNG files are convenient for slides and quick inspection. PDF files are the
preferred vector format for LaTeX and publication output. Figure 3 also has an
SVG version.

Figure 5 plots all 24 runtime values recorded in the canonical benchmark table.
Configuration provenance is retained in the corresponding derived plotting
table instead of being encoded as a separate visual style.

Figure 8 separates configured memory limits into panel A and measured peak RSS
into panels B–C. The quantities are not combined because a memory allocation
limit is not an observed peak-memory measurement. GIAB sample colors and bar
sizes are consistent across panels, while hatching identifies panel A limits.

## Other directories

- `exploratory/` — alternative views, diagnostics, and figures not selected
  for the main report.
- `archive/` — drafts and duplicate/orphaned files kept for provenance.

The current LaTeX report contains no active `\includegraphics` commands, so no
report references required path changes during this reorganization.
