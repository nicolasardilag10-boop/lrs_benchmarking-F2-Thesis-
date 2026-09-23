# Final report workspace

This directory contains the F2 report and workflow material used to support it.

## Report

| File | Role |
|---|---|
| [`final_f2_report.tex`](final_f2_report.tex) | Canonical LaTeX source |
| [`commands_latex_files.tex`](commands_latex_files.tex) | Supporting LaTeX command/reference material |
| [`final_f2_report.pdf`](final_f2_report.pdf) | Compiled report |

LaTeX auxiliary files (`.aux`, `.fls`, `.fdb_latexmk`, `.out`, `.toc`) are generated and ignored by Git.

## Workflows

- [`snakemake_aligners_benchmarking/`](snakemake_aligners_benchmarking/README.md): integrated alignment workflow, active configuration, sample sheet, environments, and validation.
- `snakemake_assemblers_benchmarking/`: historical/working assembly workflow tree. The maintained assembly entry point is documented under [`../assemblers/`](../assemblers/README.md).
- Historical Snakemake training projects have been moved to [`../archive/tutorials/`](../archive/tutorials/) so they cannot be mistaken for production workflows.

The report currently has no active `\includegraphics` statements; figure inclusion should use paths documented in the relevant analysis module.
