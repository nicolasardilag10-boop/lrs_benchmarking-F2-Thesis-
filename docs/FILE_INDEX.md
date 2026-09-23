# Canonical file index

This page is the shortest route to the files most likely to be edited, executed, inspected, or cited. Historical copies and exploratory outputs are deliberately excluded.

## Alignment benchmark

### Canonical data and figures

| Purpose | File or directory |
|---|---|
| Final 30× table | [`../alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv`](../alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv) |
| Table guide | [`../alignment_analysis/tables/30x/README.md`](../alignment_analysis/tables/30x/README.md) |
| Final figures | [`../alignment_analysis/figures/30x/final/`](../alignment_analysis/figures/30x/final/) |
| Figure guide | [`../alignment_analysis/figures/30x/README.md`](../alignment_analysis/figures/30x/README.md) |
| Derived plotting tables | [`../alignment_analysis/tables/30x/derived/plot_data/`](../alignment_analysis/tables/30x/derived/plot_data/) |
| Raw 30× CRAM statistics | [`../alignment_analysis/statistics_cram_files/`](../alignment_analysis/statistics_cram_files/) |

### Production analysis code

| Purpose | File or directory |
|---|---|
| Extract full alignment metrics | [`../alignment_analysis/scripts/30x/pipeline/quality_check_aligners.py`](../alignment_analysis/scripts/30x/pipeline/quality_check_aligners.py) |
| Extract indel metrics | [`../alignment_analysis/scripts/30x/pipeline/quality_check_aligners_indels.py`](../alignment_analysis/scripts/30x/pipeline/quality_check_aligners_indels.py) |
| Build consolidated table | [`../alignment_analysis/scripts/30x/pipeline/build_alignment_summary_30x_table.py`](../alignment_analysis/scripts/30x/pipeline/build_alignment_summary_30x_table.py) |
| Merge indels | [`../alignment_analysis/scripts/30x/pipeline/merge_indels_into_benchmark.py`](../alignment_analysis/scripts/30x/pipeline/merge_indels_into_benchmark.py) |
| Final plotting scripts | [`../alignment_analysis/scripts/30x/plots/final/`](../alignment_analysis/scripts/30x/plots/final/) |
| Script guide | [`../alignment_analysis/scripts/30x/README.md`](../alignment_analysis/scripts/30x/README.md) |

### Figure-to-script map

| Figure | Script |
|---|---|
| `01_alignment_error_rate_30x` | [`plot_alignment_error_rate_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_alignment_error_rate_30x.py) |
| `02_aligned_base_yield_30x` | [`plot_aligned_base_yield_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_aligned_base_yield_30x.py) |
| `03_read_mapping_yield_30x` | [`plot_read_mapping_yield_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_read_mapping_yield_30x.py) |
| `04_mq0_reads_30x` | [`plot_mq0_reads_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_mq0_reads_30x.py) |
| `05_runtime_30x` (all four aligners; ONT left, PacBio HiFi right) | [`plot_runtime_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_runtime_30x.py) |
| `06_cigar_composition_30x` | [`plot_cigar_composition_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_cigar_composition_30x.py) |
| `07_alignment_record_types_30x` | [`plot_alignment_record_types_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_alignment_record_types_30x.py) |
| `08_memory_30x` (vertically stacked memory panels) | [`plot_memory_30x.py`](../alignment_analysis/scripts/30x/plots/final/plot_memory_30x.py) |

### Alignment workflows

The constitution-required mapper entry points remain at repository root:

```text
ont.read_mapping.{minimap2,pbmm2,vacmap,vg}.smk
pb.read_mapping.{minimap2,pbmm2,vacmap,vg}.smk
```

The integrated alternative workflow is:

- [`../final_report_files/snakemake_aligners_benchmarking/snakemake_aligners.smk`](../final_report_files/snakemake_aligners_benchmarking/snakemake_aligners.smk)
- [`../final_report_files/snakemake_aligners_benchmarking/config.yaml`](../final_report_files/snakemake_aligners_benchmarking/config.yaml)
- [`../final_report_files/snakemake_aligners_benchmarking/samples.tsv`](../final_report_files/snakemake_aligners_benchmarking/samples.tsv)

## Assembly benchmark

| Purpose | Location |
|---|---|
| Assembly overview | [`../assemblers/README.md`](../assemblers/README.md) |
| Whole-genome workflows | [`../assemblers/whole_genome_asm/`](../assemblers/whole_genome_asm/) |
| Chromosome-21 benchmark | [`../assemblers/benchmark_chr21_real/`](../assemblers/benchmark_chr21_real/) |
| Assembly assessment | [`../assemblers/whole_genome_asm/assessment/`](../assemblers/whole_genome_asm/assessment/) |
| Analysis and figures | [`../assembly_analysis/`](../assembly_analysis/README.md) |

## Variant benchmark

| Purpose | Location |
|---|---|
| Variant-analysis module | [`../variant_calling_analysis/`](../variant_calling_analysis/README.md) |
| SV validation workspace | [`../SV aligners call/`](<../SV aligners call/README.md>) |
| Caller inventory | [`inventories/caller_tools_inventory.csv`](inventories/caller_tools_inventory.csv) |
| hap.py outputs | [`../happy_results/`](../happy_results/) |
| Truvari outputs | [`../truvari/`](../truvari/) |
| Visual reports | [`../vcf_called/snv_indel/`](../vcf_called/snv_indel/) |

## Reports and methods

| Purpose | Location |
|---|---|
| Main report source | [`../final_report_files/final_f2_report.tex`](../final_report_files/final_f2_report.tex) |
| Compiled report | [`../final_report_files/final_f2_report.pdf`](../final_report_files/final_f2_report.pdf) |
| Report workspace guide | [`../final_report_files/README.md`](../final_report_files/README.md) |
| Alignment methodology guide | [`methodology/alignment_study_guide/`](methodology/alignment_study_guide/README.md) |
| Tool-selection records | [`methodology/tool_selection/`](methodology/tool_selection/) |

## Maintenance utilities

| Purpose | File |
|---|---|
| Refresh README explorers | [`../scripts/update_repository_tree.py`](../scripts/update_repository_tree.py) |
| Audit cleanup candidates | [`../scripts/repository_cleanup_audit.py`](../scripts/repository_cleanup_audit.py) |
| Install repository hooks | [`../scripts/install_git_hooks.sh`](../scripts/install_git_hooks.sh) |

Files below `archive/`, `exploratory/`, `archived_runs/`, or with historical snapshot names are retained for provenance and should not be used as current entry points.
