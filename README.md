# Long-read sequencing benchmarking

Reproducible comparison of long-read alignment, assembly, and variant-calling workflows using Genome in a Bottle (GIAB) samples. The principal alignment comparison evaluates **minimap2, pbmm2, VACMap, and VG Giraffe** on **Oxford Nanopore (ONT)** and **PacBio HiFi** data.

> **Start here:** [`docs/FILE_INDEX.md`](docs/FILE_INDEX.md) identifies the canonical workflows, tables, figures, scripts, and report files. [`docs/REPOSITORY_TREE.md`](docs/REPOSITORY_TREE.md) explains the directory structure without listing thousands of generated files.

## Most important project outputs

| Item | Canonical location | Description |
|---|---|---|
| Alignment benchmark table | [`alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv`](alignment_analysis/tables/30x/final/alignment_benchmark_30x.tsv) | One row per sample, technology, and aligner configuration |
| Final alignment figures | [`alignment_analysis/figures/30x/final/`](alignment_analysis/figures/30x/final/) | Numbered, publication-ready PNG and PDF figures |
| Alignment plotting code | [`alignment_analysis/scripts/30x/plots/final/`](alignment_analysis/scripts/30x/plots/final/) | One clearly named script per final figure |
| Alignment metric pipeline | [`alignment_analysis/scripts/30x/pipeline/`](alignment_analysis/scripts/30x/pipeline/) | Metric extraction, indel integration, and table construction |
| Assembly workflows | [`assemblers/`](assemblers/README.md) | Whole-genome and chromosome-21 assembly benchmarking |
| Assembly analysis | [`assembly_analysis/`](assembly_analysis/README.md) | Metric aggregation and scientific figures |
| Variant analysis | [`variant_calling_analysis/`](variant_calling_analysis/README.md) | Analysis layer for SNV, indel, and SV benchmarking |
| Main report | [`final_report_files/final_f2_report.tex`](final_report_files/final_f2_report.tex) | LaTeX source; the compiled PDF is stored beside it |

## Study design

The primary 30× alignment benchmark contains:

- GIAB samples HG002, HG003, and HG004;
- ONT and PacBio HiFi sequencing data;
- minimap2, pbmm2, VACMap, and VG Giraffe;
- 24 expected sample–technology–aligner combinations.

The project evaluates mapping yield, alignment error, CIGAR-derived base composition, mapping quality, supplementary and secondary records, indels, runtime, CPU, and memory. Metric definitions and interpretation cautions are documented in [`alignment_analysis/README.md`](alignment_analysis/README.md).

## Scientific workflow

```text
FASTQ + reference
       │
       ├── alignment ──► CRAM/BAM ──► QC and run metrics ──► tables ──► figures
       │
       ├── assembly ───► contigs/scaffolds ──► QUAST/BUSCO/Merqury ──► figures
       │
       └── variants ───► VCF ──► hap.py/Truvari validation ──► figures
```

Detailed maps:

- [`docs/WORKFLOW_MAP.md`](docs/WORKFLOW_MAP.md)
- [`docs/DATA_MAP.md`](docs/DATA_MAP.md)
- [`docs/PROJECT_MAP.md`](docs/PROJECT_MAP.md)

## Running workflows

The project constitution requires Snakemake to run from the repository root:

```bash
cd /path/to/lrs_benchmarking
snakemake --snakefile WORKFLOW.smk --dry-run --printshellcmds
snakemake --snakefile WORKFLOW.smk --cores 8 --printshellcmds
```

Always inspect the relevant module README and configuration before a production run. Active mapper workflows intentionally remain at the repository root because their names and relative execution model are defined by [`CONSTITUTION.md`](CONSTITUTION.md).

## Root workflow index

| Stage | Workflow files |
|---|---|
| ONT alignment | `ont.read_mapping.{minimap2,pbmm2,vacmap,vg}.smk` |
| PacBio alignment | `pb.read_mapping.{minimap2,pbmm2,vacmap,vg}.smk` |
| Small variants | `*.snv_calling.*.smk`, `ont.snv_filtering.v3.smk`, `run_happy.smk` |
| Structural variants | `cuteSV.hg38.smk`, `pbsv.hg38.smk`, `sawfish.hg38.smk`, `sniffles2.hg38.smk`, `truvari_anno.smk` |

## Repository conventions

- `README.md` files define the canonical entry points for each scientific module.
- `scripts/` contains reusable code; `pipeline/` is production data processing and `plots/final/` creates report-ready figures.
- `tables/.../final/` contains canonical analysis tables; `derived/` contains reproducible intermediate tables; `source/` preserves imported source data.
- `figures/.../final/` contains the selected figure set; `exploratory/` contains alternatives; `archive/` contains superseded drafts.
- Generated data remain separate from source workflows. Historical snapshots are never treated as executable entry points.
- Paths in maintained Python code should be repository-relative or discovered from the script location, not hardcoded to a user home directory.

## Repository explorer

<!-- AUTO_REPOSITORY_TREE_START -->
Generated from version-control candidates. This view highlights scientific entry points instead of listing every result file. See [`docs/REPOSITORY_TREE.md`](docs/REPOSITORY_TREE.md) for the expanded map.

| Area | Location | Purpose |
|---|---|---|
| Alignment benchmark | [`alignment_analysis/`](alignment_analysis/README.md) | metrics, canonical tables, and final figures |
| Assembly workflows | [`assemblers/`](assemblers/README.md) | whole-genome and chromosome-21 assembly workflows |
| Assembly analysis | [`assembly_analysis/`](assembly_analysis/README.md) | assembly metric aggregation and figures |
| Variant analysis | [`variant_calling_analysis/`](variant_calling_analysis/README.md) | variant-calling analysis layer |
| SV workspace | [`SV aligners call/`](<SV aligners call/README.md>) | SV and aligner validation workflows |
| Final report | [`final_report_files/`](final_report_files/README.md) | LaTeX report and integrated aligner workflow |
| Documentation | [`docs/`](docs/README.md) | maps, methods, inventories, and troubleshooting |
| Project figures | [`figures/`](figures/README.md) | cross-project and variant figures |
| Results | [`results/`](results/README.md) | organized links to generated outputs |
| Tests | [`tests/`](tests/README.md) | fixtures, smoke tests, and validation |
| Archive | [`archive/`](archive/README.md) | historical code, snapshots, and retired material |

<details>
<summary><b>Constitution-required root workflows</b></summary>

- **Alignment:** [`ont.read_mapping.minimap2.smk`](ont.read_mapping.minimap2.smk), [`ont.read_mapping.pbmm2.smk`](ont.read_mapping.pbmm2.smk), [`ont.read_mapping.vacmap.smk`](ont.read_mapping.vacmap.smk), [`ont.read_mapping.vg.smk`](ont.read_mapping.vg.smk), [`pb.read_mapping.minimap2.smk`](pb.read_mapping.minimap2.smk), [`pb.read_mapping.pbmm2.smk`](pb.read_mapping.pbmm2.smk), [`pb.read_mapping.vacmap.smk`](pb.read_mapping.vacmap.smk), [`pb.read_mapping.vg.smk`](pb.read_mapping.vg.smk)
- **Small variants:** [`ilmn.snv_calling.clair3.smk`](ilmn.snv_calling.clair3.smk), [`ont.snv_calling.clair3.smk`](ont.snv_calling.clair3.smk), [`ont.snv_calling.deepvariant.smk`](ont.snv_calling.deepvariant.smk), [`ont.snv_filtering.v3.smk`](ont.snv_filtering.v3.smk), [`pb.snv_calling.clair3.smk`](pb.snv_calling.clair3.smk), [`pb.snv_calling.deepvariant.smk`](pb.snv_calling.deepvariant.smk), [`run_happy.smk`](run_happy.smk), [`wgs.snv_calling.deepvariant.smk`](wgs.snv_calling.deepvariant.smk)
- **Structural variants:** [`cuteSV.hg38.smk`](cuteSV.hg38.smk), [`pbsv.hg38.smk`](pbsv.hg38.smk), [`sawfish.hg38.smk`](sawfish.hg38.smk), [`sniffles2.hg38.smk`](sniffles2.hg38.smk), [`truvari_anno.smk`](truvari_anno.smk)

</details>

<!-- AUTO_REPOSITORY_TREE_END -->

## Reproducibility and governance

- Project rules and naming requirements: [`CONSTITUTION.md`](CONSTITUTION.md)
- Detailed file index: [`docs/FILE_INDEX.md`](docs/FILE_INDEX.md)
- Repository layout: [`docs/REPOSITORY_TREE.md`](docs/REPOSITORY_TREE.md)
- Organization changelog: [`docs/ORGANIZATION_CHANGELOG.md`](docs/ORGANIZATION_CHANGELOG.md)
- Environment and path guidance: [`docs/f2/PATH_STABILITY.md`](docs/f2/PATH_STABILITY.md)
- Navigation shortcuts: [`docs/lrs_shortcuts.sh`](docs/lrs_shortcuts.sh)

Before committing documentation changes, refresh the generated explorers:

```bash
python3 scripts/update_repository_tree.py
python3 scripts/update_repository_tree.py --check
```
