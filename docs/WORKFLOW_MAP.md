# Workflow map

Updated: 20 September 2026

## Alignment benchmarking

```text
HG002 / HG003 / HG004
ONT / PacBio HiFi
        │
        ├── minimap2
        ├── pbmm2
        ├── VACMap
        └── VG Giraffe
              │
              ▼
        BAM or CRAM alignments
              │
              ├── samtools statistics and flag summaries
              └── wall time, CPU, and peak RAM
              │
              ▼
        canonical 30× benchmark table
              │
              ▼
        statistical comparisons + numbered final figures
```

Workflow entry points are the constitution-required root mapper files and the integrated workflow at `final_report_files/snakemake_aligners_benchmarking/snakemake_aligners.smk`. Curated analysis is documented in [`../alignment_analysis/README.md`](../alignment_analysis/README.md).

## Assembly benchmarking

```text
ONT and/or PacBio HiFi reads
        │
        ├── Flye
        ├── GoldRush ──► ntLink scaffolding
        └── Verkko
              │
              ▼
        assemblies
              │
              ├── QUAST
              ├── BUSCO
              └── Merqury
              │
              ▼
        assembly_analysis/
```

Production and reduced chromosome-21 workflows are documented in [`../assemblers/README.md`](../assemblers/README.md).

## Variant benchmarking

```text
aligned reads
     │
     ├── Clair3 / DeepVariant ──► SNV and indel VCFs ──► hap.py
     └── Sniffles2 / cuteSV / pbsv / Sawfish ──► SV VCFs ──► Truvari
                                                       │
                                                       ▼
                                             variant_calling_analysis/
```

Root caller workflows remain in place for path stability. Validation-specific aligner rules are under [`../SV aligners call/`](<../SV aligners call/README.md>).

## Report flow

```text
canonical final tables
        │
        ▼
versioned plotting scripts
        │
        ▼
numbered PDF + PNG figures
        │
        ▼
final_report_files/final_f2_report.tex
```

The report currently contains no active figure inclusions. When figures are added, reference the canonical PDF in the appropriate analysis module and record the relationship in its README.
