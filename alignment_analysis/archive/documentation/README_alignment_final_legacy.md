# Long-Read Alignment Summary

This directory contains the script used to build one comparable alignment
summary table for Oxford Nanopore (ONT) and PacBio HiFi reads.

The benchmark includes the same core metrics for:

| Read technology | Aligners |
|---|---|
| ONT | minimap2, pbmm2, VACMap, VG Giraffe |
| PacBio HiFi | minimap2, pbmm2, VACMap, VG Giraffe |

The script **does not run the aligners** and **does not modify CRAM files**.
It reads existing `samtools stats` reports and combines them into one TSV.

Missing metrics are written as `NA`. `NA` must not be interpreted as zero.

---

## 1. Repository layout

Recommended repository structure:

```text
lrs_benchmarking/
├── cram/
│   ├── *.cram
│   ├── *.cram.crai
│   ├── *.cram.idxstats
│   └── *.cram.stats
│
└── alignment_analysis/
    ├── scripts/
    │   └── build_alignment_summary.py
    ├── tables/
    │   ├── alignment_summary.tsv
    │   ├── alignment_metadata.tsv        # optional
    │   └── *.samtools_stats.txt          # optional copied/test stats
    └── README.md
```

On the UKW server, the repository root is expected to be:

```text
/data/genmedbfx/schilling_m/repos/lrs_benchmarking
```

Large CRAM, CRAI and FASTQ files should remain on the server and should not be
committed to GitHub.

---

## 2. Canonical dataset and filename convention

The repository Constitution defines the dataset as:

```text
<sample>.<platform>.<coverage>
```

Examples:

```text
HG002.ont.30x
HG002.pb.30x
HG003.ont.1k
```

The canonical mapping-result filename contains:

```text
<dataset>.<reference>.<mapper_tag>
```

For this project the reference is `hg38`.

Examples:

```text
HG002.ont.30x.hg38.mm2-ont.cram.stats
HG002.pb.30x.hg38.pbmm2-pb.cram.stats
HG003.ont.30x.hg38.vacmap-ont.cram.stats
HG004.pb.30x.hg38.vg-pb.cram.stats
```

### Important: `1k` and `30x` are different datasets

The `1k` files are downsampled test datasets. They must **not** replace the
production `30x` files.

The script therefore creates one row per:

```text
sample × read technology × coverage × reference × mapper
```

For example, these are kept as separate rows:

```text
HG002.ont.1k.hg38.mm2-ont
HG002.ont.30x.hg38.mm2-ont
```

---

## 3. Supported mapper tags

Canonical mapper tags follow the repository Constitution:

| Mapper | ONT | PacBio HiFi |
|---|---|---|
| minimap2 | `mm2-ont` | `mm2-pb` |
| pbmm2 | `pbmm2-ont` | `pbmm2-pb` |
| VACMap | `vacmap-ont` | `vacmap-pb` |
| VG Giraffe | `vg-ont` | `vg-pb` |

The corresponding presets/configurations are:

| Mapper tag | Aligner | Preset |
|---|---|---|
| `mm2-ont` | minimap2 | `map-ont` |
| `mm2-pb` | minimap2 | `map-hifi` |
| `pbmm2-ont` | pbmm2 | `SUBREAD` |
| `pbmm2-pb` | pbmm2 | `CCS/HIFI` |
| `vacmap-ont` | VACMap | `vacmap-ont` |
| `vacmap-pb` | VACMap | `vacmap-pb` |
| `vg-ont` | VG Giraffe | `vg-ont` |
| `vg-pb` | VG Giraffe | `vg-pb` |

For compatibility with older files, the parser also accepts:

```text
pbmm2-subread  -> pbmm2-ont
pbmm2-ccs      -> pbmm2-pb
```

Older statistics filenames that omit `.hg38.` are also accepted because this
benchmark is fixed to the `hg38` reference. New files should use the canonical
naming convention.

---

## 4. Input files

The script searches recursively in:

```text
alignment_analysis/tables/
cram/
```

It recognizes:

```text
*.samtools_stats.txt
*.cram.stats
*.stats.txt
*.stats
```

The filename must contain:

- `HG002`, `HG003` or `HG004`
- `.ont.` or `.pb.`
- a coverage token such as `.1k.` or `.30x.`
- one supported mapper tag

---

## 5. What is extracted from `samtools stats`

The script reads three standard sections of `samtools stats`.

### `SN` summary section

Used for:

```text
raw total sequences
reads mapped
reads unmapped
reads MQ0
non-primary alignments
supplementary alignments
total length
bases mapped
bases mapped (cigar)
mismatches
error rate
average length
maximum length
```

`non-primary alignments` is reported in the TSV as
`secondary_alignments`, because it corresponds to secondary alignments.

### `MAPQ` section

Used to calculate:

```text
mapq_mean
mapq_median
```

MAPQ value `255` means mapping quality is unavailable and is excluded from the
numerical mean and median.

MAPQ should be used descriptively. MAPQ values are produced by each mapper and
are not necessarily calibrated identically between different aligners.

### `ID` indel-distribution section

Used to calculate:

```text
insertion_events
deletion_events
inserted_bases
deleted_bases
insertion_events_per_100kb
deletion_events_per_100kb
```

The per-100-kb values use `bases mapped (cigar)` as the denominator.

Official reference:

```text
https://www.htslib.org/doc/samtools-stats.html
```

---

## 6. Output table

The script writes:

```text
alignment_analysis/tables/alignment_summary.tsv
```

### Identification columns

```text
sample
read_technology
coverage
dataset
reference
aligner
preset
mapper_tag
configuration
statistics_file
```

### Mapping and alignment-quality columns

```text
raw_total_sequences
reads_mapped
reads_unmapped
mapped_reads_percent

reads_mq0
reads_mq0_percent
mapq_mean
mapq_median

secondary_alignments
secondary_alignments_per_100_mapped_reads
supplementary_alignments
supplementary_alignments_per_100_mapped_reads

total_length
bases_mapped
bases_mapped_cigar
mapped_bases_percent

mismatches
error_rate
error_percent

insertion_events
deletion_events
inserted_bases
deleted_bases
insertion_events_per_100kb
deletion_events_per_100kb

average_length
maximum_length
```

### Optional enrichment columns

The same columns are present for every aligner:

```text
soft_clipped_bases
soft_clipped_percent

mean_coverage
median_coverage
breadth_1x_percent
breadth_10x_percent
breadth_20x_percent
breadth_30x_percent

runtime_seconds
peak_ram_mb
threads

aligner_version
samtools_version
aligner_docker_image
samtools_docker_image
```

These optional values are **not guessed**. If they are not available, the
script writes `NA`.

---

## 7. Optional metadata file

Metrics that are not reliably available from the standard `samtools stats`
summary can be supplied through:

```text
alignment_analysis/tables/alignment_metadata.tsv
```

The metadata file is optional.

It must use these five columns to identify a row:

```text
sample
read_technology
coverage
reference
mapper_tag
```

It may then contain any of the optional enrichment columns.

Example:

```text
sample	read_technology	coverage	reference	mapper_tag	runtime_seconds	peak_ram_mb	threads	aligner_version	samtools_version
HG002	ONT	30x	hg38	mm2-ont	1234.5	8192	32	2.27	1.21
HG002	PacBio	30x	hg38	mm2-pb	1102.2	7900	32	2.27	1.21
```

You do **not** need to create this file before running the summary script.
Without it, the core `samtools stats` metrics are still extracted normally.

The metadata file is the appropriate place to add:

- coverage metrics from the selected coverage workflow/tool
- soft-clipping metrics if calculated separately
- Snakemake or `/usr/bin/time` runtime/RAM measurements
- thread counts
- exact tool versions
- pinned Docker image names

This keeps the main parser simple and prevents values from being inferred
incorrectly.

---

## 8. File-selection behavior

If two statistics files describe the exact same:

```text
sample × technology × coverage × reference × mapper
```

the script keeps one deterministic source.

Preference is:

1. canonical `cram/*.cram.stats`
2. another `*.cram.stats`
3. `*.samtools_stats.txt`
4. `*.stats.txt`
5. another `*.stats`

This preference is used **only for exact duplicate rows**.

A `1k` file never replaces a `30x` file.

Files with an unrecognized sample, technology or mapper tag are skipped.

A warning is printed if a selected file does not contain
`raw total sequences` or `reads mapped`.

---

## 9. Running on the UKW server

From the repository root:

```bash
cd /data/genmedbfx/schilling_m/repos/lrs_benchmarking

python3 alignment_analysis/scripts/build_alignment_summary.py \
    --project /data/genmedbfx/schilling_m/repos/lrs_benchmarking \
    --out /data/genmedbfx/schilling_m/repos/lrs_benchmarking/alignment_analysis/tables/alignment_summary.tsv
```

If the script is stored at:

```text
alignment_analysis/scripts/build_alignment_summary.py
```

the repository root can also be inferred automatically:

```bash
python3 alignment_analysis/scripts/build_alignment_summary.py
```

The script itself only reads text statistics/metadata and writes the summary
TSV. It does not execute minimap2, pbmm2, VACMap, VG, samtools or any other
benchmarking binary.

All mapper/QC workflows should continue to follow the repository Constitution,
including pinned Docker images and the existing containerization rules.

---

## 10. Checking the output

Basic checks:

```bash
test -s alignment_analysis/tables/alignment_summary.tsv

head -n 2 alignment_analysis/tables/alignment_summary.tsv

column -t -s $'\t' alignment_analysis/tables/alignment_summary.tsv | less -S
```

Check row count and detected aligners:

```bash
python3 - <<'PY'
import csv

path = "alignment_analysis/tables/alignment_summary.tsv"

with open(path, newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

print(f"Rows: {len(rows)}")
print("Aligners:", sorted({row["aligner"] for row in rows}))
print("Coverages:", sorted({row["coverage"] for row in rows}))

for row in rows:
    print(
        row["sample"],
        row["read_technology"],
        row["coverage"],
        row["aligner"],
        row["mapper_tag"],
    )
PY
```

For a complete `30x` production comparison with three samples, two
technologies and four aligners, the theoretical maximum is:

```text
3 samples × 2 technologies × 4 aligners = 24 rows
```

The actual number may be lower if some alignment/statistics files have not yet
been produced.

The `1k` test rows are additional rows and are not counted as part of those 24
production combinations.

---

## 11. Interpretation

Use the metrics together rather than ranking an aligner from one number.

- `mapped_reads_percent` describes how many reads received an alignment.
- `mapped_bases_percent` describes how much input read sequence was represented
  in mapped alignments.
- `reads_mq0_percent` identifies mapped reads with MAPQ 0.
- `mapq_mean` and `mapq_median` describe mapper-reported mapping confidence,
  but should not be treated as directly calibrated accuracy scores across
  different aligners.
- `secondary_alignments` and `supplementary_alignments` describe additional
  alignment structure.
- `mismatches` and `error_rate` are the values reported by `samtools stats`;
  the mismatch statistic is NM-derived and should not be described as a
  substitution-only error count.
- insertion/deletion metrics describe indel burden in the alignments.
- coverage metrics describe how completely and deeply the reference is covered.
- runtime, RAM and threads describe computational efficiency and must be
  interpreted together.

All aligners should be compared using the same:

```text
sample
read technology
coverage
reference
input/filtering policy
metric definitions
resource-accounting method
```

A high mapping percentage alone does not prove that an aligner is the most
accurate.

The complete TSV should be retained as the source table for downstream figures,
statistical tests and the final report.
