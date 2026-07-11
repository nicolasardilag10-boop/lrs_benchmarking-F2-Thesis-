# File Input/Output Map

## Compact lineage

```text
FASTQ.gz
   ↓
FASTQ QC tables
   ↓
SAM/BAM alignments
   ↓
flagstat / idxstats / samtools stats
   ↓
alignment_summary.tsv
   ↓
Pandas filtering, reshaping and validation
   ↓
PNG / PDF figures
   ↓
statistics and interpretation
```

## Stage-by-stage table

| Stage | Input | Operation | Output | Validation |
|---|---|---|---|---|
| Raw data | `*.fastq.gz` | Collect six long-read subsets | Six FASTQ files | Naming and compression |
| Read counts | FASTQ | Count records | `fastq_read_counts.tsv` | 1,000 reads per file |
| Length QC | FASTQ | Median, mean, maximum, N50 | TSV summaries | No missing/negative lengths |
| Quality QC | FASTQ qualities | Q20 and Q30 | TSV summaries | Correct technology labels |
| Reference | FASTA | Prepare reference/index | Reference/index files | Index generated successfully |
| Test alignment | FASTQ + reference | Trial run | SAM/BAM | Pipeline works before scaling |
| Full benchmark | Six FASTQs + reference | Run four configurations | 24 alignments | Complete 3 × 2 × 4 design |
| Alignment reports | SAM/BAM | `flagstat`, `idxstats`, `samtools stats` | TXT/TSV reports | Every alignment has reports |
| Master summary | Reports | Parse and combine | `alignment_summary.tsv` | 24 rows, no duplicates |
| Plot preparation | Master table | Filter, map, pivot | DataFrames | Complete paired values |
| Plotting | DataFrames | Matplotlib | PNG/PDF | Correct units and labels |
| Statistics | Paired arrays | `ttest_rel` | t and p values | Same samples and order |
| Reporting | Figures + tests | Interpret | Slides/docs | Avoid unsupported claims |

## Main file types

### `.fastq.gz`

Compressed sequencing reads.

### `.fa` / `.fasta`

Reference genome sequence.

### `.sam` / `.bam`

Read-to-reference alignments.

### `.flagstat.txt`

Overall mapped/unmapped read statistics.

### `.idxstats.tsv`

Reference-level alignment counts.

### `.samtools_stats.txt`

Detailed alignment statistics.

### `.tsv`

Structured tab-separated analysis tables.

### `.png` / `.pdf`

Final figures for review, slides, and reporting.

## Core columns in `alignment_summary.tsv`

```text
sample
read_technology
aligner
preset
configuration
raw_total_sequences
reads_mapped
reads_unmapped
mapped_reads_percent
total_length
bases_mapped
bases_mapped_cigar
mapped_bases_percent
mismatches
error_rate
error_percent
insertions
deletions
average_length
maximum_length
non_primary_alignments
```