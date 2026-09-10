# Long-Read Alignment Summary

This workflow produces one comparable TSV for the four long-read aligners used
in the project:

- minimap2
- pbmm2
- VACMap
- VG Giraffe

for both:

- Oxford Nanopore (ONT)
- PacBio HiFi

The same columns and definitions are used for every aligner.

## What is automatic

The script uses existing `samtools stats` reports for the standard alignment
metrics and, when an hg38 FASTA is provided, uses the corresponding CRAM files
to calculate the remaining alignment-quality metrics with a pinned samtools
Docker image.

### Automatically obtained from `samtools stats`

- total reads
- mapped reads
- unmapped reads
- mapped-read percentage
- MQ0 reads and MQ0 percentage
- mean MAPQ
- median MAPQ
- secondary alignments
- supplementary alignments
- total read bases
- mapped bases
- CIGAR-mapped bases
- mapped-bases percentage
- mismatches
- error rate / error percentage
- insertion events
- deletion events
- inserted bases
- deleted bases
- insertion events per 100 kb
- deletion events per 100 kb
- average read length
- maximum read length

`samtools stats` 1.24 defines `reads MQ0`, secondary (`non-primary`) and
supplementary alignments in its SN section, MAPQ distributions in the MAPQ
section, and indel-size distributions in the ID section.

### Automatically calculated from each CRAM

When `--reference-fasta` is supplied:

- soft-clipped bases
- soft-clipped percentage
- mean coverage
- median coverage
- breadth >=1x
- breadth >=10x
- breadth >=20x
- breadth >=30x
- samtools version
- samtools Docker image

Coverage is calculated with `samtools depth -aa`, so zero-depth reference
positions are included. For comparability, coverage and clipping exclude
unmapped, secondary, supplementary, QC-failed and duplicate records.

The pinned default utility image is:

```text
quay.io/biocontainers/samtools:1.24--h9dcdb79_1
```

This respects the repository requirement that external bioinformatics tools run
inside pinned Docker images.

## Metrics that must have been recorded during the mapper run

These values describe the original mapper execution:

- runtime_seconds
- peak_ram_mb
- threads
- aligner_version
- aligner_docker_image

Runtime, peak RAM, thread count and the exact mapper image/version **cannot be
reconstructed scientifically after an alignment has already finished if the
workflow never recorded them**.

The script therefore:

1. automatically looks for a uniquely matching Snakemake benchmark file and
   reads `s` as runtime and `max_rss` as peak RAM when available;
2. reads exact run provenance from
   `alignment_analysis/tables/alignment_run_metadata.tsv` when present;
3. writes `NA` rather than inventing any missing value.

This is intentional.

## Canonical datasets

Coverage is part of the dataset identity:

```text
<sample>.<platform>.<coverage>
```

Examples:

```text
HG002.ont.1k
HG002.ont.30x
HG002.pb.30x
```

`1k` and `30x` are never deduplicated against each other.

The final row identity is:

```text
sample × technology × coverage × reference × mapper
```

For the production 30x benchmark:

```text
3 samples × 2 technologies × 4 aligners = 24 possible rows
```

## Canonical mapper tags

| Aligner | ONT | PacBio HiFi |
|---|---|---|
| minimap2 | `mm2-ont` | `mm2-pb` |
| pbmm2 | `pbmm2-ont` | `pbmm2-pb` |
| VACMap | `vacmap-ont` | `vacmap-pb` |
| VG Giraffe | `vg-ont` | `vg-pb` |

Legacy aliases are accepted:

```text
pbmm2-subread -> pbmm2-ont
pbmm2-ccs     -> pbmm2-pb
```

## Canonical filenames

New files should include `hg38`, for example:

```text
HG002.ont.30x.hg38.mm2-ont.cram
HG002.ont.30x.hg38.mm2-ont.cram.stats

HG002.pb.30x.hg38.pbmm2-pb.cram
HG002.pb.30x.hg38.pbmm2-pb.cram.stats
```

Older recognized files without `.hg38.` remain readable for compatibility.

## Files searched

Statistics are searched recursively under:

```text
alignment_analysis/tables/
cram/
samtools_stats_30x_Christian/
```

The repository currently contains committed 30x SN extracts such as:

```text
samtools_stats_30x_Christian/HG002_ont_30x.hg38.mm2-ont.cram.stats.SN.txt
```

These files are sufficient for the SN summary metrics (mapped/unmapped reads,
MQ0, secondary/supplementary counts, mapped bases, mismatches, error rate and
read lengths). When the full server-side `cram/*.cram.stats` file exists, it is
preferred automatically because the full report also contains the MAPQ and ID
sections required for MAPQ-distribution and detailed indel metrics.

Recognized statistics suffixes:

```text
*.samtools_stats.txt
*.cram.stats.SN.txt
*.cram.stats
*.stats.txt
*.stats
```

CRAM files are searched under:

```text
cram/
```

## Quick test in this clean branch

The clean branch may not contain the large `cram/` directory locally. You can
still verify the parser against the committed 30x SN reports:

```bash
cd ~/lrs_benchmarking_clean_pr

python3 alignment_analysis/scripts/30x/quality_check_aligners.py     --project "$PWD"
```

This should create:

```text
alignment_analysis/tables/alignment_summary.tsv
```

Rows built only from `*.cram.stats.SN.txt` will correctly leave metrics that
require the full MAPQ/ID sections or the CRAM itself as `NA`.

## Recommended full server command

From the repository root:

```bash
cd /data/genmedbfx/schilling_m/repos/lrs_benchmarking
```

Run:

```bash
python3 alignment_analysis/scripts/30x/quality_check_aligners.py \
    --project /data/genmedbfx/schilling_m/repos/lrs_benchmarking \
    --reference-fasta /PATH/TO/YOUR/hg38.fa \
    --threads 8
```

Output:

```text
alignment_analysis/tables/alignment_summary.tsv
```

Replace `/PATH/TO/YOUR/hg38.fa` with the exact hg38 FASTA already used by the
mapping workflow. Do not use a different reference for the summary analysis.

If `--reference-fasta` is omitted, the script still builds all metrics that can
be obtained from the existing `samtools stats` reports, but CRAM-derived
coverage/clipping fields remain `NA`.

## Optional run metadata

If the mapper workflows already record runtime/provenance, place the exact
values in:

```text
alignment_analysis/tables/alignment_run_metadata.tsv
```

Minimal header:

```text
sample	read_technology	coverage	reference	mapper_tag	runtime_seconds	peak_ram_mb	threads	aligner_version	aligner_docker_image
```

Example:

```text
HG002	ONT	30x	hg38	mm2-ont	1234.5	8192	32	2.27	your-pinned-minimap2-image:tag
HG002	PacBio	30x	hg38	pbmm2-pb	1400.2	9100	32	1.13.1	your-pinned-pbmm2-image:tag
```

Do not copy the example values into real results. Use the values recorded by
your actual workflow.

## Final output columns

The TSV contains:

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
cram_file

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

soft_clipped_bases
soft_clipped_percent

average_length
maximum_length

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

metrics_complete
missing_metrics
```

## How to know whether everything is present

Two final columns make this explicit:

```text
metrics_complete
missing_metrics
```

If a row has every requested metric:

```text
metrics_complete = YES
missing_metrics   =
```

If something was never recorded, for example mapper RAM:

```text
metrics_complete = NO
missing_metrics   = peak_ram_mb
```

This prevents incomplete rows from looking complete.

## Interpretation notes

Use the same dataset, reference and filtering rules for all four aligners.

MAPQ is useful descriptively, but different aligners may calibrate MAPQ
differently. Do not rank mapper accuracy solely by mean or median MAPQ.

The `mismatches` and `error_rate` fields are those reported by `samtools stats`.
The mismatch count is NM-derived, so it should not be described as a pure
substitution-only count.

Runtime should always be interpreted together with the thread count.

The master TSV should be the single source used for downstream plots,
statistics and the final report.
