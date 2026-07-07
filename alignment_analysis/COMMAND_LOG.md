# Alignment analysis command log

## FASTQ read-count analysis

### Objective

Count the number of reads in each FASTQ file.

### Input files

- HG002 ONT
- HG002 PacBio
- HG003 ONT
- HG003 PacBio
- HG004 ONT
- HG004 PacBio

### Command

```bash
SOURCE="/mnt/c/Users/nicol/OneDrive/Documentos/Maestría/Third Semester/F2 Practicals/Material to study/giab(Samples)"

find "$SOURCE" \
-maxdepth 1 \
-type f \
-iname "*.fastq*" \
-print
```

### Result

Each FASTQ subset contains 1,000 reads.

Output table:

```text
alignment_analysis/tables/fastq_read_counts.tsv
```

### Interpretation

The six files are 1,000-read subsets intended for testing and learning.
They are not the complete sequencing datasets.