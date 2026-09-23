# Chromosome 21 assembly benchmark contract

## Samples

- HG002
- HG003
- HG004

## Technologies

- PacBio HiFi
- Oxford Nanopore

## Standard input policy

- Region: chromosome 21
- Primary alignments only
- Minimum read quality: Q20
- Target depth: 30x relative to the same chr21 reference sequence
- Deterministic downsampling with one recorded seed per dataset
- Whole-read selection
- FASTQ output for assembly
- No silent overwriting

## Required input validation

Each normalized dataset must pass:

1. Source BAM exists and is nonempty.
2. samtools quickcheck passes.
3. Reference header contains chromosome 21.
4. Primary read count is greater than zero.
5. Read names are unique in the normalized FASTQ.
6. FASTQ structure is valid.
7. Total sequence bases are recorded.
8. Observed coverage is within the accepted tolerance.
9. SHA256 checksum is recorded.
10. Validation status is PASS.

## Benchmark groups

### Single-technology de novo assembly

- Flye ONT
- Flye HiFi
- GoldRush ONT
- GoldRush HiFi

### Hybrid assembly

- Verkko using paired HiFi and ONT inputs

### Scaffolding

- Flye assemblies before and after ntLink

## Exclusion policy

A condition marked FAILED must not enter summary statistics, paired
tests or performance figures. It must be reported as FAILED or NA,
not as zero performance.
