#!/usr/bin/env bash

set -euo pipefail

# ------------------------------------------------------------
# Usage
# ------------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path/to/sample.fastq.gz"
    exit 1
fi


# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

PROJECT="$HOME/lrs_benchmarking"

SUBREAD_INDEX="$PROJECT/reference/GRCh38_GIABv3_pbmm2_subread.mmi"
FASTQ_GZ="$1"

ALIGN_DIR="$PROJECT/alignment_analysis/alignments"
TABLE_DIR="$PROJECT/alignment_analysis/tables"
LOG_DIR="$PROJECT/alignment_analysis/logs"
TMP_ROOT="$PROJECT/alignment_analysis/tmp"

mkdir -p \
    "$ALIGN_DIR" \
    "$TABLE_DIR" \
    "$LOG_DIR" \
    "$TMP_ROOT"


# ------------------------------------------------------------
# Validate programs and input files
# ------------------------------------------------------------

if ! command -v pbmm2 >/dev/null 2>&1; then
    echo "ERROR: pbmm2 is not available"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available"
    exit 1
fi

if [ ! -s "$SUBREAD_INDEX" ]; then
    echo "ERROR: SUBREAD index is missing or empty:"
    echo "$SUBREAD_INDEX"
    exit 1
fi

if [ ! -s "$FASTQ_GZ" ]; then
    echo "ERROR: FASTQ file is missing or empty:"
    echo "$FASTQ_GZ"
    exit 1
fi


# ------------------------------------------------------------
# Determine the sample name
# ------------------------------------------------------------

FILENAME=$(basename "$FASTQ_GZ")
STEM="${FILENAME%.fastq.gz}"

case "$FILENAME" in
    *.ont.*|*.pb.*)
        ;;
    *)
        echo "ERROR: unrecognised FASTQ filename:"
        echo "$FILENAME"
        exit 1
        ;;
esac


# ------------------------------------------------------------
# Output paths
# ------------------------------------------------------------

BAM="$ALIGN_DIR/${STEM}.pbmm2-subread.sorted.bam"
BAI="$BAM.bai"

FLAGSTAT="$TABLE_DIR/${STEM}.pbmm2-subread.flagstat.txt"
STATS="$TABLE_DIR/${STEM}.pbmm2-subread.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/${STEM}.pbmm2-subread.idxstats.tsv"

LOG="$LOG_DIR/${STEM}.pbmm2-subread.log"
TMP_DIR="$TMP_ROOT/${STEM}.pbmm2-subread"


# ------------------------------------------------------------
# Prepare temporary directory
# ------------------------------------------------------------

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

trap 'rm -rf "$TMP_DIR"' EXIT

rm -f \
    "$BAM" \
    "$BAI" \
    "$FLAGSTAT" \
    "$STATS" \
    "$IDXSTATS" \
    "$LOG"


# ------------------------------------------------------------
# Run pbmm2 SUBREAD alignment
# ------------------------------------------------------------

echo "Sample: $STEM"
echo "Reads:  $FASTQ_GZ"
echo "Index:  $SUBREAD_INDEX"
echo "Output: $BAM"
echo
echo "Starting pbmm2 SUBREAD alignment..."

TMPDIR="$TMP_DIR" pbmm2 align \
    "$SUBREAD_INDEX" \
    "$FASTQ_GZ" \
    "$BAM" \
    --preset SUBREAD \
    --sort \
    --unmapped \
    -j 2 \
    -J 1 \
    2> "$LOG"


# ------------------------------------------------------------
# Validate and index the BAM
# ------------------------------------------------------------

samtools quickcheck -v "$BAM"

if [ ! -s "$BAI" ]; then
    samtools index -@ 2 "$BAM"
fi


# ------------------------------------------------------------
# Generate samtools reports
# ------------------------------------------------------------

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"
samtools stats -@ 2 "$BAM" > "$STATS"
samtools idxstats "$BAM" > "$IDXSTATS"


# ------------------------------------------------------------
# Print completion summary
# ------------------------------------------------------------

PRIMARY_MAPPED=$(samtools view -c -F 2308 "$BAM")

echo
echo "SUCCESS: pbmm2 SUBREAD alignment completed"
echo "Primary mapped reads: $PRIMARY_MAPPED"
echo

ls -lh "$BAM" "$BAI"

echo
echo "Mapping summary:"
cat "$FLAGSTAT"
