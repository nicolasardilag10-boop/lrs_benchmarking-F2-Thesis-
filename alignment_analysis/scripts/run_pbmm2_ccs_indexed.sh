#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 sample.pb.1k.fastq.gz"
    exit 1
fi

PROJECT="$HOME/lrs_benchmarking"

CCS_INDEX="$PROJECT/reference/GRCh38_GIABv3_pbmm2_ccs.mmi"
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

if ! command -v pbmm2 >/dev/null 2>&1; then
    echo "ERROR: pbmm2 is not available"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available"
    exit 1
fi

if [ ! -s "$CCS_INDEX" ]; then
    echo "ERROR: CCS index is missing or empty:"
    echo "$CCS_INDEX"
    exit 1
fi

if [ ! -s "$FASTQ_GZ" ]; then
    echo "ERROR: FASTQ is missing or empty:"
    echo "$FASTQ_GZ"
    exit 1
fi

FILENAME=$(basename "$FASTQ_GZ")
STEM="${FILENAME%.fastq.gz}"

case "$FILENAME" in
    *.pb.*)
        ;;
    *)
        echo "ERROR: expected a PacBio FASTQ file"
        exit 1
        ;;
esac

BAM="$ALIGN_DIR/${STEM}.pbmm2-ccs.sorted.bam"
BAI="$BAM.bai"

FLAGSTAT="$TABLE_DIR/${STEM}.pbmm2-ccs.flagstat.txt"
STATS="$TABLE_DIR/${STEM}.pbmm2-ccs.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/${STEM}.pbmm2-ccs.idxstats.tsv"

LOG="$LOG_DIR/${STEM}.pbmm2-ccs.log"
TMP_DIR="$TMP_ROOT/${STEM}.pbmm2-ccs"

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

echo "Sample: $STEM"
echo "Index:  $CCS_INDEX"
echo "Reads:  $FASTQ_GZ"
echo "Output: $BAM"
echo
echo "Starting indexed pbmm2 CCS alignment..."

TMPDIR="$TMP_DIR" pbmm2 align \
    "$CCS_INDEX" \
    "$FASTQ_GZ" \
    "$BAM" \
    --preset CCS \
    --sort \
    -j 2 \
    -J 1 \
    2> "$LOG"

samtools quickcheck -v "$BAM"

if [ ! -s "$BAI" ]; then
    samtools index -@ 2 "$BAM"
fi

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"
samtools stats -@ 2 "$BAM" > "$STATS"
samtools idxstats "$BAM" > "$IDXSTATS"

PRIMARY_MAPPED=$(samtools view -c -F 2308 "$BAM")

echo
echo "SUCCESS: indexed pbmm2 CCS alignment completed"
echo "Primary mapped reads: $PRIMARY_MAPPED"
echo

ls -lh "$BAM" "$BAI"
