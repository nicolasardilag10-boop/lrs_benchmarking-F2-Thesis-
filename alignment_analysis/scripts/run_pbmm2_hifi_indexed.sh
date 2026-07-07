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

# This index was created with the CCS preset.
# CCS and HIFI represent the same pbmm2 preset family.
HIFI_INDEX="$PROJECT/reference/GRCh38_GIABv3_pbmm2_ccs.mmi"

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
# Validate programs and inputs
# ------------------------------------------------------------

if ! command -v pbmm2 >/dev/null 2>&1; then
    echo "ERROR: pbmm2 is not available"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available"
    exit 1
fi

if [ ! -s "$HIFI_INDEX" ]; then
    echo "ERROR: CCS/HIFI index is missing or empty:"
    echo "$HIFI_INDEX"
    exit 1
fi

if [ ! -s "$FASTQ_GZ" ]; then
    echo "ERROR: FASTQ file is missing or empty:"
    echo "$FASTQ_GZ"
    exit 1
fi


# ------------------------------------------------------------
# Determine sample name
# ------------------------------------------------------------

FILENAME=$(basename "$FASTQ_GZ")
STEM="${FILENAME%.fastq.gz}"

case "$FILENAME" in
    *.ont.*|*.pb.*)
        ;;
    *)
        echo "ERROR: expected an ONT or PacBio FASTQ file"
        exit 1
        ;;
esac


# ------------------------------------------------------------
# Output names
# ------------------------------------------------------------

# Keep the existing pbmm2-ccs filename convention.
# In the final benchmark this represents pbmm2-pb.
BAM="$ALIGN_DIR/${STEM}.pbmm2-ccs.sorted.bam"
BAI="$BAM.bai"

FLAGSTAT="$TABLE_DIR/${STEM}.pbmm2-ccs.flagstat.txt"
STATS="$TABLE_DIR/${STEM}.pbmm2-ccs.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/${STEM}.pbmm2-ccs.idxstats.tsv"

LOG="$LOG_DIR/${STEM}.pbmm2-ccs.log"
TMP_DIR="$TMP_ROOT/${STEM}.pbmm2-ccs"


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
# Run pbmm2 HIFI alignment
# ------------------------------------------------------------

echo "Sample:  $STEM"
echo "Reads:   $FASTQ_GZ"
echo "Index:   $HIFI_INDEX"
echo "Output:  $BAM"
echo "Preset:  HIFI"
echo
echo "Starting pbmm2 HIFI alignment..."

TMPDIR="$TMP_DIR" pbmm2 align \
    "$HIFI_INDEX" \
    "$FASTQ_GZ" \
    "$BAM" \
    --preset HIFI \
    --sort \
    --unmapped \
    --rg "@RG\tID:${STEM}\tSM:${STEM}" \
    -j 2 \
    -J 1 \
    2> "$LOG"


# ------------------------------------------------------------
# Validate and index BAM
# ------------------------------------------------------------

samtools quickcheck -v "$BAM"

if [ ! -s "$BAI" ]; then
    samtools index -@ 2 "$BAM"
fi


# ------------------------------------------------------------
# Generate reports
# ------------------------------------------------------------

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"
samtools stats -@ 2 "$BAM" > "$STATS"
samtools idxstats "$BAM" > "$IDXSTATS"


# ------------------------------------------------------------
# Completion summary
# ------------------------------------------------------------

PRIMARY_MAPPED=$(samtools view -c -F 2308 "$BAM")

echo
echo "SUCCESS: pbmm2 HIFI alignment completed"
echo "Primary mapped reads: $PRIMARY_MAPPED"
echo

ls -lh "$BAM" "$BAI"

echo
echo "Mapping summary:"
cat "$FLAGSTAT"
