#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   run_minimap2_cross_preset.sh reads.fastq.gz mm2-ont
#   run_minimap2_cross_preset.sh reads.fastq.gz mm2-pb

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 path/to/sample.fastq.gz {mm2-ont|mm2-pb}"
    exit 1
fi

PROJECT="$HOME/lrs_benchmarking"

REF="$PROJECT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

FASTQ_GZ="$1"
CONFIGURATION="$2"

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
# Validate tools and inputs
# ------------------------------------------------------------

if ! command -v minimap2 >/dev/null 2>&1; then
    echo "ERROR: minimap2 is not available"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools is not available"
    exit 1
fi

if [ ! -s "$REF" ]; then
    echo "ERROR: reference FASTA is missing:"
    echo "$REF"
    exit 1
fi

if [ ! -s "$FASTQ_GZ" ]; then
    echo "ERROR: FASTQ is missing:"
    echo "$FASTQ_GZ"
    exit 1
fi

# ------------------------------------------------------------
# Translate benchmark label into minimap2 preset
# ------------------------------------------------------------

case "$CONFIGURATION" in
    mm2-ont)
        PRESET="map-ont"
        ;;
    mm2-pb)
        PRESET="map-hifi"
        ;;
    *)
        echo "ERROR: configuration must be mm2-ont or mm2-pb"
        exit 1
        ;;
esac

# ------------------------------------------------------------
# Output filenames
# ------------------------------------------------------------

FILENAME=$(basename "$FASTQ_GZ")
STEM="${FILENAME%.fastq.gz}"

BAM="$ALIGN_DIR/${STEM}.${CONFIGURATION}.sorted.bam"
BAI="$BAM.bai"

FLAGSTAT="$TABLE_DIR/${STEM}.${CONFIGURATION}.flagstat.txt"
STATS="$TABLE_DIR/${STEM}.${CONFIGURATION}.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/${STEM}.${CONFIGURATION}.idxstats.tsv"

LOG="$LOG_DIR/${STEM}.${CONFIGURATION}.log"
TMP_DIR="$TMP_ROOT/${STEM}.${CONFIGURATION}"

SORT_PREFIX="$TMP_DIR/sort"
SPLIT_PREFIX="$TMP_DIR/minimap2-split"

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
# Alignment
# ------------------------------------------------------------

echo "Sample:        $STEM"
echo "Configuration: $CONFIGURATION"
echo "Preset:        $PRESET"
echo "Reads:         $FASTQ_GZ"
echo "Output:        $BAM"
echo
echo "Starting minimap2 alignment..."

minimap2 \
    -t 2 \
    -ax "$PRESET" \
    -I 1G \
    -K 50M \
    --split-prefix="$SPLIT_PREFIX" \
    "$REF" \
    "$FASTQ_GZ" \
    2> "$LOG" |
samtools sort \
    -@ 1 \
    -m 512M \
    -T "$SORT_PREFIX" \
    -o "$BAM" \
    - \
    2>> "$LOG"

# ------------------------------------------------------------
# Validation and reports
# ------------------------------------------------------------

samtools quickcheck -v "$BAM"

samtools index -@ 2 "$BAM"

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"
samtools stats -@ 2 "$BAM" > "$STATS"
samtools idxstats "$BAM" > "$IDXSTATS"

PRIMARY_MAPPED=$(samtools view -c -F 2308 "$BAM")

echo
echo "SUCCESS: minimap2 cross-preset alignment completed"
echo "Primary mapped reads: $PRIMARY_MAPPED"
echo

ls -lh "$BAM" "$BAI"

echo
echo "Mapping summary:"
cat "$FLAGSTAT"
