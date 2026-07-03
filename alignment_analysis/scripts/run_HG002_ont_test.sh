#!/usr/bin/env bash

# Stop when an undefined variable is used.
set -u

PROJECT="$HOME/lrs_benchmarking"

REF="$PROJECT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

READ="$PROJECT/samples_try/HG002.ont.1k.fastq.gz"

ALIGN_DIR="$PROJECT/alignment_analysis/alignments"
LOG_DIR="$PROJECT/alignment_analysis/logs"
TMP_DIR="$PROJECT/alignment_analysis/tmp/HG002_ont"

SAM="$ALIGN_DIR/HG002.ont.test.sam"
LOG="$LOG_DIR/HG002.ont.test.minimap2.log"

# Create all required directories.
mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$TMP_DIR"

# Validate input files.
if [ ! -s "$REF" ]; then
    echo "ERROR: reference FASTA not found:"
    echo "$REF"
    exit 1
fi

if [ ! -s "$READ" ]; then
    echo "ERROR: FASTQ file not found:"
    echo "$READ"
    exit 1
fi

# Check that the output directory is writable.
if ! touch "$ALIGN_DIR/write_test.txt"; then
    echo "ERROR: alignment directory is not writable"
    exit 1
fi

rm -f "$ALIGN_DIR/write_test.txt"

# Remove old or incomplete output.
rm -f "$SAM" "$LOG"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

echo "Reference: $REF"
echo "Reads:     $READ"
echo "SAM:       $SAM"
echo "Log:       $LOG"
echo
echo "Starting minimap2..."

# Run a low-memory ONT alignment.
/usr/bin/time -v minimap2 \
    -t 2 \
    -I 1G \
    -K 50M \
    -ax map-ont \
    --split-prefix="$TMP_DIR/mm2" \
    -o "$SAM" \
    "$REF" \
    "$READ" \
    2> "$LOG"

STATUS=$?

echo
echo "Minimap2 exit code: $STATUS"

if [ "$STATUS" -eq 0 ] && [ -s "$SAM" ]; then
    echo "SUCCESS: SAM file created"
    ls -lh "$SAM"
else
    echo "ERROR: alignment failed"
    echo
    tail -n 30 "$LOG"
    exit 1
fi
