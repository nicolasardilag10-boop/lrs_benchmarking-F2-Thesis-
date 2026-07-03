#!/usr/bin/env bash

# Stop on errors, undefined variables, or failed pipelines.
set -euo pipefail

# ------------------------------------------------------------
# 1. Check the command-line argument
# ------------------------------------------------------------

if [ "$#" -ne 1 ]; then
    echo "Usage:"
    echo "  $0 path/to/sample.fastq.gz"
    exit 1
fi


# ------------------------------------------------------------
# 2. Define project paths
# ------------------------------------------------------------

PROJECT="$HOME/lrs_benchmarking"

REF="$PROJECT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

FASTQ="$1"

ALIGN_DIR="$PROJECT/alignment_analysis/alignments"
LOG_DIR="$PROJECT/alignment_analysis/logs"
TABLE_DIR="$PROJECT/alignment_analysis/tables"
TMP_ROOT="$PROJECT/alignment_analysis/tmp"

mkdir -p \
    "$ALIGN_DIR" \
    "$LOG_DIR" \
    "$TABLE_DIR" \
    "$TMP_ROOT"


# ------------------------------------------------------------
# 3. Validate the input files
# ------------------------------------------------------------

if [ ! -s "$REF" ]; then
    echo "ERROR: reference FASTA missing:"
    echo "$REF"
    exit 1
fi

if [ ! -s "$FASTQ" ]; then
    echo "ERROR: FASTQ file missing:"
    echo "$FASTQ"
    exit 1
fi


# ------------------------------------------------------------
# 4. Determine sample name and technology
# ------------------------------------------------------------

FILENAME=$(basename "$FASTQ")
STEM="${FILENAME%.fastq.gz}"

case "$FILENAME" in

    *.ont.*)
        PRESET="map-ont"
        TECHNOLOGY="ONT"
        ;;

    *.pb.*)
        PRESET="map-hifi"
        TECHNOLOGY="PacBio-HiFi"
        ;;

    *)
        echo "ERROR: technology could not be detected from:"
        echo "$FILENAME"
        exit 1
        ;;

esac


# ------------------------------------------------------------
# 5. Define output paths
# ------------------------------------------------------------

BAM="$ALIGN_DIR/${STEM}.sorted.bam"
BAI="$BAM.bai"

LOG="$LOG_DIR/${STEM}.minimap2.log"

FLAGSTAT="$TABLE_DIR/${STEM}.flagstat.txt"
STATS="$TABLE_DIR/${STEM}.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/${STEM}.idxstats.tsv"

TMP="$TMP_ROOT/$STEM"


# ------------------------------------------------------------
# 6. Prepare a clean temporary directory
# ------------------------------------------------------------

rm -rf "$TMP"
mkdir -p "$TMP"

rm -f \
    "$BAM" \
    "$BAI" \
    "$LOG" \
    "$FLAGSTAT" \
    "$STATS" \
    "$IDXSTATS"


# ------------------------------------------------------------
# 7. Display the configuration
# ------------------------------------------------------------

echo "============================================"
echo "Sample:      $STEM"
echo "Technology:  $TECHNOLOGY"
echo "Preset:      $PRESET"
echo "FASTQ:       $FASTQ"
echo "Reference:   $REF"
echo "Output BAM:  $BAM"
echo "============================================"
echo


# ------------------------------------------------------------
# 8. Align and sort directly into BAM
# ------------------------------------------------------------

# -I 1G and --split-prefix reduce reference-index memory use.
#
# The minimap2 SAM output is sent directly into samtools sort.
# This avoids creating another large intermediate SAM file.

/usr/bin/time -v minimap2 \
    -t 2 \
    -I 1G \
    -K 50M \
    -ax "$PRESET" \
    --split-prefix="$TMP/mm2" \
    "$REF" \
    "$FASTQ" \
    2> "$LOG" |
samtools sort \
    -@ 1 \
    -m 512M \
    -T "$TMP/sort" \
    -o "$BAM" \
    -


# ------------------------------------------------------------
# 9. Index and validate the BAM
# ------------------------------------------------------------

samtools index -@ 2 "$BAM"

samtools quickcheck -v "$BAM"

echo "BAM validation: OK"


# ------------------------------------------------------------
# 10. Generate alignment statistics
# ------------------------------------------------------------

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"

samtools stats -@ 2 "$BAM" > "$STATS"

samtools idxstats "$BAM" > "$IDXSTATS"


# ------------------------------------------------------------
# 11. Count primary mapped reads
# ------------------------------------------------------------

PRIMARY_MAPPED=$(samtools view -c -F 2308 "$BAM")


# ------------------------------------------------------------
# 12. Report completion
# ------------------------------------------------------------

echo
echo "Alignment completed successfully"
echo "Primary mapped reads: $PRIMARY_MAPPED"
echo

ls -lh "$BAM" "$BAI"

echo
echo "Reports:"
echo "$FLAGSTAT"
echo "$STATS"
echo "$IDXSTATS"


# ------------------------------------------------------------
# 13. Remove temporary files
# ------------------------------------------------------------

rm -rf "$TMP"
