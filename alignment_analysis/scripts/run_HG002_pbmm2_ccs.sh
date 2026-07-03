#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking"

REF="$PROJECT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

READ_GZ="$PROJECT/samples_try/HG002.pb.1k.fastq.gz"

ALIGN_DIR="$PROJECT/alignment_analysis/alignments"
TABLE_DIR="$PROJECT/alignment_analysis/tables"
LOG_DIR="$PROJECT/alignment_analysis/logs"
TMP_DIR="$PROJECT/alignment_analysis/tmp/HG002_pbmm2_ccs"

READ="$TMP_DIR/HG002.pb.1k.fastq"

BAM="$ALIGN_DIR/HG002.pb.1k.pbmm2-ccs.sorted.bam"
LOG="$LOG_DIR/HG002.pb.1k.pbmm2-ccs.log"

FLAGSTAT="$TABLE_DIR/HG002.pb.1k.pbmm2-ccs.flagstat.txt"
STATS="$TABLE_DIR/HG002.pb.1k.pbmm2-ccs.samtools_stats.txt"
IDXSTATS="$TABLE_DIR/HG002.pb.1k.pbmm2-ccs.idxstats.tsv"

mkdir -p \
    "$ALIGN_DIR" \
    "$TABLE_DIR" \
    "$LOG_DIR"

rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

if [ ! -s "$REF" ]; then
    echo "ERROR: reference missing"
    exit 1
fi

if [ ! -s "$READ_GZ" ]; then
    echo "ERROR: PacBio FASTQ missing"
    exit 1
fi

echo "Decompressing the 1,000-read PacBio FASTQ..."
gzip -dc "$READ_GZ" > "$READ"

echo "Running pbmm2 with the CCS preset..."

if ! /usr/bin/time -v pbmm2 align \
    "$REF" \
    "$READ" \
    "$BAM" \
    --preset CCS \
    --sort \
    --rg '@RG\tID:HG002_pbmm2_ccs\tSM:HG002\tPL:PACBIO' \
    -j 2 \
    -J 1 \
    2> "$LOG"
then
    echo "ERROR: pbmm2 alignment failed"
    tail -n 40 "$LOG"
    exit 1
fi

samtools index -@ 2 "$BAM"

samtools quickcheck -v "$BAM"

samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT"
samtools stats -@ 2 "$BAM" > "$STATS"
samtools idxstats "$BAM" > "$IDXSTATS"

echo
echo "SUCCESS: pbmm2 CCS alignment completed"
ls -lh "$BAM" "$BAM.bai"

echo
echo "Mapping summary:"
cat "$FLAGSTAT"

rm -rf "$TMP_DIR"
