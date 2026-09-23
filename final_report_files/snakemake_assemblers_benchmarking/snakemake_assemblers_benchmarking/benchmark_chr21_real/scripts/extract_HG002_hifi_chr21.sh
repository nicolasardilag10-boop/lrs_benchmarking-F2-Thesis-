#!/usr/bin/env bash

set -euo pipefail


ROOT="benchmark_chr21_real"

BAM_URL="https://downloads.pacbcloud.com/public/dataset/HG002-CpG-methylation-202202/HG002.GRCh38.haplotagged.bam"

OUTPUT_DIR="${ROOT}/inputs/pacbio/HG002"

BAM="${OUTPUT_DIR}/HG002.hifi.chr21.primary.q20.bam"

FASTQ="${OUTPUT_DIR}/HG002.hifi.chr21.primary.q20.fastq.gz"

HEADER="${OUTPUT_DIR}/HG002.remote.header.sam"

IDXSTATS="${OUTPUT_DIR}/HG002.remote.idxstats.tsv"

SUMMARY="${ROOT}/results/validation/HG002.hifi.chr21.summary.tsv"

PROVENANCE="${ROOT}/results/validation/HG002.hifi.chr21.provenance.tsv"

REMOTE_LOG="${ROOT}/results/logs/extraction/HG002.hifi.remote.log"

FASTQ_LOG="${ROOT}/results/logs/extraction/HG002.hifi.fastq.log"

FILTER_SCRIPT="${ROOT}/scripts/filter_sam_progress.py"

STREAM_SCRIPT="${ROOT}/scripts/stream_file_progress.py"

CHR21_LENGTH=46709983

# Remove unmapped, secondary and supplementary alignments.
EXCLUDED_FLAGS=2308

MINIMUM_MAPQ=20


progress() {
    printf '[%3d%%] %s\n' "$1" "$2"
}


fail() {
    echo "ERROR: $1" >&2
    exit 1
}


for program in samtools python3 gzip
do
    command -v "$program" >/dev/null 2>&1 \
        || fail "required program not found: $program"
done


[[ -f "$FILTER_SCRIPT" ]] \
    || fail "missing helper: $FILTER_SCRIPT"

[[ -f "$STREAM_SCRIPT" ]] \
    || fail "missing helper: $STREAM_SCRIPT"


if pgrep -f \
    'flye-modules|minimap2.*benchmark_1k|snakemake.*benchmark_1k' \
    >/dev/null 2>&1
then
    fail "an old benchmark process is still running"
fi


mkdir -p \
    "$OUTPUT_DIR" \
    "$(dirname "$SUMMARY")" \
    "$(dirname "$REMOTE_LOG")"


progress 0 \
    "Starting real HG002 PacBio HiFi chromosome-21 extraction."

progress 3 \
    "Testing the indexed remote BAM and reading its header."


samtools view \
    -H \
    "$BAM_URL" \
    > "$HEADER"


if grep -q $'SN:chr21\t' "$HEADER"
then
    CHROMOSOME="chr21"
elif grep -q $'SN:21\t' "$HEADER"
then
    CHROMOSOME="21"
else
    fail "chromosome 21 was not found in the remote BAM header"
fi


progress 7 \
    "Detected chromosome name: ${CHROMOSOME}."

progress 9 \
    "Reading the remote BAM index."


samtools idxstats \
    "$BAM_URL" \
    > "$IDXSTATS"


TOTAL_RECORDS="$(
    awk \
        -v chromosome="$CHROMOSOME" \
        '$1 == chromosome {print $3}' \
        "$IDXSTATS"
)"


[[ -n "$TOTAL_RECORDS" ]] \
    || fail "no chr21 count was returned by the BAM index"

[[ "$TOTAL_RECORDS" -gt 0 ]] \
    || fail "the BAM index reports zero chr21 records"


progress 12 \
    "Remote index reports ${TOTAL_RECORDS} chr21 records."

TMP_BAM="${BAM}.partial"

rm -f \
    "$TMP_BAM" \
    "$BAM" \
    "${BAM}.bai"


progress 15 \
    "Streaming chr21; retaining primary alignments with MAPQ >= 20."


set +e

samtools view \
    -h \
    "$BAM_URL" \
    "$CHROMOSOME" \
    2> "$REMOTE_LOG" \
| python3 "$FILTER_SCRIPT" \
    "$TOTAL_RECORDS" \
    15 \
    70 \
    "$EXCLUDED_FLAGS" \
    "$MINIMUM_MAPQ" \
| samtools view \
    -@ 2 \
    -b \
    -o "$TMP_BAM" \
    - \
    2>> "$REMOTE_LOG"

PIPE_STATUS=("${PIPESTATUS[@]}")

set -e


for status in "${PIPE_STATUS[@]}"
do
    if [[ "$status" -ne 0 ]]
    then
        fail \
            "chr21 streaming pipeline failed; inspect $REMOTE_LOG"
    fi
done


mv "$TMP_BAM" "$BAM"


progress 72 \
    "Checking the extracted BAM."


samtools quickcheck "$BAM" \
    || fail "the extracted BAM failed samtools quickcheck"


SELECTED_READS="$(
    samtools view -c "$BAM"
)"


[[ "$SELECTED_READS" -gt 0 ]] \
    || fail "no reads remained after filtering"


progress 75 \
    "Indexing local BAM with ${SELECTED_READS} selected reads."


samtools index \
    -@ 2 \
    "$BAM"


TMP_FASTQ="${FASTQ}.partial"

rm -f \
    "$TMP_FASTQ" \
    "$FASTQ"


progress 80 \
    "Converting selected chromosome-21 reads to FASTQ."


set +e

python3 "$STREAM_SCRIPT" \
    "$BAM" \
    80 \
    92 \
| samtools fastq \
    -@ 2 \
    -n \
    - \
    2> "$FASTQ_LOG" \
| gzip -1 \
    > "$TMP_FASTQ"

PIPE_STATUS=("${PIPESTATUS[@]}")

set -e


for status in "${PIPE_STATUS[@]}"
do
    if [[ "$status" -ne 0 ]]
    then
        fail \
            "BAM-to-FASTQ conversion failed; inspect $FASTQ_LOG"
    fi
done


mv "$TMP_FASTQ" "$FASTQ"


progress 94 \
    "Validating the compressed FASTQ."


gzip -t "$FASTQ" \
    || fail "gzip validation failed"


progress 96 \
    "Calculating read lengths and chromosome coverage."


python3 - \
    "$FASTQ" \
    "$BAM" \
    "$SUMMARY" \
    "$CHR21_LENGTH" \
    "$BAM_URL" \
<<'PY'
from __future__ import annotations

import gzip
import statistics
import subprocess
import sys
from pathlib import Path


fastq_path = Path(sys.argv[1])
bam_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
chromosome_length = int(sys.argv[4])
source_url = sys.argv[5]


lengths: list[int] = []
first_header = ""

with gzip.open(
    fastq_path,
    "rt",
    encoding="utf-8",
) as handle:
    while True:
        header = handle.readline()

        if not header:
            break

        sequence = handle.readline().rstrip("\r\n")
        plus = handle.readline()
        quality = handle.readline().rstrip("\r\n")

        if not sequence or not plus or not quality:
            raise RuntimeError(
                f"Truncated FASTQ: {fastq_path}"
            )

        if not header.startswith("@"):
            raise RuntimeError(
                f"Invalid FASTQ header: {fastq_path}"
            )

        if not plus.startswith("+"):
            raise RuntimeError(
                f"Invalid FASTQ separator: {fastq_path}"
            )

        if len(sequence) != len(quality):
            raise RuntimeError(
                f"Sequence-quality mismatch: {fastq_path}"
            )

        if not first_header:
            first_header = header.rstrip("\r\n")

        lengths.append(len(sequence))


def n50(values: list[int]) -> int:
    target = sum(values) / 2
    accumulated = 0

    for value in sorted(values, reverse=True):
        accumulated += value

        if accumulated >= target:
            return value

    return 0


read_count = len(lengths)
total_bases = sum(lengths)

mean_length = (
    statistics.mean(lengths)
    if lengths
    else 0
)

coverage = (
    total_bases / chromosome_length
    if chromosome_length
    else 0
)

coverage_status = (
    "READY_FOR_30X_DOWNSAMPLING"
    if coverage >= 30
    else "BELOW_30X"
)


fields = [
    "sample",
    "technology",
    "chromosome",
    "filter",
    "reads",
    "total_read_bases",
    "mean_read_length",
    "read_n50",
    "minimum_read_length",
    "maximum_read_length",
    "chromosome_length",
    "estimated_input_coverage",
    "target_coverage",
    "coverage_status",
    "first_read_header",
    "source",
    "fastq",
    "bam",
]


values = [
    "HG002",
    "PacBio_HiFi",
    "chr21",
    "primary_MAPQ20",
    str(read_count),
    str(total_bases),
    f"{mean_length:.2f}",
    str(n50(lengths)),
    str(min(lengths) if lengths else 0),
    str(max(lengths) if lengths else 0),
    str(chromosome_length),
    f"{coverage:.4f}",
    "30",
    coverage_status,
    first_header or "NA",
    source_url,
    str(fastq_path),
    str(bam_path),
]


summary_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

with summary_path.open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write("\t".join(fields) + "\n")
    handle.write("\t".join(values) + "\n")


print(f"Reads: {read_count:,}")
print(f"Total read bases: {total_bases:,}")
print(f"Read N50: {n50(lengths):,} bp")
print(f"Estimated chr21 coverage: {coverage:.4f}x")
print(f"Status: {coverage_status}")
PY


cat > "$PROVENANCE" <<EOF
sample	technology	source_dataset	source_bam	chromosome	filter	target_coverage
HG002	PacBio_HiFi	HPRC_HG002_Data_Freeze_v1.0	${BAM_URL}	${CHROMOSOME}	primary_MAPQ20	30x
EOF


progress 99 \
    "Final BAM and FASTQ integrity checks."


samtools quickcheck "$BAM"
gzip -t "$FASTQ"


progress 100 \
    "HG002 HiFi chromosome-21 extraction completed."


echo
echo "=== COVERAGE SUMMARY ==="

column -t -s $'\t' "$SUMMARY"

echo
echo "=== OUTPUT FILES ==="

ls -lh \
    "$BAM" \
    "${BAM}.bai" \
    "$FASTQ" \
    "$SUMMARY" \
    "$PROVENANCE"
