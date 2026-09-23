#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"

ORIGINAL_REFERENCE="$PROJECT/data/reference/genome.fasta"

READS="$PROJECT/results/vacmap_diagnostic/hg002_hifi_input_test/HG002.pb.first20.fastq"

OUT="$PROJECT/results/vacmap_diagnostic/hg002_hifi_fresh_whole_genome"

FRESH_REFERENCE="$OUT/reference_fresh.fasta"
BAM="$OUT/HG002.pb.first20.fresh_index.vacmap.sorted.bam"
LOG="$OUT/HG002.pb.first20.fresh_index.vacmap.log"
SUMMARY="$OUT/test_summary.tsv"

GLOBAL_START=$(date +%s)

# Approximate estimates:
# index creation: about 30 minutes
# mapping 20 reads: about 7 minutes
EXPECTED_INDEX_SECONDS=1800
EXPECTED_MAPPING_SECONDS=420

format_time() {
    local seconds="${1:-0}"

    if (( seconds < 0 )); then
        seconds=0
    fi

    printf '%02d:%02d:%02d' \
        $((seconds / 3600)) \
        $(((seconds % 3600) / 60)) \
        $((seconds % 60))
}

show_stage() {
    local percent="$1"
    local text="$2"
    local now elapsed

    now=$(date +%s)
    elapsed=$((now - GLOBAL_START))

    printf '\n[%3d%%] %-50s Elapsed: %s\n' \
        "$percent" \
        "$text" \
        "$(format_time "$elapsed")"
}

fail() {
    echo
    echo "================================================"
    echo "VACMAP FRESH WHOLE-GENOME TEST: FAIL"
    echo "Reason: $1"
    echo "Log: $LOG"
    echo "================================================"
    exit 1
}

available_memory_mb() {
    awk '
        /^MemAvailable:/ {
            print int($2 / 1024)
        }
    ' /proc/meminfo
}

free_swap_mb() {
    awk '
        /^SwapFree:/ {
            print int($2 / 1024)
        }
    ' /proc/meminfo
}

free_disk_gb() {
    df -BG "$OUT" |
    awk '
        NR == 2 {
            gsub(/G/, "", $4)
            print $4
        }
    '
}

cd "$PROJECT"

mkdir -p "$OUT"

show_stage 0 "Starting fresh whole-genome-index diagnostic"

echo
echo "Purpose:"
echo "  Test the same 20 HG002 HiFi reads using a new GRCh38 VACmap index."
echo
echo "Estimated duration: 20–60 minutes"
echo "Progress type: stage-based estimate"
echo

show_stage 5 "Checking files, programs and resources"

command -v samtools >/dev/null 2>&1 ||
    fail "samtools was not found"

command -v conda >/dev/null 2>&1 ||
    fail "conda was not found"

[[ -s "$ORIGINAL_REFERENCE" ]] ||
    fail "Reference is missing: $ORIGINAL_REFERENCE"

[[ -s "$READS" ]] ||
    fail "20-read FASTQ is missing: $READS"

READ_COUNT=$(
    awk '
        END {
            print int(NR / 4)
        }
    ' "$READS"
)

[[ "$READ_COUNT" -eq 20 ]] ||
    fail "Expected 20 input reads but found $READ_COUNT"

echo "Input reads:          $READ_COUNT"
echo "Available RAM:        $(available_memory_mb) MiB"
echo "Free swap:            $(free_swap_mb) MiB"
echo "Free project disk:    $(free_disk_gb) GiB"

show_stage 10 "Preparing isolated reference location"

# Remove only previous diagnostic outputs.
# The original reference and original VACmap index are not touched.

/bin/rm -f \
    "$FRESH_REFERENCE" \
    "${FRESH_REFERENCE}.fai" \
    "$BAM" \
    "${BAM}.bai" \
    "$LOG" \
    "$SUMMARY"

find "$OUT" \
    -maxdepth 1 \
    -type f \
    -name 'reference_fresh.fasta.w*_k*.mmi' \
    -print \
    -delete

# A hard link gives VACmap a new path without duplicating the
# multi-gigabyte reference file. VACmap will write its diagnostic
# index beside this new path.

if ln "$ORIGINAL_REFERENCE" "$FRESH_REFERENCE" 2>/dev/null
then
    echo "Fresh reference path created using a hard link."
else
    echo "Hard link unavailable; creating an independent reference copy."

    cp --reflink=auto \
        "$ORIGINAL_REFERENCE" \
        "$FRESH_REFERENCE"
fi

[[ -s "$FRESH_REFERENCE" ]] ||
    fail "Fresh reference path was not created"

ORIGINAL_SHA=$(
    sha256sum "$ORIGINAL_REFERENCE" |
    awk '{print $1}'
)

FRESH_SHA=$(
    sha256sum "$FRESH_REFERENCE" |
    awk '{print $1}'
)

[[ "$ORIGINAL_SHA" == "$FRESH_SHA" ]] ||
    fail "Fresh reference does not match the original reference"

echo "Original reference SHA256: $ORIGINAL_SHA"
echo "Fresh reference SHA256:    $FRESH_SHA"

show_stage 15 "Confirming that no diagnostic index exists"

EXISTING_INDEX=$(
    find "$OUT" \
        -maxdepth 1 \
        -type f \
        -name 'reference_fresh.fasta.w*_k*.mmi' \
        -print \
        -quit
)

if [[ -n "$EXISTING_INDEX" ]]; then
    fail "A diagnostic VACmap index already exists: $EXISTING_INDEX"
fi

echo "Fresh-index condition: PASS"

show_stage 20 "Starting VACmap fresh-index run"

touch "$LOG"

RUN_START=$(date +%s)
MAPPING_START=0

set +e

conda run -n vacmap_env vacmap \
    -ref "$FRESH_REFERENCE" \
    -read "$READS" \
    -mode L \
    -k 15 \
    -w 10 \
    -t 4 \
    --rg-id HG002_fresh_index_test \
    --rg-sm HG002 \
    --rg-pl PACBIO \
    --force \
    -o "$BAM" \
    > "$LOG" 2>&1 &

VACMAP_PID=$!

while kill -0 "$VACMAP_PID" 2>/dev/null
do
    NOW=$(date +%s)
    TOTAL_ELAPSED=$((NOW - GLOBAL_START))
    RUN_ELAPSED=$((NOW - RUN_START))

    MEMORY_MB=$(available_memory_mb)
    SWAP_MB=$(free_swap_mb)
    DISK_GB=$(free_disk_gb)

    if grep -qF "Reading $READS" "$LOG" 2>/dev/null
    then
        if (( MAPPING_START == 0 )); then
            MAPPING_START=$NOW
        fi

        MAPPING_ELAPSED=$((NOW - MAPPING_START))

        PERCENT=$(
            awk \
                -v elapsed="$MAPPING_ELAPSED" \
                -v expected="$EXPECTED_MAPPING_SECONDS" \
                'BEGIN {
                    p = 60 + int((elapsed / expected) * 30)

                    if (p > 90) {
                        p = 90
                    }

                    print p
                }'
        )

        REMAINING=$((EXPECTED_MAPPING_SECONDS - MAPPING_ELAPSED))

        if (( REMAINING > 0 )); then
            ETA=$(format_time "$REMAINING")
        else
            ETA="re-estimating"
        fi

        STAGE="Mapping 20 HiFi reads"
    else
        PERCENT=$(
            awk \
                -v elapsed="$RUN_ELAPSED" \
                -v expected="$EXPECTED_INDEX_SECONDS" \
                'BEGIN {
                    p = 20 + int((elapsed / expected) * 35)

                    if (p > 55) {
                        p = 55
                    }

                    print p
                }'
        )

        REMAINING=$(
            (
                EXPECTED_INDEX_SECONDS
                + EXPECTED_MAPPING_SECONDS
                - RUN_ELAPSED
            )
        )

        if (( REMAINING > 0 )); then
            ETA=$(format_time "$REMAINING")
        else
            ETA="re-estimating"
        fi

        STAGE="Building fresh whole-genome index"
    fi

    printf \
        '\r[%3d%%] %s | Elapsed: %s | ETA: %s | RAM available: %s MiB | Swap free: %s MiB | Disk: %s GiB   ' \
        "$PERCENT" \
        "$STAGE" \
        "$(format_time "$TOTAL_ELAPSED")" \
        "$ETA" \
        "$MEMORY_MB" \
        "$SWAP_MB" \
        "$DISK_GB"

    # Emergency protection against complete memory exhaustion.
    if (( MEMORY_MB + SWAP_MB < 384 )); then
        echo
        echo "ERROR: available RAM plus free swap is critically low."

        kill "$VACMAP_PID" 2>/dev/null || true
        wait "$VACMAP_PID" 2>/dev/null || true

        fail "VACmap was stopped to avoid memory exhaustion"
    fi

    sleep 10
done

wait "$VACMAP_PID"
VACMAP_STATUS=$?

set -e

echo

RUN_END=$(date +%s)
RUN_SECONDS=$((RUN_END - RUN_START))

show_stage 92 "VACmap process finished"

echo "VACmap exit code: $VACMAP_STATUS"
echo "VACmap runtime:   $(format_time "$RUN_SECONDS")"

if [[ "$VACMAP_STATUS" -ne 0 ]]; then
    echo
    echo "=== VACMAP LOG ==="
    cat "$LOG"

    fail "VACmap returned a non-zero exit code"
fi

show_stage 94 "Checking the newly created index"

INDEX_FILE=$(
    find "$OUT" \
        -maxdepth 1 \
        -type f \
        -name 'reference_fresh.fasta.w10_k15.mmi' \
        -print \
        -quit
)

if [[ -n "$INDEX_FILE" ]]; then
    echo "Fresh index created:"
    ls -lh "$INDEX_FILE"
else
    echo "WARNING: expected saved index was not found."
    echo "The BAM will still be validated."
fi

show_stage 96 "Validating the VACmap BAM"

[[ -s "$BAM" ]] ||
    fail "VACmap did not create a non-empty BAM file"

samtools quickcheck "$BAM" ||
    fail "VACmap BAM failed samtools quickcheck"

TOTAL_RECORDS=$(samtools view -c "$BAM")
PRIMARY_RECORDS=$(samtools view -c -F 0x900 "$BAM")
MAPPED_PRIMARY=$(samtools view -c -F 0x904 "$BAM")

echo
echo "Input reads:            $READ_COUNT"
echo "Total BAM records:      $TOTAL_RECORDS"
echo "Primary records:        $PRIMARY_RECORDS"
echo "Mapped primary records: $MAPPED_PRIMARY"

show_stage 98 "Writing diagnostic summary"

{
    printf 'input_reads\ttotal_records\tprimary_records\tmapped_primary\truntime_seconds\tbam\tindex\n'

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$READ_COUNT" \
        "$TOTAL_RECORDS" \
        "$PRIMARY_RECORDS" \
        "$MAPPED_PRIMARY" \
        "$RUN_SECONDS" \
        "$BAM" \
        "${INDEX_FILE:-NOT_SAVED}"
} > "$SUMMARY"

column -t -s $'\t' "$SUMMARY"

echo
echo "=== IMPORTANT VACMAP LOG LINES ==="

grep -E \
    'CMD:|Mode:|Threads:|Kmer size:|Loading index:|Reading |sequences processed|User time|WARNING|ERROR' \
    "$LOG" \
    || true

GLOBAL_END=$(date +%s)
TOTAL_SECONDS=$((GLOBAL_END - GLOBAL_START))

echo
echo "================================================"
echo "[100%] VACMAP FRESH WHOLE-GENOME TEST COMPLETE"
echo "================================================"
echo "Total elapsed time:     $(format_time "$TOTAL_SECONDS")"
echo "Input reads:            $READ_COUNT"
echo "Total BAM records:      $TOTAL_RECORDS"
echo "Primary records:        $PRIMARY_RECORDS"
echo "Mapped primary records: $MAPPED_PRIMARY"
echo "Summary:                $SUMMARY"
echo "Log:                    $LOG"
echo

if (( MAPPED_PRIMARY > 0 )); then
    echo "CONCLUSION:"
    echo "The fresh whole-genome index restores VACmap mapping."
    echo "The previous whole-genome VACmap index is stale, incompatible or corrupted."
    echo "STATUS: FRESH INDEX TEST PASS"
else
    echo "CONCLUSION:"
    echo "A fresh whole-genome index still produced zero mapped primary reads."
    echo "The index alone does not explain the problem."
    echo "STATUS: FRESH INDEX TEST ZERO ALIGNMENTS"
fi

echo "================================================"
