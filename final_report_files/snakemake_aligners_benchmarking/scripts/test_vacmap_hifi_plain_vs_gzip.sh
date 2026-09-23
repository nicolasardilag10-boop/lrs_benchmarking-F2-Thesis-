#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"
REFERENCE="$PROJECT/data/reference/genome.fasta"
SOURCE_FASTQ="$PROJECT/data/samples/HG002.pb.1k.fastq.gz"
OUT="$PROJECT/results/vacmap_diagnostic/hg002_hifi_input_test"

PLAIN_FASTQ="$OUT/HG002.pb.first20.fastq"
GZIP_FASTQ="$OUT/HG002.pb.first20.fastq.gz"

PLAIN_BAM="$OUT/HG002.pb.plain.vacmap.sorted.bam"
GZIP_BAM="$OUT/HG002.pb.gzip.vacmap.sorted.bam"

PLAIN_LOG="$OUT/HG002.pb.plain.vacmap.log"
GZIP_LOG="$OUT/HG002.pb.gzip.vacmap.log"

INITIAL_EXPECTED_SECONDS=900
GLOBAL_START=$(date +%s)

format_time() {
    local seconds="${1:-0}"

    (( seconds < 0 )) && seconds=0

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

    printf '\n[%3d%%] %-48s Elapsed: %s\n' \
        "$percent" \
        "$text" \
        "$(format_time "$elapsed")"
}

fail() {
    echo
    echo "================================================"
    echo "VACMAP INPUT TEST: FAIL"
    echo "Reason: $1"
    echo "================================================"
    exit 1
}

monitor_process() {
    local pid="$1"
    local start_percent="$2"
    local end_percent="$3"
    local expected_seconds="$4"
    local label="$5"
    local run_start="$6"

    while kill -0 "$pid" 2>/dev/null
    do
        local now run_elapsed global_elapsed span percent remaining

        now=$(date +%s)
        run_elapsed=$((now - run_start))
        global_elapsed=$((now - GLOBAL_START))
        span=$((end_percent - start_percent))

        percent=$(
            awk \
                -v start="$start_percent" \
                -v span="$span" \
                -v elapsed="$run_elapsed" \
                -v expected="$expected_seconds" \
                'BEGIN {
                    value = start + int(span * elapsed / expected)

                    if (value >= start + span) {
                        value = start + span - 1
                    }

                    if (value < start) {
                        value = start
                    }

                    print value
                }'
        )

        remaining=$((expected_seconds - run_elapsed))

        if (( remaining > 0 )); then
            eta="$(format_time "$remaining")"
        else
            eta="estimating..."
        fi

        printf '\r[%3d%%] %s | Elapsed: %s | ETA: %s   ' \
            "$percent" \
            "$label" \
            "$(format_time "$global_elapsed")" \
            "$eta"

        sleep 10
    done

    echo
}

count_bam() {
    local bam="$1"
    local prefix="$2"

    if [[ ! -s "$bam" ]]; then
        printf -v "${prefix}_TOTAL" '%s' "0"
        printf -v "${prefix}_PRIMARY" '%s' "0"
        printf -v "${prefix}_MAPPED" '%s' "0"
        return
    fi

    samtools quickcheck "$bam" ||
        fail "BAM integrity failure: $bam"

    printf -v "${prefix}_TOTAL" \
        '%s' \
        "$(samtools view -c "$bam")"

    printf -v "${prefix}_PRIMARY" \
        '%s' \
        "$(samtools view -c -F 0x900 "$bam")"

    printf -v "${prefix}_MAPPED" \
        '%s' \
        "$(samtools view -c -F 0x904 "$bam")"
}

cd "$PROJECT"
mkdir -p "$OUT"

show_stage 0 "Starting controlled VACmap HiFi test"

echo
echo "Purpose:"
echo "  Compare identical reads as plain FASTQ and FASTQ.GZ"
echo
echo "Initial estimated duration: 20–45 minutes"
echo "Progress type: stage-based estimate"
echo

show_stage 5 "Checking programs and input files"

command -v samtools >/dev/null ||
    fail "samtools was not found"

command -v gzip >/dev/null ||
    fail "gzip was not found"

[[ -s "$REFERENCE" ]] ||
    fail "Reference is missing: $REFERENCE"

[[ -s "$SOURCE_FASTQ" ]] ||
    fail "Input FASTQ is missing: $SOURCE_FASTQ"

show_stage 10 "Creating identical 20-read test inputs"

/bin/rm -f \
    "$PLAIN_FASTQ" \
    "$GZIP_FASTQ" \
    "$PLAIN_BAM" \
    "${PLAIN_BAM}.bai" \
    "$GZIP_BAM" \
    "${GZIP_BAM}.bai" \
    "$PLAIN_LOG" \
    "$GZIP_LOG"

gzip -cd "$SOURCE_FASTQ" |
awk 'NR <= 80 {print}' \
    > "$PLAIN_FASTQ"

gzip -c "$PLAIN_FASTQ" \
    > "$GZIP_FASTQ"

PLAIN_READS=$(awk 'END {print int(NR / 4)}' "$PLAIN_FASTQ")
GZIP_READS=$(gzip -cd "$GZIP_FASTQ" | awk 'END {print int(NR / 4)}')

echo "Plain FASTQ reads:      $PLAIN_READS"
echo "Compressed FASTQ reads: $GZIP_READS"

[[ "$PLAIN_READS" -eq 20 ]] ||
    fail "Plain FASTQ does not contain exactly 20 reads"

[[ "$GZIP_READS" -eq 20 ]] ||
    fail "Compressed FASTQ does not contain exactly 20 reads"

PLAIN_SHA=$(
    sha256sum "$PLAIN_FASTQ" |
    awk '{print $1}'
)

GZIP_CONTENT_SHA=$(
    gzip -cd "$GZIP_FASTQ" |
    sha256sum |
    awk '{print $1}'
)

echo "Plain content SHA256:   $PLAIN_SHA"
echo "Gzip content SHA256:    $GZIP_CONTENT_SHA"

[[ "$PLAIN_SHA" == "$GZIP_CONTENT_SHA" ]] ||
    fail "The two FASTQ inputs are not identical"

show_stage 15 "Starting plain-FASTQ VACmap run"

PLAIN_START=$(date +%s)

set +e

conda run -n vacmap_env vacmap \
    -ref "$REFERENCE" \
    -read "$PLAIN_FASTQ" \
    -mode L \
    -t 4 \
    --rg-id HG002 \
    --rg-pl PACBIO \
    --force \
    -o "$PLAIN_BAM" \
    > "$PLAIN_LOG" 2>&1 &

PLAIN_PID=$!

monitor_process \
    "$PLAIN_PID" \
    15 \
    50 \
    "$INITIAL_EXPECTED_SECONDS" \
    "VACmap plain FASTQ running" \
    "$PLAIN_START"

wait "$PLAIN_PID"
PLAIN_STATUS=$?

set -e

PLAIN_END=$(date +%s)
PLAIN_SECONDS=$((PLAIN_END - PLAIN_START))

show_stage 50 "Plain-FASTQ VACmap run finished"

echo "Exit code: $PLAIN_STATUS"
echo "Duration:  $(format_time "$PLAIN_SECONDS")"

if [[ "$PLAIN_STATUS" -ne 0 ]]; then
    echo
    echo "=== PLAIN FASTQ VACMAP LOG ==="
    cat "$PLAIN_LOG"

    fail "Plain FASTQ VACmap run failed"
fi

show_stage 55 "Validating plain-FASTQ BAM"

count_bam "$PLAIN_BAM" PLAIN

echo "Total records:          $PLAIN_TOTAL"
echo "Primary records:        $PLAIN_PRIMARY"
echo "Mapped primary records: $PLAIN_MAPPED"

SECOND_EXPECTED=$PLAIN_SECONDS

if (( SECOND_EXPECTED < 60 )); then
    SECOND_EXPECTED=60
fi

show_stage 60 "Starting compressed-FASTQ VACmap run"

echo
echo "Second-run ETA based on first run: $(format_time "$SECOND_EXPECTED")"
echo

GZIP_START=$(date +%s)

set +e

conda run -n vacmap_env vacmap \
    -ref "$REFERENCE" \
    -read "$GZIP_FASTQ" \
    -mode L \
    -t 4 \
    --rg-id HG002 \
    --rg-pl PACBIO \
    --force \
    -o "$GZIP_BAM" \
    > "$GZIP_LOG" 2>&1 &

GZIP_PID=$!

monitor_process \
    "$GZIP_PID" \
    60 \
    90 \
    "$SECOND_EXPECTED" \
    "VACmap compressed FASTQ running" \
    "$GZIP_START"

wait "$GZIP_PID"
GZIP_STATUS=$?

set -e

GZIP_END=$(date +%s)
GZIP_SECONDS=$((GZIP_END - GZIP_START))

show_stage 90 "Compressed-FASTQ VACmap run finished"

echo "Exit code: $GZIP_STATUS"
echo "Duration:  $(format_time "$GZIP_SECONDS")"

if [[ "$GZIP_STATUS" -ne 0 ]]; then
    echo
    echo "=== COMPRESSED FASTQ VACMAP LOG ==="
    cat "$GZIP_LOG"

    fail "Compressed FASTQ VACmap run failed"
fi

show_stage 94 "Validating compressed-FASTQ BAM"

count_bam "$GZIP_BAM" GZIP

echo "Total records:          $GZIP_TOTAL"
echo "Primary records:        $GZIP_PRIMARY"
echo "Mapped primary records: $GZIP_MAPPED"

show_stage 97 "Comparing plain and compressed results"

echo
echo "================================================"
printf '%-22s %12s %12s\n' \
    "Metric" \
    "Plain FASTQ" \
    "FASTQ.GZ"

printf '%-22s %12s %12s\n' \
    "Input reads" \
    "$PLAIN_READS" \
    "$GZIP_READS"

printf '%-22s %12s %12s\n' \
    "Total BAM records" \
    "$PLAIN_TOTAL" \
    "$GZIP_TOTAL"

printf '%-22s %12s %12s\n' \
    "Primary records" \
    "$PLAIN_PRIMARY" \
    "$GZIP_PRIMARY"

printf '%-22s %12s %12s\n' \
    "Mapped primary" \
    "$PLAIN_MAPPED" \
    "$GZIP_MAPPED"

printf '%-22s %12s %12s\n' \
    "Runtime" \
    "$(format_time "$PLAIN_SECONDS")" \
    "$(format_time "$GZIP_SECONDS")"

echo "================================================"

if (( PLAIN_MAPPED > 0 && GZIP_MAPPED == 0 )); then
    CONCLUSION="Compressed FASTQ handling is the likely cause."
elif (( PLAIN_MAPPED == 0 && GZIP_MAPPED == 0 )); then
    CONCLUSION="Compression is not the main cause; inspect mode, reference index, or VACmap parameters."
elif (( PLAIN_MAPPED > 0 && GZIP_MAPPED > 0 )); then
    CONCLUSION="Both input formats work; the original workflow command or output handling is the likely cause."
else
    CONCLUSION="The results are inconsistent and require log inspection."
fi

TOTAL_END=$(date +%s)
TOTAL_SECONDS=$((TOTAL_END - GLOBAL_START))

echo
echo "================================================"
echo "[100%] VACMAP HIFI INPUT TEST COMPLETE"
echo "================================================"
echo "Total elapsed time: $(format_time "$TOTAL_SECONDS")"
echo
echo "Conclusion:"
echo "$CONCLUSION"
echo
echo "Plain log:      $PLAIN_LOG"
echo "Compressed log: $GZIP_LOG"
echo "STATUS: PASS"
echo "================================================"
