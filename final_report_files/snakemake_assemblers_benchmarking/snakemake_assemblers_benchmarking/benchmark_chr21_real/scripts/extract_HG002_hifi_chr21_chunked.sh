#!/usr/bin/env bash

set -euo pipefail


ROOT="benchmark_chr21_real"

BAM_URL="https://downloads.pacbcloud.com/public/dataset/HG002-CpG-methylation-202202/HG002.GRCh38.haplotagged.bam"

OUTPUT_DIR="${ROOT}/inputs/pacbio/HG002"
CHUNK_DIR="${OUTPUT_DIR}/chr21_chunks"

LOG_DIR="${ROOT}/results/logs/extraction/HG002_hifi_chr21_chunks"

FINAL_BAM="${OUTPUT_DIR}/HG002.hifi.chr21.primary.q20.bam"

SUMMARY="${ROOT}/results/validation/HG002.hifi.chr21.extraction.tsv"

FILTER_SCRIPT="${ROOT}/scripts/filter_sam_window.py"

CHROMOSOME="chr21"
CHROMOSOME_LENGTH=46709983
WINDOW_SIZE=1000000

# Exclude unmapped, secondary and supplementary alignments.
EXCLUDED_FLAGS=2308
MINIMUM_MAPQ=20

MAX_ATTEMPTS=5


format_seconds() {
    local total_seconds="$1"

    local hours=$((total_seconds / 3600))
    local minutes=$(((total_seconds % 3600) / 60))
    local seconds=$((total_seconds % 60))

    printf '%02d:%02d:%02d' \
        "$hours" \
        "$minutes" \
        "$seconds"
}


fail() {
    echo
    echo "ERROR: $1" >&2
    exit 1
}


for program in samtools python3 awk
do
    command -v "$program" >/dev/null 2>&1 \
        || fail "Required program not found: $program"
done


[[ -f "$FILTER_SCRIPT" ]] \
    || fail "Missing filter script: $FILTER_SCRIPT"


if pgrep -f \
    '[f]lye-modules|[s]nakemake.*benchmark_1k/Snakefile' \
    >/dev/null 2>&1
then
    fail "An old benchmark process is still running."
fi


mkdir -p \
    "$OUTPUT_DIR" \
    "$CHUNK_DIR" \
    "$LOG_DIR" \
    "$(dirname "$SUMMARY")"


TOTAL_WINDOWS=$(
    (
        CHROMOSOME_LENGTH
        + WINDOW_SIZE
        - 1
    ) / WINDOW_SIZE
)

RUN_START=$SECONDS
COMPLETED_WINDOWS=0
TOTAL_SELECTED=0


echo "[  0%] Starting resumable HG002 HiFi chr21 extraction."
echo "       Chromosome length: ${CHROMOSOME_LENGTH} bp"
echo "       Number of windows: ${TOTAL_WINDOWS}"
echo "       Window size: ${WINDOW_SIZE} bp"
echo "       Maximum attempts per window: ${MAX_ATTEMPTS}"
echo "       ETA will stabilize after approximately 3–5 windows."
echo


rm -f \
    "${FINAL_BAM}.partial" \
    "${OUTPUT_DIR}/HG002.hifi.chr21.primary.q20.bam.partial"


for ((window_number = 1; window_number <= TOTAL_WINDOWS; window_number++))
do
    window_start=$(
        (window_number - 1) * WINDOW_SIZE + 1
    )

    window_end=$(
        window_number * WINDOW_SIZE
    )

    if ((window_end > CHROMOSOME_LENGTH))
    then
        window_end="$CHROMOSOME_LENGTH"
    fi

    region="${CHROMOSOME}:${window_start}-${window_end}"

    chunk="$(
        printf \
            '%s/chunk_%03d_%09d_%09d.bam' \
            "$CHUNK_DIR" \
            "$window_number" \
            "$window_start" \
            "$window_end"
    )"

    partial="${chunk}.partial"


    if [[ -s "$chunk" ]] \
        && samtools quickcheck "$chunk" 2>/dev/null
    then
        selected="$(samtools view -c "$chunk")"

        TOTAL_SELECTED=$((TOTAL_SELECTED + selected))
        COMPLETED_WINDOWS=$((COMPLETED_WINDOWS + 1))

        elapsed=$((SECONDS - RUN_START))

        if ((COMPLETED_WINDOWS > 0))
        then
            estimated_total=$(
                elapsed * TOTAL_WINDOWS
                / COMPLETED_WINDOWS
            )

            eta=$((estimated_total - elapsed))

            if ((eta < 0))
            then
                eta=0
            fi
        else
            eta=0
        fi

        overall_percent=$(
            10
            + COMPLETED_WINDOWS * 70
            / TOTAL_WINDOWS
        )

        printf \
            '[%3d%%] Reused window %d/%d | selected %d | elapsed %s | ETA %s\n' \
            "$overall_percent" \
            "$window_number" \
            "$TOTAL_WINDOWS" \
            "$selected" \
            "$(format_seconds "$elapsed")" \
            "$(format_seconds "$eta")"

        continue
    fi


    success=0

    for ((attempt = 1; attempt <= MAX_ATTEMPTS; attempt++))
    do
        rm -f "$partial"

        log="$(
            printf \
                '%s/window_%03d_attempt_%d.log' \
                "$LOG_DIR" \
                "$window_number" \
                "$attempt"
        )"

        elapsed=$((SECONDS - RUN_START))

        if ((COMPLETED_WINDOWS > 0))
        then
            estimated_total=$(
                elapsed * TOTAL_WINDOWS
                / COMPLETED_WINDOWS
            )

            eta=$((estimated_total - elapsed))

            if ((eta < 0))
            then
                eta=0
            fi
        else
            eta=0
        fi

        overall_percent=$(
            10
            + COMPLETED_WINDOWS * 70
            / TOTAL_WINDOWS
        )

        printf \
            '[%3d%%] Window %d/%d: %s | attempt %d/%d | elapsed %s | ETA %s\n' \
            "$overall_percent" \
            "$window_number" \
            "$TOTAL_WINDOWS" \
            "$region" \
            "$attempt" \
            "$MAX_ATTEMPTS" \
            "$(format_seconds "$elapsed")" \
            "$(format_seconds "$eta")"


        set +e

        samtools view \
            -h \
            "$BAM_URL" \
            "$region" \
            2> "$log" \
        | python3 "$FILTER_SCRIPT" \
            "$CHROMOSOME" \
            "$window_start" \
            "$window_end" \
            "$EXCLUDED_FLAGS" \
            "$MINIMUM_MAPQ" \
            2>> "$log" \
        | samtools view \
            -@ 2 \
            -b \
            -o "$partial" \
            - \
            2>> "$log"

        statuses=("${PIPESTATUS[@]}")

        set -e


        if [[
            "${statuses[0]}" -eq 0
            &&
            "${statuses[1]}" -eq 0
            &&
            "${statuses[2]}" -eq 0
            &&
            -s "$partial"
        ]] \
            && samtools quickcheck "$partial" 2>/dev/null
        then
            mv "$partial" "$chunk"
            success=1
            break
        fi


        rm -f "$partial"

        delay=$((attempt * 10))

        echo "       Window failed."
        echo "       Retrying in ${delay} seconds."
        echo "       Log: ${log}"

        sleep "$delay"
    done


    if [[ "$success" -ne 1 ]]
    then
        fail \
            "Window ${window_number} failed after ${MAX_ATTEMPTS} attempts. Rerun the same command to resume."
    fi


    selected="$(samtools view -c "$chunk")"

    TOTAL_SELECTED=$((TOTAL_SELECTED + selected))
    COMPLETED_WINDOWS=$((COMPLETED_WINDOWS + 1))

    elapsed=$((SECONDS - RUN_START))

    estimated_total=$(
        elapsed * TOTAL_WINDOWS
        / COMPLETED_WINDOWS
    )

    eta=$((estimated_total - elapsed))

    if ((eta < 0))
    then
        eta=0
    fi

    overall_percent=$(
        10
        + COMPLETED_WINDOWS * 70
        / TOTAL_WINDOWS
    )

    printf \
        '[%3d%%] Completed window %d/%d | selected %d | total %d | elapsed %s | ETA %s\n' \
        "$overall_percent" \
        "$window_number" \
        "$TOTAL_WINDOWS" \
        "$selected" \
        "$TOTAL_SELECTED" \
        "$(format_seconds "$elapsed")" \
        "$(format_seconds "$eta")"
done


echo
echo "[ 82%] All chromosome windows completed."
echo "[ 85%] Concatenating chromosome-21 BAM chunks."


mapfile -t CHUNKS < <(
    find "$CHUNK_DIR" \
        -maxdepth 1 \
        -type f \
        -name 'chunk_*.bam' \
        | sort
)


if [[ "${#CHUNKS[@]}" -ne "$TOTAL_WINDOWS" ]]
then
    fail \
        "Expected ${TOTAL_WINDOWS} BAM chunks but found ${#CHUNKS[@]}."
fi


rm -f \
    "${FINAL_BAM}.partial" \
    "$FINAL_BAM" \
    "${FINAL_BAM}.bai"


samtools cat \
    -o "${FINAL_BAM}.partial" \
    "${CHUNKS[@]}"


samtools quickcheck "${FINAL_BAM}.partial" \
    || fail "Concatenated BAM failed integrity validation."


mv \
    "${FINAL_BAM}.partial" \
    "$FINAL_BAM"


echo "[ 90%] Indexing final chromosome-21 BAM."

samtools index \
    -@ 2 \
    "$FINAL_BAM"


echo "[ 94%] Calculating read count and estimated input coverage."


READS="$(samtools view -c "$FINAL_BAM")"


TOTAL_BASES="$(
    samtools view "$FINAL_BAM" \
    | awk '
        $10 != "*" {
            total += length($10)
        }

        END {
            print total + 0
        }
    '
)"


COVERAGE="$(
    python3 - \
        "$TOTAL_BASES" \
        "$CHROMOSOME_LENGTH" \
    <<'PY'
import sys

total_bases = int(sys.argv[1])
chromosome_length = int(sys.argv[2])

print(f"{total_bases / chromosome_length:.4f}")
PY
)"


COVERAGE_STATUS="$(
    python3 - "$COVERAGE" <<'PY'
import sys

coverage = float(sys.argv[1])

if coverage >= 30:
    print("READY_FOR_30X_DOWNSAMPLING")
else:
    print("BELOW_30X")
PY
)"


cat > "$SUMMARY" <<EOF
sample	technology	chromosome	filter	reads	total_read_bases	chromosome_length	estimated_input_coverage	target_coverage	coverage_status	bam
HG002	PacBio_HiFi	chr21	primary_MAPQ20	${READS}	${TOTAL_BASES}	${CHROMOSOME_LENGTH}	${COVERAGE}	30	${COVERAGE_STATUS}	${FINAL_BAM}
EOF


echo "[ 98%] Running final BAM integrity check."

samtools quickcheck "$FINAL_BAM" \
    || fail "Final BAM integrity check failed."


TOTAL_ELAPSED=$((SECONDS - RUN_START))


echo "[100%] HG002 HiFi chromosome-21 extraction completed."
echo
echo "Total elapsed time: $(format_seconds "$TOTAL_ELAPSED")"
echo

column -t -s $'\t' "$SUMMARY"

echo
echo "=== OUTPUT FILES ==="

ls -lh \
    "$FINAL_BAM" \
    "${FINAL_BAM}.bai" \
    "$SUMMARY"
