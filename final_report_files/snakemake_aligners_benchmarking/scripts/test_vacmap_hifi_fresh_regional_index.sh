#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"

REFERENCE="$PROJECT/data/reference/genome.fasta"
SOURCE_FASTQ="$PROJECT/data/samples/HG002.pb.1k.fastq.gz"
MINIMAP_BAM="$PROJECT/results/bam/HG002.pb.minimap2.sorted.bam"

OUT="$PROJECT/results/vacmap_diagnostic/hg002_hifi_regional_test"

READ_FASTQ="$OUT/selected_read.fastq"
REGION_FASTA="$OUT/reference_region.fasta"
RESULTS="$OUT/results.tsv"

GLOBAL_START=$(date +%s)

# Initial estimate for each small regional run.
EXPECTED_CASE_SECONDS=180

declare -A TOTAL
declare -A PRIMARY
declare -A MAPPED
declare -A RUNTIME

format_time() {
    local seconds="${1:-0}"

    (( seconds < 0 )) && seconds=0

    printf '%02d:%02d:%02d' \
        $((seconds / 3600)) \
        $(((seconds % 3600) / 60)) \
        $((seconds % 60))
}

show_stage() {
    local percentage="$1"
    local message="$2"
    local now elapsed

    now=$(date +%s)
    elapsed=$((now - GLOBAL_START))

    printf '\n[%3d%%] %-48s Elapsed: %s\n' \
        "$percentage" \
        "$message" \
        "$(format_time "$elapsed")"
}

fail() {
    echo
    echo "================================================"
    echo "VACMAP REGIONAL DIAGNOSTIC: FAIL"
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
    local case_start="$6"

    while kill -0 "$pid" 2>/dev/null
    do
        local now case_elapsed global_elapsed span percentage
        local remaining eta

        now=$(date +%s)

        case_elapsed=$((now - case_start))
        global_elapsed=$((now - GLOBAL_START))
        span=$((end_percent - start_percent))

        percentage=$(
            awk \
                -v start="$start_percent" \
                -v span="$span" \
                -v elapsed="$case_elapsed" \
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

        remaining=$((expected_seconds - case_elapsed))

        if (( remaining > 0 )); then
            eta=$(format_time "$remaining")
        else
            eta="estimating..."
        fi

        printf '\r[%3d%%] %s | Elapsed: %s | ETA: %s   ' \
            "$percentage" \
            "$label" \
            "$(format_time "$global_elapsed")" \
            "$eta"

        sleep 5
    done

    echo
}

run_case() {
    local name="$1"
    local mode="$2"
    local use_quality_ignore="$3"
    local start_percent="$4"
    local end_percent="$5"

    local bam="$OUT/${name}.sorted.bam"
    local log="$OUT/${name}.log"
    local case_reference="$OUT/reference_region.${name}.fasta"

    local case_start case_end status
    local extra_arguments=()

    if [[ "$use_quality_ignore" == "yes" ]]; then
        extra_arguments+=(--Q)
    fi

    cp "$REGION_FASTA" "$case_reference"

    /bin/rm -f \
        "$bam" \
        "${bam}.bai" \
        "$log"

    find "$OUT" \
        -maxdepth 1 \
        -type f \
        -name "$(basename "$case_reference").w*_k*.mmi" \
        -delete

    show_stage "$start_percent" "Starting $name"

    echo "Mode: $mode"
    echo "Ignore base quality: $use_quality_ignore"
    echo "Reference: $case_reference"
    echo

    case_start=$(date +%s)

    set +e

    conda run -n vacmap_env vacmap \
        -ref "$case_reference" \
        -read "$READ_FASTQ" \
        -mode "$mode" \
        -t 4 \
        --rg-id HG002 \
        --rg-pl PACBIO \
        --nowriteindex \
        --force \
        "${extra_arguments[@]}" \
        -o "$bam" \
        > "$log" 2>&1 &

    local pid=$!

    monitor_process \
        "$pid" \
        "$start_percent" \
        "$end_percent" \
        "$EXPECTED_CASE_SECONDS" \
        "$name running" \
        "$case_start"

    wait "$pid"
    status=$?

    set -e

    case_end=$(date +%s)
    RUNTIME["$name"]=$((case_end - case_start))

    if [[ "$status" -ne 0 ]]; then
        echo
        echo "=== $name LOG ==="
        cat "$log"

        fail "$name returned exit code $status"
    fi

    if [[ ! -s "$bam" ]]; then
        fail "$name did not produce a BAM"
    fi

    samtools quickcheck "$bam" ||
        fail "$name BAM failed samtools quickcheck"

    TOTAL["$name"]=$(samtools view -c "$bam")
    PRIMARY["$name"]=$(samtools view -c -F 0x900 "$bam")
    MAPPED["$name"]=$(samtools view -c -F 0x904 "$bam")

    show_stage "$end_percent" "$name finished"

    echo "Runtime:                $(format_time "${RUNTIME[$name]}")"
    echo "Total records:          ${TOTAL[$name]}"
    echo "Primary records:        ${PRIMARY[$name]}"
    echo "Mapped primary records: ${MAPPED[$name]}"

    echo
    echo "Important log lines:"

    grep -E \
        'CMD:|Mode:|Kmer size:|Loading index:|Building|Reading |sequences processed|User time|WARNING|ERROR' \
        "$log" \
        || true
}

cd "$PROJECT"

source "$HOME/miniconda3/etc/profile.d/conda.sh"

mkdir -p "$OUT"

/bin/rm -f \
    "$READ_FASTQ" \
    "$REGION_FASTA" \
    "$RESULTS"

show_stage 0 "Starting VACmap regional diagnostic"

echo
echo "Estimated total duration: 3–12 minutes"
echo "Progress method: three controlled mapping stages"
echo

show_stage 5 "Checking inputs and programs"

command -v samtools >/dev/null ||
    fail "samtools was not found"

command -v python3 >/dev/null ||
    fail "python3 was not found"

command -v conda >/dev/null ||
    fail "conda was not found"

[[ -s "$REFERENCE" ]] ||
    fail "Reference missing: $REFERENCE"

[[ -s "$SOURCE_FASTQ" ]] ||
    fail "FASTQ missing: $SOURCE_FASTQ"

[[ -s "$MINIMAP_BAM" ]] ||
    fail "Minimap2 BAM missing: $MINIMAP_BAM"

samtools quickcheck "$MINIMAP_BAM" ||
    fail "Minimap2 BAM failed quickcheck"

if [[ ! -s "${REFERENCE}.fai" ]]; then
    samtools faidx "$REFERENCE"
fi

show_stage 10 "Selecting a confidently mapped HiFi read"

SELECTED=$(
    samtools view -F 0x904 "$MINIMAP_BAM" |
    awk -F '\t' '
        $5 >= 50 && !found {
            selected = $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6
            found = 1
        }

        END {
            if (found) {
                print selected
            }
        }
    '
)

if [[ -z "$SELECTED" ]]; then
    SELECTED=$(
        samtools view -F 0x904 "$MINIMAP_BAM" |
        awk -F '\t' '
            !found {
                selected = $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6
                found = 1
            }

            END {
                if (found) {
                    print selected
                }
            }
        '
    )
fi

[[ -n "$SELECTED" ]] ||
    fail "No mapped primary read was found in the minimap2 BAM"

IFS=$'\t' read -r QNAME CONTIG POSITION MAPQ CIGAR <<< "$SELECTED"

echo "Read name:       $QNAME"
echo "Reference:       $CONTIG"
echo "Position:        $POSITION"
echo "Minimap2 MAPQ:   $MAPQ"
echo "Minimap2 CIGAR:  $CIGAR"

show_stage 15 "Extracting the selected FASTQ record"

python3 - \
    "$SOURCE_FASTQ" \
    "$QNAME" \
    "$READ_FASTQ" <<'PY'
import gzip
import sys
from pathlib import Path

source = Path(sys.argv[1])
target_name = sys.argv[2]
output = Path(sys.argv[3])

found = False

with gzip.open(source, "rt", encoding="utf-8") as source_handle:
    with output.open("w", encoding="utf-8") as output_handle:
        while True:
            header = source_handle.readline()

            if not header:
                break

            sequence = source_handle.readline()
            separator = source_handle.readline()
            quality = source_handle.readline()

            if not quality:
                raise SystemExit("ERROR: incomplete FASTQ record")

            observed_name = header[1:].split()[0]

            if observed_name == target_name:
                output_handle.write(header)
                output_handle.write(sequence)
                output_handle.write(separator)
                output_handle.write(quality)
                found = True
                break

if not found:
    raise SystemExit(
        f"ERROR: read {target_name!r} was not found in the FASTQ"
    )
PY

READ_LENGTH=$(awk 'NR == 2 {print length($0)}' "$READ_FASTQ")
QUALITY_LENGTH=$(awk 'NR == 4 {print length($0)}' "$READ_FASTQ")

[[ "$READ_LENGTH" -gt 0 ]] ||
    fail "Extracted read has zero length"

[[ "$READ_LENGTH" -eq "$QUALITY_LENGTH" ]] ||
    fail "Sequence and quality lengths differ"

echo "Read length:      $READ_LENGTH"
echo "Quality length:   $QUALITY_LENGTH"

show_stage 20 "Creating a fresh regional reference"

CONTIG_LENGTH=$(
    awk \
        -v contig="$CONTIG" \
        '$1 == contig {print $2}' \
        "${REFERENCE}.fai"
)

[[ -n "$CONTIG_LENGTH" ]] ||
    fail "Contig $CONTIG was not found in the reference index"

REGION_START=$((POSITION - 100000))

if (( REGION_START < 1 )); then
    REGION_START=1
fi

REGION_END=$((POSITION + READ_LENGTH + 100000))

if (( REGION_END > CONTIG_LENGTH )); then
    REGION_END=$CONTIG_LENGTH
fi

REGION="${CONTIG}:${REGION_START}-${REGION_END}"

samtools faidx "$REFERENCE" "$REGION" \
    > "$REGION_FASTA"

[[ -s "$REGION_FASTA" ]] ||
    fail "Regional reference extraction failed"

REGION_LENGTH=$(
    awk '
        !/^>/ {
            length_sum += length($0)
        }

        END {
            print length_sum + 0
        }
    ' "$REGION_FASTA"
)

echo "Region:           $REGION"
echo "Region length:    $REGION_LENGTH"
echo "Index policy:     fresh temporary index for every test"

run_case \
    "mode_L" \
    "L" \
    "no" \
    25 \
    45

run_case \
    "mode_L_ignore_quality" \
    "L" \
    "yes" \
    50 \
    70

run_case \
    "mode_H" \
    "H" \
    "no" \
    75 \
    90

show_stage 95 "Comparing results"

{
    printf 'test\tmode\tignore_quality\ttotal_records\tprimary_records\tmapped_primary\truntime_seconds\n'

    printf 'mode_L\tL\tno\t%s\t%s\t%s\t%s\n' \
        "${TOTAL[mode_L]}" \
        "${PRIMARY[mode_L]}" \
        "${MAPPED[mode_L]}" \
        "${RUNTIME[mode_L]}"

    printf 'mode_L_ignore_quality\tL\tyes\t%s\t%s\t%s\t%s\n' \
        "${TOTAL[mode_L_ignore_quality]}" \
        "${PRIMARY[mode_L_ignore_quality]}" \
        "${MAPPED[mode_L_ignore_quality]}" \
        "${RUNTIME[mode_L_ignore_quality]}"

    printf 'mode_H\tH\tno\t%s\t%s\t%s\t%s\n' \
        "${TOTAL[mode_H]}" \
        "${PRIMARY[mode_H]}" \
        "${MAPPED[mode_H]}" \
        "${RUNTIME[mode_H]}"
} > "$RESULTS"

echo
column -t -s $'\t' "$RESULTS"

L_MAPPED=${MAPPED[mode_L]}
LQ_MAPPED=${MAPPED[mode_L_ignore_quality]}
H_MAPPED=${MAPPED[mode_H]}

if (( L_MAPPED > 0 )); then
    CONCLUSION="Mode L works with a fresh regional index. The existing whole-genome VACmap index is the main suspect."
elif (( L_MAPPED == 0 && LQ_MAPPED > 0 )); then
    CONCLUSION="Ignoring FASTQ quality restores mapping. Base-quality handling is the likely cause."
elif (( L_MAPPED == 0 && LQ_MAPPED == 0 && H_MAPPED > 0 )); then
    CONCLUSION="Mode H maps the read but Mode L does not. Mode-L sensitivity or a Mode-L implementation issue is likely."
else
    CONCLUSION="All three tests produced zero mapped primary records. A deeper VACmap/read-reference compatibility investigation is required."
fi

GLOBAL_END=$(date +%s)
TOTAL_SECONDS=$((GLOBAL_END - GLOBAL_START))

echo
echo "================================================"
echo "[100%] VACMAP REGIONAL DIAGNOSTIC COMPLETE"
echo "================================================"
echo "Total elapsed time: $(format_time "$TOTAL_SECONDS")"
echo
echo "Conclusion:"
echo "$CONCLUSION"
echo
echo "Results table:"
echo "$RESULTS"
echo
echo "STATUS: DIAGNOSTIC COMPLETE"
echo "================================================"
