#!/usr/bin/env bash

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
ROOT="$PROJECT/benchmark_chr21_real"
RESULTS="$ROOT/results/assembler_benchmark"
ASSEMBLIES="$RESULTS/assemblies"
CHECKPOINTS="$RESULTS/checkpoints"
LOGROOT="$RESULTS/logs"
CONSOLE="$LOGROOT/real_run.console.log"

TOTAL_ASSEMBLIES=21

count_outputs()
{
    local directory="$1"

    find "$directory" \
        -type f \
        -name 'assembly.fasta' \
        -size +0c \
        2>/dev/null \
    | wc -l
}

format_seconds()
{
    local total="${1:-0}"

    printf "%02d:%02d:%02d" \
        $((total / 3600)) \
        $(((total % 3600) / 60)) \
        $((total % 60))
}

draw_bar()
{
    local completed="$1"
    local total="$2"
    local width=40

    local filled=$((completed * width / total))
    local empty=$((width - filled))

    local left=""
    local right=""

    printf -v left "%*s" "$filled" ""
    printf -v right "%*s" "$empty" ""

    left="${left// /#}"
    right="${right// /-}"

    local percentage
    percentage="$(
        awk \
            -v completed="$completed" \
            -v total="$total" \
            'BEGIN {
                printf "%.2f", 100 * completed / total
            }'
    )"

    printf "[%s%s] %6.2f%%" \
        "$left" \
        "$right" \
        "$percentage"
}

FLYE="$(count_outputs "$ASSEMBLIES/flye")"
GOLDRUSH="$(count_outputs "$ASSEMBLIES/goldrush")"
VERKKO="$(count_outputs "$ASSEMBLIES/verkko")"
NTLINK="$(count_outputs "$ASSEMBLIES/ntlink_goldrush")"

COMPLETED=$(
    (
        FLYE
        + GOLDRUSH
        + VERKKO
        + NTLINK
    )
)

NOW="$(date +%s)"
ELAPSED=0
ETA_TEXT="estimating after first completed assembly"

if [[ -s "$CHECKPOINTS/run_session_started.epoch" ]]
then
    STARTED="$(cat "$CHECKPOINTS/run_session_started.epoch")"

    if [[ "$STARTED" =~ ^[0-9]+$ ]]
    then
        ELAPSED=$((NOW - STARTED))
    fi
fi

if [[ "$COMPLETED" -gt 0 ]]
then
    AVERAGE=$((ELAPSED / COMPLETED))
    REMAINING=$((TOTAL_ASSEMBLIES - COMPLETED))
    ETA=$((AVERAGE * REMAINING))
    ETA_TEXT="$(format_seconds "$ETA")"
fi

echo "============================================================================"
echo "CHROMOSOME 21 ASSEMBLER BENCHMARK"
echo "============================================================================"

draw_bar "$COMPLETED" "$TOTAL_ASSEMBLIES"

echo
echo
printf "Completed assemblies:  %d / %d\n" \
    "$COMPLETED" \
    "$TOTAL_ASSEMBLIES"

printf "Flye:                 %d / 6\n" "$FLYE"
printf "GoldRush:             %d / 6\n" "$GOLDRUSH"
printf "Verkko:               %d / 3\n" "$VERKKO"
printf "ntLink + GoldRush:    %d / 6\n" "$NTLINK"

echo
echo "Session elapsed:       $(format_seconds "$ELAPSED")"
echo "Rough ETA:             $ETA_TEXT"

echo
echo "=== ACTIVE ASSEMBLER PROCESSES ==="

ps -eo pid,etime,%cpu,%mem,rss,args \
| grep -E \
    '[f]lye |[g]oldrush run|[v]erkko |[n]tLink scaffold|[s]nakemake.*Snakefile.assemblers' \
|| echo "No active assembler process detected."

echo
echo "=== MEMORY ==="
free -h

echo
echo "=== DISK SPACE ==="
df -h "$ROOT" | tail -n 1

echo
echo "=== NEWEST TOOL LOG ==="

LATEST_LOG="$(
    find "$LOGROOT" \
        -type f \
        -name '*.log' \
        ! -name 'dry_run.log' \
        ! -name 'real_run.console.log' \
        -printf '%T@ %p\n' \
        2>/dev/null \
    | sort -nr \
    | head -n 1 \
    | cut -d' ' -f2-
)"

if [[ -n "$LATEST_LOG" ]] && [[ -s "$LATEST_LOG" ]]
then
    echo "$LATEST_LOG"
    echo
    tail -n 12 "$LATEST_LOG"
else
    echo "Waiting for the first assembler log."
    echo
    echo "Latest Snakemake output:"
    tail -n 12 "$CONSOLE" 2>/dev/null
fi

echo
echo "=== FINAL CHECKPOINT ==="

if [[ -s "$CHECKPOINTS/workflow_complete.ok" ]]
then
    echo "COMPLETE: $CHECKPOINTS/workflow_complete.ok"
elif pgrep -af '[s]nakemake.*Snakefile.assemblers' >/dev/null
then
    echo "RUNNING"
elif [[ -s "$CHECKPOINTS/run_exit_status.txt" ]]
then
    STATUS="$(cat "$CHECKPOINTS/run_exit_status.txt")"
    echo "STOPPED — exit status: $STATUS"
else
    echo "Workflow process not detected."
fi

echo
echo "Updated: $(date --iso-8601=seconds)"
echo "============================================================================"
