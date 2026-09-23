#!/usr/bin/env bash

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
ROOT="$PROJECT/benchmark_chr21_real"
RESULTS="$ROOT/results/assembler_benchmark"
ASSEMBLIES="$RESULTS/assemblies"
LOGROOT="$RESULTS/logs"
CHECKPOINTS="$RESULTS/checkpoints"

count_assemblies()
{
    local directory="$1"

    find "$directory" \
        -type f \
        -name 'assembly.fasta' \
        -size +0c \
        2>/dev/null \
    | wc -l
}

FLYE="$(count_assemblies "$ASSEMBLIES/flye")"
GOLDRUSH="$(count_assemblies "$ASSEMBLIES/goldrush")"
VERKKO="$(count_assemblies "$ASSEMBLIES/verkko")"
NTLINK="$(count_assemblies "$ASSEMBLIES/ntlink_goldrush")"

COMPLETED=$((FLYE + GOLDRUSH + VERKKO + NTLINK))
TOTAL=21

PERCENT="$(
    awk \
        -v completed="$COMPLETED" \
        -v total="$TOTAL" \
        'BEGIN {
            printf "%.2f", 100 * completed / total
        }'
)"

WIDTH=40
FILLED=$((COMPLETED * WIDTH / TOTAL))
EMPTY=$((WIDTH - FILLED))

LEFT=""
RIGHT=""

printf -v LEFT "%*s" "$FILLED" ""
printf -v RIGHT "%*s" "$EMPTY" ""

LEFT="${LEFT// /#}"
RIGHT="${RIGHT// /-}"

echo "============================================================================"
echo "CHROMOSOME 21 ASSEMBLER BENCHMARK"
echo "============================================================================"
echo
printf "[%s%s] %s%%\n" "$LEFT" "$RIGHT" "$PERCENT"

echo
echo "Completed final assemblies: $COMPLETED / $TOTAL"
echo
printf "Flye:                 %d / 6\n" "$FLYE"
printf "GoldRush:             %d / 6\n" "$GOLDRUSH"
printf "Verkko:               %d / 3\n" "$VERKKO"
printf "ntLink + GoldRush:    %d / 6\n" "$NTLINK"

echo
echo "=== WORKFLOW STATUS ==="

if pgrep -f '[s]nakemake.*Snakefile\.assemblers' >/dev/null
then
    echo "RUNNING: Snakemake is active."
elif [[ -s "$CHECKPOINTS/workflow_complete.ok" ]]
then
    echo "COMPLETE: final workflow checkpoint exists."
elif [[ -s "$CHECKPOINTS/run_exit_status.txt" ]]
then
    STATUS="$(cat "$CHECKPOINTS/run_exit_status.txt")"
    echo "STOPPED: Snakemake exit status $STATUS"
else
    echo "WARNING: Snakemake process not detected."
fi

echo
echo "=== ACTIVE WORKFLOW AND TOOL PROCESSES ==="

ps -eo pid,ppid,stat,etime,%cpu,%mem,rss,args \
| grep -Ei \
    '[s]nakemake.*Snakefile\.assemblers|[g]oldrush|[f]lye|[v]erkko|[n]t[Ll]ink' \
|| echo "No active assembler subprocess detected."

echo
echo "=== CURRENT TOOL ==="

CURRENT_TOOL="$(
    ps -eo args= \
    | grep -Ei \
        '[g]oldrush|[f]lye|[v]erkko|[n]t[Ll]ink' \
    | grep -v 'monitor_chr21_assemblers' \
    | head -n 1
)"

if [[ -n "$CURRENT_TOOL" ]]
then
    case "$CURRENT_TOOL" in
        *goldrush*)
            echo "GoldRush is currently active."
            ;;
        *flye*)
            echo "Flye is currently active."
            ;;
        *verkko*)
            echo "Verkko is currently active."
            ;;
        *ntLink*|*ntlink*)
            echo "ntLink is currently active."
            ;;
        *)
            echo "Assembler subprocess detected."
            ;;
    esac
else
    echo "Snakemake may be preparing or switching between jobs."
fi

echo
echo "=== MEMORY ==="
free -h

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

if [[ -n "$LATEST_LOG" ]] && [[ -f "$LATEST_LOG" ]]
then
    echo "$LATEST_LOG"
    stat \
        -c 'Modified: %y | Size: %s bytes' \
        "$LATEST_LOG"

    echo
    tail -n 15 "$LATEST_LOG"
else
    echo "No individual tool log found yet."

    echo
    echo "Latest Snakemake output:"

    tail -n 15 \
        "$LOGROOT/real_run.console.log" \
        2>/dev/null
fi

echo
echo "Updated: $(date --iso-8601=seconds)"
echo "============================================================================"
