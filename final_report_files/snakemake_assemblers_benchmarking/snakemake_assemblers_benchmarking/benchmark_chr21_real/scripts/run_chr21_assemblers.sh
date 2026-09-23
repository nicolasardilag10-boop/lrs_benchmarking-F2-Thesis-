#!/usr/bin/env bash

set -uo pipefail

source "/home/nicolas/miniconda3/etc/profile.d/conda.sh"
conda activate snakemake

PROJECT="/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
ROOT="/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking/benchmark_chr21_real"
SNAKEFILE="/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking/benchmark_chr21_real/Snakefile.assemblers"
LOG="/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking/benchmark_chr21_real/results/assembler_benchmark/logs/real_run.console.log"
STATUS_FILE="/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking/benchmark_chr21_real/results/assembler_benchmark/checkpoints/run_exit_status.txt"

cd "$PROJECT" || exit 1

mkdir -p     "$ROOT/results/assembler_benchmark/logs"     "$ROOT/results/assembler_benchmark/checkpoints"

date +%s     > "$ROOT/results/assembler_benchmark/checkpoints/run_session_started.epoch"

{
    echo
    echo "============================================================"
    echo "REAL CHROMOSOME 21 ASSEMBLER BENCHMARK"
    echo "Started: $(date --iso-8601=seconds)"
    echo "============================================================"
    echo
} | tee -a "$LOG"

snakemake     --snakefile "$SNAKEFILE"     --cores 4     --use-conda     --resources mem_mb=9500     --rerun-incomplete     --latency-wait 60     --printshellcmds     --show-failed-logs     --keep-going     2>&1 | tee -a "$LOG"

STATUS=${PIPESTATUS[0]}

printf '%s\n' "$STATUS"     > "$STATUS_FILE"

{
    echo
    echo "============================================================"
    echo "Snakemake exit status: $STATUS"
    echo "Finished: $(date --iso-8601=seconds)"
    echo "============================================================"
} | tee -a "$LOG"

exit "$STATUS"
