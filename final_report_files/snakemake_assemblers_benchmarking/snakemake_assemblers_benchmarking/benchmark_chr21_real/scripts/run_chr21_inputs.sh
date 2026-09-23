#!/usr/bin/env bash

set -euo pipefail


cd "$(
    dirname "$0"
)/../.."

ROOT="benchmark_chr21_real"

SOURCE_CHECKPOINT="$(
    printf \
        '%s/results/checkpoints/source/all_sources_validated.ok' \
        "$ROOT"
)"

SOURCE_MANIFEST="$(
    printf \
        '%s/config/sources.resolved.tsv' \
        "$ROOT"
)"

[[ -s "$SOURCE_CHECKPOINT" ]] || {
    echo "ERROR: source checkpoint is missing."
    exit 1
}

[[ -s "$SOURCE_MANIFEST" ]] || {
    echo "ERROR: source manifest is missing."
    exit 1
}

mkdir -p \
    "$ROOT/results/logs/workflow" \
    "$ROOT/results/provenance"

date +%s \
    > "$ROOT/results/provenance/current_input_run_started.epoch"

date --iso-8601=seconds \
    > "$ROOT/results/provenance/current_input_run_started.txt"

timestamp="$(date +%Y%m%d_%H%M%S)"

workflow_log="$(
    printf \
        '%s/results/logs/workflow/chr21_inputs_%s.log' \
        "$ROOT" \
        "$timestamp"
)"

echo "============================================================"
echo "MATCHED CHROMOSOME-21 INPUT WORKFLOW"
echo "============================================================"
echo "Started: $(date --iso-8601=seconds)"
echo "Cores: 4"
echo "Maximum simultaneous remote requests: 2"
echo "Log: $workflow_log"
echo
echo "Stopping safely: Ctrl+C once"
echo "Resuming: run this same script again"
echo

set -o pipefail

snakemake \
    --snakefile "$ROOT/Snakefile.inputs" \
    --cores 4 \
    --resources remote_slots=2 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --keep-going \
    --retries 2 \
    --latency-wait 120 \
    --printshellcmds \
    --show-failed-logs \
    2>&1 \
| tee -a "$workflow_log"

status="${PIPESTATUS[0]}"

echo
echo "Snakemake exit status: $status"
echo "Finished: $(date --iso-8601=seconds)"

exit "$status"
