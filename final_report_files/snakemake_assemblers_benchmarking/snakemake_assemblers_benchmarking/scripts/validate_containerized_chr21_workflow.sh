#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
ROOT="$PROJECT/benchmark_chr21_real"
SNAKEFILE="$ROOT/Snakefile.assemblers"
CONTAINER_CONFIG="$ROOT/config/containers.yaml"
LOG="$ROOT/results/assembler_benchmark/logs/containerized_dry_run.log"

cd "$PROJECT"

FAILURES=0

echo "============================================================"
echo "CONTAINERIZED CHR21 WORKFLOW VALIDATION"
echo "============================================================"

echo
echo "=== STATIC DIRECTIVE CHECK ==="

CONTAINER_COUNT="$(
    grep -Ec \
        '^[[:space:]]+container:[[:space:]]*$' \
        "$SNAKEFILE" \
    || true
)"

CONDA_COUNT="$(
    grep -Ec \
        '^[[:space:]]+conda:[[:space:]]*$' \
        "$SNAKEFILE" \
    || true
)"

echo "Container rule directives: $CONTAINER_COUNT"
echo "Conda rule directives:     $CONDA_COUNT"

if [[ "$CONTAINER_COUNT" -eq 4 ]]
then
    echo "PASS: four assembler container directives detected."
else
    echo "FAIL: expected four container directives."
    FAILURES=$((FAILURES + 1))
fi

if [[ "$CONDA_COUNT" -eq 0 ]]
then
    echo "PASS: no Conda rule directives remain."
else
    echo "FAIL: Conda directives remain."
    FAILURES=$((FAILURES + 1))
fi

for KEY in flye2 goldrush ntlink verkko2
do
    if grep -Fq "CONTAINERS[\"$KEY\"]" "$SNAKEFILE"
    then
        echo "PASS: Snakefile uses $KEY."
    else
        echo "FAIL: Snakefile does not use $KEY."
        FAILURES=$((FAILURES + 1))
    fi

    if grep -Eq \
        "^[[:space:]]+$KEY:[[:space:]]+\"docker://" \
        "$CONTAINER_CONFIG"
    then
        echo "PASS: configuration defines $KEY."
    else
        echo "FAIL: configuration is missing $KEY."
        FAILURES=$((FAILURES + 1))
    fi
done

echo
echo "=== DRY-RUN ONLY ==="
echo "No software deployment method is enabled."
echo "No container will be pulled or executed."
echo

rm -f "$LOG"

set +e

snakemake \
    --snakefile "$SNAKEFILE" \
    --dry-run \
    --forceall \
    --printshellcmds \
    --cores 4 \
    --resources mem_mb=9500 \
    --nocolor \
    2>&1 | tee "$LOG"

DRY_STATUS="${PIPESTATUS[0]}"

set -e

echo
echo "Dry-run exit status: $DRY_STATUS"

if [[ "$DRY_STATUS" -eq 0 ]]
then
    echo "PASS: Snakemake dry-run succeeded."
else
    echo "FAIL: Snakemake dry-run failed."
    FAILURES=$((FAILURES + 1))
fi

echo
echo "=== EXPECTED JOB COUNTS ==="

declare -A EXPECTED_COUNTS=(
    [validate_inputs]=1
    [flye_assembly]=6
    [goldrush_assembly]=6
    [verkko_assembly]=3
    [ntlink_scaffold]=6
    [workflow_complete]=1
    [all]=1
    [total]=24
)

for RULE in \
    validate_inputs \
    flye_assembly \
    goldrush_assembly \
    verkko_assembly \
    ntlink_scaffold \
    workflow_complete \
    all \
    total
do
    COUNT="${EXPECTED_COUNTS[$RULE]}"

    if grep -Eq \
        "^[[:space:]]*${RULE}[[:space:]]+${COUNT}[[:space:]]*$" \
        "$LOG"
    then
        printf 'PASS: %-22s %s\n' "$RULE" "$COUNT"
    else
        printf 'FAIL: %-22s expected %s\n' "$RULE" "$COUNT"
        FAILURES=$((FAILURES + 1))
    fi
done

echo
echo "=== RUNTIME SAFETY CHECK ==="

if pgrep -f \
    '[g]oldrush|[f]lye|[v]erkko|[n]t[Ll]ink scaffold' \
    >/dev/null
then
    echo "FAIL: an assembler process was unexpectedly detected."

    pgrep -af \
        '[g]oldrush|[f]lye|[v]erkko|[n]t[Ll]ink scaffold'

    FAILURES=$((FAILURES + 1))
else
    echo "PASS: no assembler process was started."
fi

echo
echo "============================================================"

if [[ "$FAILURES" -eq 0 ]]
then
    echo "CONTAINERIZED WORKFLOW VALIDATION PASSED"
    exit 0
else
    echo "CONTAINERIZED WORKFLOW VALIDATION FAILED: $FAILURES problem(s)"
    exit 1
fi
