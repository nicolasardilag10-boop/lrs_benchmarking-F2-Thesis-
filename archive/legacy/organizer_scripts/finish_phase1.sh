#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

echo "========================================"
echo "FINISHING REPOSITORY ORGANIZATION PHASE 1"
echo "========================================"

# ------------------------------------------------------------
# 1. Find the latest migration snapshot
# ------------------------------------------------------------

SNAPSHOT="$(
    find archive/migration_snapshots \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        2>/dev/null \
        | sort \
        | tail -1
)"

if [[ -z "$SNAPSHOT" ]]; then
    echo "ERROR: no migration snapshot found."
    exit 1
fi

echo
echo "Using snapshot:"
echo "  $SNAPSHOT"


# ------------------------------------------------------------
# 2. Verify root FASTQ files
# ------------------------------------------------------------

echo
echo "=== FASTQ SAFETY CHECK ==="

if [[ ! -f "$SNAPSHOT/fastq_metadata_before.tsv" ]]; then
    echo "ERROR: FASTQ metadata snapshot not found."
    exit 2
fi

find fastq \
    -maxdepth 1 \
    -type f \
    -printf '%f\t%s bytes\n' \
    | sort \
    > /tmp/lrs_fastq_metadata_now.tsv

if diff -q \
    "$SNAPSHOT/fastq_metadata_before.tsv" \
    /tmp/lrs_fastq_metadata_now.tsv \
    >/dev/null
then
    echo "FASTQ CHECK: PASS"
    echo "Root FASTQ filenames and sizes are unchanged."
else
    echo
    echo "FASTQ CHECK: FAIL"
    echo "A difference was detected."
    echo "NO FURTHER ORGANIZATION WILL BE PERFORMED."
    echo
    diff -u \
        "$SNAPSHOT/fastq_metadata_before.tsv" \
        /tmp/lrs_fastq_metadata_now.tsv || true
    exit 99
fi


# ------------------------------------------------------------
# 3. Make sure organizational directories exist
# ------------------------------------------------------------

mkdir -p \
    workflows/aligners \
    workflows/assemblers \
    workflows/variant_callers \
    workflows/upstream_reference \
    analysis/aligners \
    analysis/assemblers \
    analysis/variants \
    analysis/tool_selection \
    analysis/notebooks \
    results/alignment \
    results/assembly \
    results/variants \
    results/benchmarking \
    reports/f2 \
    reports/figures \
    reports/tables \
    notes/daily \
    notes/aligners \
    notes/assemblers \
    notes/variant_callers \
    notes/debugging \
    notes/notebooks/experiments \
    notes/notebooks/exploration \
    notes/notebooks/tutorials \
    tests/smoke \
    tests/fixtures \
    tests/container_validation \
    tests/workflow_validation \
    docs/maps \
    docs/architecture \
    docs/workflows \
    docs/containers \
    docs/troubleshooting \
    docs/governance \
    archive/backups \
    archive/installers \
    archive/tutorials \
    archive/legacy \
    archive/external \
    archive/unclassified \
    archive/legacy/organizer_scripts


# ------------------------------------------------------------
# 4. Notes documentation
# ------------------------------------------------------------

cat > notes/README.md <<'EOF'
# Research Notes

Working notes for the long-read benchmarking project.

## Areas

- `daily/` — daily internship/research notes
- `aligners/` — minimap2, pbmm2, VACmap, vg giraffe
- `assemblers/` — Flye, GoldRush, Verkko, ntLink
- `variant_callers/` — structural and small variant callers
- `debugging/` — errors, causes, fixes, and lessons
- `notebooks/` — experiments, exploration, tutorials

Large sequencing files do not belong in this directory.
EOF


# ------------------------------------------------------------
# 5. Migration status
# ------------------------------------------------------------

cat > docs/MIGRATION_STATUS.md <<'EOF'
# Repository Migration Status

The repository is being reorganized incrementally so that active
Snakemake workflows are not broken.

## Protected data

These locations are not moved during the migration:

- `fastq/`
- `reference/`

## Temporarily location-sensitive

These currently require dependency/path review before relocation:

- `assemblers/`
- `alignment_analysis/`
- `SV aligners call/`
- `samples_try/`
- `final_report_files/`
- root-level Snakemake `.smk` files
- generated benchmark/result directories

## Target architecture

```text
lrs_benchmarking/
├── README.md
├── fastq/
├── reference/
├── workflows/
├── analysis/
├── results/
├── reports/
├── notes/
├── tests/
├── docs/
└── archive/
