#!/usr/bin/env bash

set -Eeuo pipefail

# ============================================================
# LRS BENCHMARKING REPOSITORY ORGANIZER — PHASE 1
#
# PURPOSE
# -------
# Safely organize low-risk repository content WITHOUT moving
# sequencing data or active Snakemake workflows.
#
# Default:
#   DRY RUN
#
# Apply:
#   ./organize_repository_phase1.sh --apply
#
# HARD RULE:
#   FASTQ files are NEVER moved by this script.
# ============================================================


# ------------------------------------------------------------
# 1. Resolve repository root
# ------------------------------------------------------------

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"

if [[ -z "$ROOT" ]]; then
    echo "ERROR: This does not appear to be a Git repository."
    exit 1
fi

cd "$ROOT"


# ------------------------------------------------------------
# 2. Select mode
# ------------------------------------------------------------

MODE="dry-run"

if [[ "${1:-}" == "--apply" ]]; then
    MODE="apply"
elif [[ -n "${1:-}" ]]; then
    echo "Usage:"
    echo "  $0"
    echo "  $0 --apply"
    exit 1
fi


echo
echo "============================================================"
echo " LRS BENCHMARKING REPOSITORY ORGANIZER"
echo "============================================================"
echo
echo "Repository:"
echo "  $ROOT"
echo
echo "Mode:"
echo "  $MODE"
echo


# ------------------------------------------------------------
# 3. HARD PROTECTION
# ------------------------------------------------------------
#
# None of these may be moved during Phase 1.
#
# Some are protected because they contain real data.
# Others are protected because current Snakemake workflows
# depend on their location.
# ------------------------------------------------------------

is_protected() {

    local path="$1"

    case "$path" in

        fastq|fastq/*)
            return 0
            ;;

        reference|reference/*)
            return 0
            ;;

        data|data/*)
            return 0
            ;;

        samples_try|samples_try/*)
            return 0
            ;;

        assemblers|assemblers/*)
            return 0
            ;;

        assemblies|assemblies/*)
            return 0
            ;;

        "SV aligners call"|"SV aligners call"/*)
            return 0
            ;;

        alignment_analysis|alignment_analysis/*)
            return 0
            ;;

        final_report_files|final_report_files/*)
            return 0
            ;;

        workflow|workflow/*)
            return 0
            ;;

        results|results/*)
            return 0
            ;;

        vcf_called|vcf_called/*)
            return 0
            ;;

        happy_results|happy_results/*)
            return 0
            ;;

        local_sv_validation|local_sv_validation/*)
            return 0
            ;;

        variant_caller_validation|variant_caller_validation/*)
            return 0
            ;;

        sawfish|sawfish/*)
            return 0
            ;;

        truvari|truvari/*)
            return 0
            ;;

        *.smk)
            return 0
            ;;

        header.smk)
            return 0
            ;;

        header_assembler.smk)
            return 0
            ;;

    esac

    return 1
}


# ------------------------------------------------------------
# 4. Command helpers
# ------------------------------------------------------------

run() {

    if [[ "$MODE" == "dry-run" ]]; then

        printf "DRY RUN: "

        printf "%q " "$@"

        echo

    else

        "$@"

    fi
}


move_path() {

    local src="$1"
    local dst="$2"

    if [[ ! -e "$src" && ! -L "$src" ]]; then
        echo "SKIP: $src does not exist."
        return 0
    fi

    if is_protected "$src"; then

        echo
        echo "============================================================"
        echo "SAFETY BLOCK"
        echo "Attempted to move protected path:"
        echo
        echo "  $src"
        echo
        echo "Nothing was moved."
        echo "============================================================"
        echo

        exit 50
    fi

    if [[ -e "$dst" || -L "$dst" ]]; then
        echo "SKIP: destination already exists:"
        echo "  $dst"
        return 0
    fi

    run mkdir -p "$(dirname "$dst")"

    # Preserve Git rename history when possible.
    if git ls-files -- "$src" | grep -q .; then
        run git mv -- "$src" "$dst"
    else
        run mv -- "$src" "$dst"
    fi
}


# ------------------------------------------------------------
# 5. FASTQ metadata signature
# ------------------------------------------------------------
#
# IMPORTANT:
# This does NOT read FASTQ contents.
#
# It hashes only:
#     filename + file size
#
# This allows us to verify that the organizer did not alter
# the root FASTQ collection.
# ------------------------------------------------------------

fastq_signature() {

    if [[ ! -d "$ROOT/fastq" ]]; then
        echo "FASTQ_DIRECTORY_MISSING"
        return
    fi

    find "$ROOT/fastq" \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s\n' \
        | sort \
        | sha256sum \
        | awk '{print $1}'
}


FASTQ_BEFORE="$(fastq_signature)"


echo "FASTQ metadata signature BEFORE:"
echo "  $FASTQ_BEFORE"
echo


# ------------------------------------------------------------
# 6. Snapshot current repository state
# ------------------------------------------------------------

STAMP="$(date +%Y%m%d_%H%M%S)"

SNAPSHOT_DIR="archive/migration_snapshots/$STAMP"

if [[ "$MODE" == "apply" ]]; then

    mkdir -p "$SNAPSHOT_DIR"

    git status --short > \
        "$SNAPSHOT_DIR/git_status_before.txt"

    git diff > \
        "$SNAPSHOT_DIR/git_diff_before.patch" || true

    git diff --cached > \
        "$SNAPSHOT_DIR/git_diff_cached_before.patch" || true

    find fastq \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s bytes\n' \
        | sort \
        > "$SNAPSHOT_DIR/fastq_metadata_before.tsv"

fi


# ------------------------------------------------------------
# 7. Create clean repository architecture
# ------------------------------------------------------------

echo "Creating repository architecture..."
echo

run mkdir -p \
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
    archive/migration_snapshots


# ------------------------------------------------------------
# 8. Move LOW-RISK documentation / analysis files
# ------------------------------------------------------------

echo
echo "Organizing low-risk project files..."
echo


# Project governance
move_path \
    "CONSTITUTION.md" \
    "docs/governance/CONSTITUTION.md"


# Tool-selection analysis
move_path \
    "lrs_benchmarking_tools_selection.R.Rmd" \
    "analysis/tool_selection/lrs_benchmarking_tools_selection.R.Rmd"


move_path \
    "lrs_benchmarking_tools_selection.xlsx" \
    "analysis/tool_selection/lrs_benchmarking_tools_selection.xlsx"


# Large installer should not live in repository root
move_path \
    "miniconda.sh" \
    "archive/installers/miniconda.sh"


# ------------------------------------------------------------
# 9. Archive obvious root artifacts
# ------------------------------------------------------------

ROOT_ARTIFACTS="archive/unclassified/root_artifacts_${STAMP}"

for artifact in \
    "%1" \
    "FETCH_HEAD" \
    "Final_reports" \
    "echo" \
    "file.sh" \
    "pbmm2" \
    ".DS_Store" \
    ".tmp"
do

    if [[ -e "$artifact" || -L "$artifact" ]]; then

        move_path \
            "$artifact" \
            "$ROOT_ARTIFACTS/$artifact"

    fi

done


# ------------------------------------------------------------
# 10. Add Notes README
# ------------------------------------------------------------

if [[ "$MODE" == "apply" && ! -f notes/README.md ]]; then

cat > notes/README.md <<'EOF'
# Research Notes

Personal working documentation for the long-read benchmarking project.

## Structure

- `daily/` — daily research and internship notes
- `aligners/` — minimap2, pbmm2, VACmap, vg giraffe
- `assemblers/` — Flye, GoldRush, Verkko, ntLink
- `variant_callers/` — SV and small-variant caller notes
- `debugging/` — problems, errors, explanations, and solutions
- `notebooks/`
  - `experiments/`
  - `exploration/`
  - `tutorials/`

Notes are not workflow outputs and should not contain large sequencing data.
EOF

fi


# ------------------------------------------------------------
# 11. Create migration status document
# ------------------------------------------------------------

if [[ "$MODE" == "apply" ]]; then

cat > docs/MIGRATION_STATUS.md <<'EOF'
# Repository Migration Status

This repository is being reorganized incrementally to avoid breaking
working Snakemake pipelines.

## Protected during Phase 1

The following locations MUST NOT be moved yet:

- `fastq/`
- `reference/`
- `data/`
- `samples_try/`
- `assemblers/`
- `assemblies/`
- `alignment_analysis/`
- `SV aligners call/`
- `final_report_files/`
- `workflow/`
- current root-level `.smk` workflows
- workflow result directories

## Reason

Some current workflows use location-dependent paths, including
repository-root `fastq/` references.

They will be migrated only after their path handling has been
refactored and validated with Snakemake dry-runs.

## Intended final logical structure

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
