#!/usr/bin/env bash

# LRS Benchmarking navigation
# Canonical file:
#   ~/lrs_benchmarking/docs/lrs_shortcuts.sh

_LRS_SHORTCUTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LRS="$(cd "${_LRS_SHORTCUTS_DIR}/.." && pwd)"

# Protected data
export FASTQ="$LRS/fastq"
export REF="$LRS/reference"
export SAMPLES_TRY="$LRS/samples_try"

# Active assembler project
export ASM="$LRS/assemblers"
export WGS="$ASM/whole_genome_asm"
export CHR21="$ASM/benchmark_chr21_real"
export ASMCONFIG="$ASM/config"
export ASMCONT="$ASM/containers"
export ASMSCRIPTS="$ASM/scripts"
export ASMRULES="$ASM/rules"
export ASMRESULTS="$ASM/results"

# Aligner work
export ALIGN="$LRS/final_report_files/snakemake_aligners_benchmarking"
export ALIGNAN="$LRS/alignment_analysis"

# Variant-caller/SV work
export SV="$LRS/SV aligners call"
export SVVALID="$LRS/variant_caller_validation"
export LOCALSV="$LRS/local_sv_validation"

# Organized repository areas
export WORKFLOWS="$LRS/workflows"
export ANALYSIS="$LRS/analysis"
export RESULTS="$LRS/results"
export REPORTS="$LRS/reports"
export NOTES="$LRS/notes"
export DOCS="$LRS/docs"
export TESTS="$LRS/tests"
export ARCHIVE="$LRS/archive"
export LEGACY="$ARCHIVE/legacy"

# Core navigation
alias go-lrs='cd "$LRS"'
alias go-fastq='cd "$FASTQ"'
alias go-ref='cd "$REF"'
alias go-workflows='cd "$WORKFLOWS"'
alias go-legacy='cd "$LEGACY"'
# Compatibility alias retained for older shell habits.
alias go-upstream='cd "$LEGACY"'
alias go-analysis='cd "$ANALYSIS"'
alias go-results='cd "$RESULTS"'
alias go-reports='cd "$REPORTS"'
alias go-notes='cd "$NOTES"'
alias go-docs='cd "$DOCS"'
alias go-tests='cd "$TESTS"'
alias go-archive='cd "$ARCHIVE"'

# Assemblers
alias go-asm='cd "$ASM"'
alias go-wgs='cd "$WGS"'
alias go-chr21='cd "$CHR21"'
alias go-config='cd "$ASMCONFIG"'
alias go-asm-config='cd "$ASMCONFIG"'
alias go-cont='cd "$ASMCONT"'
alias go-asm-cont='cd "$ASMCONT"'
alias go-scripts='cd "$ASMSCRIPTS"'
alias go-asm-scripts='cd "$ASMSCRIPTS"'
alias go-rules='cd "$ASMRULES"'
alias go-asm-results='cd "$ASMRESULTS"'

alias go-flye='cd "$ASMCONT/flye2"'
alias go-goldrush='cd "$ASMCONT/goldrush"'
alias go-ntlink='cd "$ASMCONT/ntlink"'
alias go-verkko='cd "$ASMCONT/verkko2"'

# Aligners
alias go-align='cd "$ALIGN"'
alias go-align-analysis='cd "$ALIGNAN"'

# Variant callers
alias go-sv='cd "$SV"'
alias go-sv-validation='cd "$SVVALID"'
alias go-local-sv='cd "$LOCALSV"'

# Maps
lrs-map() {
    less "$DOCS/PROJECT_MAP.md"
}

lrs-data-map() {
    less "$DOCS/DATA_MAP.md"
}

lrs-workflow-map() {
    less "$DOCS/WORKFLOW_MAP.md"
}

lrs-map-code() {
    code \
        "$DOCS/PROJECT_MAP.md" \
        "$DOCS/DATA_MAP.md" \
        "$DOCS/WORKFLOW_MAP.md"
}

# Path checker
_lrs_path_status() {
    local name="$1"
    local value="$2"

    if [[ -e "$value" ]]; then
        printf '  %-14s [OK]      %s\n' "$name" "$value"
    else
        printf '  %-14s [MISSING] %s\n' "$name" "$value"
    fi
}

lrs-where() {
    printf '\nLRS BENCHMARKING PATHS\n'
    printf '============================================================\n'

    _lrs_path_status "LRS" "$LRS"

    printf '\nPROTECTED DATA\n'
    _lrs_path_status "FASTQ" "$FASTQ"
    _lrs_path_status "REF" "$REF"
    _lrs_path_status "SAMPLES_TRY" "$SAMPLES_TRY"

    printf '\nASSEMBLERS\n'
    _lrs_path_status "ASM" "$ASM"
    _lrs_path_status "WGS" "$WGS"
    _lrs_path_status "CHR21" "$CHR21"
    _lrs_path_status "ASMCONFIG" "$ASMCONFIG"
    _lrs_path_status "ASMCONT" "$ASMCONT"

    printf '\nALIGNERS / VARIANTS\n'
    _lrs_path_status "ALIGN" "$ALIGN"
    _lrs_path_status "ALIGNAN" "$ALIGNAN"
    _lrs_path_status "SV" "$SV"

    printf '\nORGANIZATION\n'
    _lrs_path_status "WORKFLOWS" "$WORKFLOWS"
    _lrs_path_status "LEGACY" "$LEGACY"
    _lrs_path_status "ANALYSIS" "$ANALYSIS"
    _lrs_path_status "RESULTS" "$RESULTS"
    _lrs_path_status "REPORTS" "$REPORTS"
    _lrs_path_status "NOTES" "$NOTES"
    _lrs_path_status "DOCS" "$DOCS"
    _lrs_path_status "TESTS" "$TESTS"
    _lrs_path_status "ARCHIVE" "$ARCHIVE"

    printf '\nCURRENT DIRECTORY\n  %s\n\n' "$PWD"
}

lrs-tree() {
    printf '\nLRS ROOT — COMPACT VIEW\n'
    printf '============================================================\n'

    find "$LRS" \
        -mindepth 1 \
        -maxdepth 2 \
        -type d \
        ! -path "$LRS/.git*" \
        ! -path "$LRS/.snakemake*" \
        ! -path "$LRS/assemblers/results*" \
        ! -path "$LRS/final_report_files/snakemake_aligners_benchmarking/results*" \
        | sed "s#^$LRS/##" \
        | sort

    printf '\n'
}

# Git inspection
lrs-status() {
    git -C "$LRS" status --short
}

lrs-asm-status() {
    git -C "$LRS" status --short -- assemblers
}

lrs-align-status() {
    git -C "$LRS" status --short -- \
        final_report_files/snakemake_aligners_benchmarking
}

lrs-sv-status() {
    git -C "$LRS" status --short -- "SV aligners call"
}

# FASTQ inventory
lrs-fastq() {
    printf '\nFASTQ STORAGE\n'
    printf '============================================================\n'
    printf 'Path: %s\n\n' "$FASTQ"

    find "$FASTQ" \
        -maxdepth 1 \
        -type f \
        -name '*.fastq.gz' \
        -printf '  %f\t%s bytes\n' \
        | sort

    printf '\n'
}

lrs-help() {
cat <<'HELP_EOF'
LRS BENCHMARKING — NAVIGATION
============================================================

CORE
  go-lrs              repository root
  go-fastq            protected FASTQ storage
  go-ref              protected reference storage

ASSEMBLERS
  go-asm              active assembler project
  go-wgs              whole-genome assembler rules
  go-chr21            chromosome-21 benchmark
  go-config           assembler config
  go-cont             assembler containers
  go-scripts          assembler scripts
  go-rules            assembler rules
  go-asm-results      assembler results
  go-flye
  go-goldrush
  go-ntlink
  go-verkko

ALIGNERS
  go-align             main Snakemake aligner benchmark
  go-align-analysis    earlier alignment-analysis workspace

VARIANT CALLERS
  go-sv
  go-sv-validation
  go-local-sv

ORGANIZATION
  go-workflows
  go-legacy             historical workflows and maintenance scripts
  go-upstream           compatibility alias for go-legacy
  go-analysis
  go-results
  go-reports
  go-notes
  go-docs
  go-tests
  go-archive

MAPS
  lrs-map
  lrs-data-map
  lrs-workflow-map
  lrs-map-code

INSPECTION
  lrs-where
  lrs-tree
  lrs-fastq
  lrs-status
  lrs-asm-status
  lrs-align-status
  lrs-sv-status
  lrs-help

MEMORY RULE
  go-X   = navigate
  lrs-X  = inspect / understand

IMPORTANT
  fastq/ and reference/ are protected.
  Active assembler code still lives in assemblers/.
  workflows/assemblers/ is currently only the target organizational area.
HELP_EOF
}

unset _LRS_SHORTCUTS_DIR
