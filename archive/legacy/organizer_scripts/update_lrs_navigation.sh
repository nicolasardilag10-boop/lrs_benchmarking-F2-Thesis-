#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(git rev-parse --show-toplevel 2>/dev/null || true)"

if [[ -z "$ROOT" ]]; then
    echo "ERROR: run this from inside ~/lrs_benchmarking"
    exit 1
fi

cd "$ROOT"

STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$ROOT/archive/backups/navigation/$STAMP"

mkdir -p "$BACKUP" "$ROOT/docs" "$ROOT/docs/maps" "$ROOT/docs_maps"

echo "============================================================"
echo "LRS NAVIGATION + MAP UPDATE"
echo "============================================================"
echo "Repository: $ROOT"
echo "Backup:     $BACKUP"
echo

# ------------------------------------------------------------
# Backup current navigation/map files
# ------------------------------------------------------------
for file in \
    docs/lrs_shortcuts.sh \
    docs/PROJECT_MAP.md \
    docs/DATA_MAP.md \
    docs/WORKFLOW_MAP.md \
    docs_maps/lrs_shortcuts.sh \
    docs_maps/PROJECT_MAP.md
do
    if [[ -f "$file" ]]; then
        mkdir -p "$BACKUP/$(dirname "$file")"
        cp -a "$file" "$BACKUP/$file"
        echo "BACKUP: $file"
    fi
done

# ------------------------------------------------------------
# Canonical shortcuts
# ------------------------------------------------------------
cat > "$ROOT/docs/lrs_shortcuts.sh" <<'SHORTCUTS_EOF'
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
export UPSTREAM="$WORKFLOWS/upstream_reference/legacy_root"
export ANALYSIS="$LRS/analysis"
export RESULTS="$LRS/results"
export REPORTS="$LRS/reports"
export NOTES="$LRS/notes"
export DOCS="$LRS/docs"
export TESTS="$LRS/tests"
export ARCHIVE="$LRS/archive"

# Core navigation
alias go-lrs='cd "$LRS"'
alias go-fastq='cd "$FASTQ"'
alias go-ref='cd "$REF"'
alias go-workflows='cd "$WORKFLOWS"'
alias go-upstream='cd "$UPSTREAM"'
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
    _lrs_path_status "UPSTREAM" "$UPSTREAM"
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
  go-upstream
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
SHORTCUTS_EOF

# ------------------------------------------------------------
# PROJECT_MAP.md
# ------------------------------------------------------------
cat > "$ROOT/docs/PROJECT_MAP.md" <<'PROJECT_MAP_EOF'
# LRS Benchmarking — Project Map

**Updated:** 21 August 2026

This map describes the current repository during migration.

## Mental model

```text
~/lrs_benchmarking
│
├── PROTECTED DATA
│   ├── fastq/
│   └── reference/
│
├── ACTIVE WORK
│   ├── assemblers/
│   ├── final_report_files/snakemake_aligners_benchmarking/
│   └── SV aligners call/
│
├── ORGANIZED AREAS
│   ├── workflows/
│   ├── analysis/
│   ├── results/
│   ├── reports/
│   ├── notes/
│   ├── tests/
│   ├── docs/
│   └── archive/
│
└── LEGACY / MIGRATION AREAS
    ├── alignment_analysis/
    ├── aligmnt_analysis/
    ├── final_report_files/
    ├── docs_maps/
    ├── assemblies/
    ├── happy_results/
    ├── truvari/
    └── vcf_called/
```

## Active assemblers

```text
assemblers/
├── Snakefile
├── benchmark_chr21_real/
├── whole_genome_asm/
├── config/
├── containers/
├── rules/
├── scripts/
├── results/
└── debug/
```

Shortcuts:

```bash
go-asm
go-wgs
go-chr21
go-config
go-cont
go-rules
```

`header_assembler.smk` remains at repository root because active WGS
rules still include it through relative paths.

## Aligner benchmark

Main current benchmark:

```text
final_report_files/snakemake_aligners_benchmarking/
```

Shortcut:

```bash
go-align
```

Earlier analysis workspace:

```text
alignment_analysis/
```

Shortcut:

```bash
go-align-analysis
```

## Variant/SV workspace

Primary current workspace:

```text
SV aligners call/
```

Related areas:

```text
variant_caller_validation/
local_sv_validation/
truvari/
sawfish/
vcf_called/
```

Shortcut:

```bash
go-sv
```

## Legacy upstream rules

Root-level legacy Snakemake rules moved during Phase 2A now live in:

```text
workflows/upstream_reference/legacy_root/
```

Shortcut:

```bash
go-upstream
```

## Migration rule

```text
dependency audit
    ↓
path refactor
    ↓
Snakemake dry-run
    ↓
move
    ↓
Snakemake dry-run again
```

Never move active workflow code only for visual tidiness.
PROJECT_MAP_EOF

# ------------------------------------------------------------
# DATA_MAP.md
# ------------------------------------------------------------
cat > "$ROOT/docs/DATA_MAP.md" <<'DATA_MAP_EOF'
# LRS Benchmarking — Data Map

**Updated:** 21 August 2026

## Protected FASTQ anchor

```text
~/lrs_benchmarking/fastq/
```

Shortcut:

```bash
go-fastq
```

This directory must not be moved or renamed during repository cleanup.

Several WGS assembler rules still construct paths using repository-root
`fastq/`.

## Protected reference anchor

```text
~/lrs_benchmarking/reference/
```

Shortcut:

```bash
go-ref
```

Treat this directory as protected until all consuming workflows have been
audited.

## Whole-genome assembler relationship

```text
fastq/
  │
  ├── ONT 30x
  └── PacBio/HiFi 30x
          │
          ▼
assemblers/whole_genome_asm/
          │
          ├── Flye
          ├── GoldRush
          ├── ntLink
          └── Verkko
```

## Chromosome-21 relationship

The chr21 workflow is more portable because input locations are controlled
through configuration / `input_root`.

```text
configured input_root
    │
    ├── HG002/
    ├── HG003/
    └── HG004/
          │
          ▼
assemblers/benchmark_chr21_real/
```

## samples_try

```text
samples_try/
```

Still referenced by older alignment/SV validation code. Do not rename or
archive it until those references have been refactored.

## Generated FASTQ-like files

FASTQ/FQ files inside `results/`, `work/`, `debug/`, or assembler output
directories may be generated intermediates. They are not automatically
authoritative source reads.

Use this distinction:

```text
SOURCE DATA
  fastq/

CONFIGURED TEST INPUTS
  samples_try/ or input_root-controlled locations

GENERATED INTERMEDIATES
  results/ work/ debug/
```
DATA_MAP_EOF

# ------------------------------------------------------------
# WORKFLOW_MAP.md
# ------------------------------------------------------------
cat > "$ROOT/docs/WORKFLOW_MAP.md" <<'WORKFLOW_MAP_EOF'
# LRS Benchmarking — Workflow Map

**Updated:** 21 August 2026

## Whole-genome assembly

```text
fastq/
  │
  ▼
header_assembler.smk
  │
  ▼
assemblers/whole_genome_asm/
  │
  ├── Flye ───────────────► assembly
  │
  ├── GoldRush ───────────► draft assembly
  │                            │
  │                            │ + long reads
  │                            ▼
  │                          ntLink
  │                            │
  │                            ▼
  │                         scaffolds
  │
  └── Verkko ─────────────► assembly
```

## Chromosome-21 assembly

```text
source data
   │
   ▼
chr21 extraction
   │
   ▼
primary + MAPQ >= 20
   │
   ▼
coverage normalization
   │
   ▼
HG002/HG003/HG004 × ONT/HiFi FASTQ
   │
   ▼
assemblers/benchmark_chr21_real/
   │
   ├── Flye
   ├── GoldRush ─► ntLink
   └── Verkko
```

## Aligner benchmark

```text
FASTQ
  │
  ▼
minimap2 / pbmm2 / VACmap / vg Giraffe
  │
  ▼
BAM
  │
  ▼
alignment metrics
  │
  ▼
summary tables
  │
  ▼
statistical comparisons
  │
  ▼
figures / report
```

Current implementation:

```text
final_report_files/snakemake_aligners_benchmarking/
```

## Variant calling / SV

```text
FASTQ
  │
  ▼
alignment
  │
  ▼
BAM
  │
  ▼
Sniffles2 / cuteSV / pbsv / Sawfish / Clair3 / DeepVariant
  │
  ▼
VCF
  │
  ▼
Truvari / validation / benchmarking
```

Primary current workspace:

```text
SV aligners call/
```

## Repository migration

```text
CURRENT ACTIVE CODE
├── assemblers/
├── final_report_files/snakemake_aligners_benchmarking/
└── SV aligners call/
        │
        ▼
dependency + path refactoring
        │
        ▼
FUTURE ORGANIZED CODE
└── workflows/
    ├── assemblers/
    ├── aligners/
    ├── variant_callers/
    └── upstream_reference/
```
WORKFLOW_MAP_EOF

# ------------------------------------------------------------
# Compatibility shim for old docs_maps source line
# ------------------------------------------------------------
cat > "$ROOT/docs_maps/lrs_shortcuts.sh" <<'COMPAT_EOF'
#!/usr/bin/env bash

# Compatibility shim.
# Canonical shortcuts now live at:
#   ~/lrs_benchmarking/docs/lrs_shortcuts.sh

_LRS_COMPAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_LRS_COMPAT_DIR}/../docs/lrs_shortcuts.sh"
unset _LRS_COMPAT_DIR
COMPAT_EOF

cat > "$ROOT/docs_maps/PROJECT_MAP.md" <<'COMPAT_MAP_EOF'
# Compatibility pointer

The canonical project map is now:

```text
../docs/PROJECT_MAP.md
```

Use:

```bash
lrs-map
```

The `docs_maps/` directory is retained temporarily only for compatibility.
COMPAT_MAP_EOF

chmod +x "$ROOT/docs/lrs_shortcuts.sh"
chmod +x "$ROOT/docs_maps/lrs_shortcuts.sh"

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------
echo
echo "VALIDATION"
echo "------------------------------------------------------------"

bash -n "$ROOT/docs/lrs_shortcuts.sh"
echo "PASS: docs/lrs_shortcuts.sh"

bash -n "$ROOT/docs_maps/lrs_shortcuts.sh"
echo "PASS: docs_maps/lrs_shortcuts.sh"

for file in \
    docs/PROJECT_MAP.md \
    docs/DATA_MAP.md \
    docs/WORKFLOW_MAP.md
do
    test -s "$file"
    echo "PASS: $file"
done

test -d "$ROOT/fastq"
echo "PASS: protected fastq/ still exists"

test -d "$ROOT/reference"
echo "PASS: protected reference/ still exists"

echo
echo "============================================================"
echo "NAVIGATION UPDATE: PASS"
echo "============================================================"
echo
echo "Canonical shortcuts:"
echo "  $ROOT/docs/lrs_shortcuts.sh"
echo
echo "Canonical maps:"
echo "  $ROOT/docs/PROJECT_MAP.md"
echo "  $ROOT/docs/DATA_MAP.md"
echo "  $ROOT/docs/WORKFLOW_MAP.md"
echo
echo "Backup:"
echo "  $BACKUP"
echo
echo "Activate now with:"
echo "  source \"$ROOT/docs/lrs_shortcuts.sh\""
