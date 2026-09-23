#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

DEST="workflows/upstream_reference/legacy_root"

echo "======================================================"
echo " PHASE 2A — SAFE ROOT WORKFLOW MIGRATION"
echo "======================================================"
echo

# ======================================================
# PROTECTED ROOT FILES
# ======================================================
#
# header_assembler.smk:
#   REQUIRED by active whole-genome assembler workflows.
#
# header.smk:
#   Kept at root temporarily for maximum compatibility.
#
# ======================================================

PROTECTED=(
    "header.smk"
    "header_assembler.smk"
)

echo "PROTECTED:"
printf '  %s\n' "${PROTECTED[@]}"

echo


# ======================================================
# Determine which .smk files can move
# ======================================================

mapfile -t FILES < <(
    find . \
        -maxdepth 1 \
        -type f \
        -name '*.smk' \
        ! -name 'header.smk' \
        ! -name 'header_assembler.smk' \
        -printf '%f\n' \
        | sort
)

if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "No movable root .smk files remain."
    exit 0
fi

echo "FILES TO MOVE:"
printf '  %s\n' "${FILES[@]}"

echo
echo "Destination:"
echo "  $DEST"


# ======================================================
# Verify assembler header dependency remains intact
# ======================================================

echo
echo "======================================================"
echo " ACTIVE ASSEMBLER HEADER CHECK"
echo "======================================================"

if [[ ! -f header_assembler.smk ]]; then
    echo "ERROR: header_assembler.smk is missing."
    exit 10
fi

grep -RInF \
    'header_assembler.smk' \
    assemblers \
    --include='*.smk' \
    --include='Snakefile*' \
    --exclude-dir=.snakemake \
    --exclude-dir=results \
    --exclude-dir=debug \
    2>/dev/null || true

echo
echo "header_assembler.smk WILL NOT BE MOVED."


# ======================================================
# Search ACTIVE code for dependencies on files to move
# ======================================================

echo
echo "======================================================"
echo " ACTIVE CODE DEPENDENCY CHECK"
echo "======================================================"

BLOCK=0

for FILE in "${FILES[@]}"; do

    HITS="$(
        grep -RInF "$FILE" \
            assemblers \
            alignment_analysis \
            "SV aligners call/workflow" \
            "SV aligners call/scripts" \
            workflow \
            scripts \
            config \
            tests \
            2>/dev/null \
            --include='*.smk' \
            --include='Snakefile*' \
            --include='*.py' \
            --include='*.sh' \
            --include='*.yaml' \
            --include='*.yml' \
            --exclude-dir=.snakemake \
            --exclude-dir=logs \
            --exclude-dir=results \
            --exclude-dir=debug \
            || true
    )"

    if [[ -n "$HITS" ]]; then

        echo
        echo "DEPENDENCY FOUND FOR:"
        echo "  $FILE"
        echo
        echo "$HITS"

        BLOCK=1

    fi

done

if [[ "$BLOCK" -ne 0 ]]; then

    echo
    echo "======================================================"
    echo "STOPPED"
    echo "An active-code dependency was detected."
    echo "Nothing moved."
    echo "======================================================"

    exit 20
fi

echo "ACTIVE CODE DEPENDENCY CHECK: PASS"


# ======================================================
# FASTQ safety signature
# ======================================================

echo
echo "======================================================"
echo " FASTQ PROTECTION"
echo "======================================================"

FASTQ_BEFORE="$(
    find fastq \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s\n' \
        | sort \
        | sha256sum \
        | awk '{print $1}'
)"

echo "FASTQ BEFORE:"
echo "  $FASTQ_BEFORE"


# ======================================================
# Dry run
# ======================================================

if [[ "${1:-}" != "--apply" ]]; then

    echo
    echo "======================================================"
    echo " DRY RUN"
    echo "======================================================"

    for FILE in "${FILES[@]}"; do
        echo
        echo "$FILE"
        echo "  -> $DEST/$FILE"
    done

    echo
    echo "Protected and staying at root:"
    echo "  header.smk"
    echo "  header_assembler.smk"

    echo
    echo "Nothing moved."
    echo
    echo "Apply with:"
    echo
    echo "  ./phase2a_move_root_workflows.sh --apply"

    exit 0
fi


# ======================================================
# Create destination
# ======================================================

mkdir -p "$DEST"


# ======================================================
# Move files
# ======================================================

echo
echo "======================================================"
echo " MOVING LEGACY ROOT WORKFLOWS"
echo "======================================================"

for FILE in "${FILES[@]}"; do

    if git ls-files --error-unmatch "$FILE" >/dev/null 2>&1; then

        git mv -- \
            "$FILE" \
            "$DEST/$FILE"

    else

        mv -- \
            "$FILE" \
            "$DEST/$FILE"

    fi

    echo "MOVED: $FILE"

done


# ======================================================
# Fix includes to the root header.smk
#
# From:
#
# workflows/upstream_reference/legacy_root/
#
# ../../../header.smk points to repository root.
# ======================================================

echo
echo "======================================================"
echo " FIXING HEADER INCLUDES"
echo "======================================================"

python3 <<'PY'
from pathlib import Path

dest = Path("workflows/upstream_reference/legacy_root")

for f in dest.glob("*.smk"):

    text = f.read_text(errors="replace")

    old = 'include: "header.smk"'
    new = 'include: "../../../header.smk"'

    if old in text:

        text = text.replace(old, new)
        f.write_text(text)

        print(f"UPDATED INCLUDE: {f}")
PY


# ======================================================
# Update README + workflow inventory
# ======================================================

echo
echo "======================================================"
echo " UPDATING DOCUMENTATION PATHS"
echo "======================================================"

python3 - "$DEST" "${FILES[@]}" <<'PY'
from pathlib import Path
import sys

dest = sys.argv[1]
filenames = sys.argv[2:]

targets = [
    Path("README.md"),
    Path("SV aligners call/docs/workflow_tool_inventory.tsv"),
]

for target in targets:

    if not target.exists():
        continue

    text = target.read_text(errors="replace")
    original = text

    for filename in filenames:

        new_path = f"{dest}/{filename}"

        # Avoid repeated replacement if script is ever inspected/reused.
        if new_path not in text:
            text = text.replace(filename, new_path)

    if text != original:

        target.write_text(text)

        print(f"UPDATED: {target}")
PY


# ======================================================
# Verify headers
# ======================================================

echo
echo "======================================================"
echo " HEADER VALIDATION"
echo "======================================================"

test -f header.smk
echo "PASS: header.smk remains at root"

test -f header_assembler.smk
echo "PASS: header_assembler.smk remains at root"


echo
echo "Assembler references:"
grep -RInF \
    'include: "../../header_assembler.smk"' \
    assemblers/whole_genome_asm \
    --include='*.smk' \
    2>/dev/null || true


echo
echo "Moved mapping header references:"
grep -RInF \
    'include: "../../../header.smk"' \
    "$DEST" \
    --include='*.smk' \
    2>/dev/null || true


# ======================================================
# Verify FASTQ again
# ======================================================

FASTQ_AFTER="$(
    find fastq \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s\n' \
        | sort \
        | sha256sum \
        | awk '{print $1}'
)"

echo
echo "FASTQ AFTER:"
echo "  $FASTQ_AFTER"

if [[ "$FASTQ_BEFORE" != "$FASTQ_AFTER" ]]; then

    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "ERROR: FASTQ METADATA CHANGED"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

    exit 99
fi

echo
echo "FASTQ CHECK: PASS"


# ======================================================
# Root check
# ======================================================

echo
echo "======================================================"
echo " ROOT-LEVEL .smk FILES NOW"
echo "======================================================"

find . \
    -maxdepth 1 \
    -type f \
    -name '*.smk' \
    -printf '%f\n' \
    | sort


echo
echo "======================================================"
echo " MOVED LEGACY WORKFLOWS"
echo "======================================================"

find "$DEST" \
    -maxdepth 1 \
    -type f \
    -name '*.smk' \
    -printf '%f\n' \
    | sort


echo
echo "======================================================"
echo " PHASE 2A: PASS"
echo "======================================================"
