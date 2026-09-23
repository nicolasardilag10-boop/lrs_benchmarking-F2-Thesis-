#!/usr/bin/env bash

set -Eeuo pipefail

ROOT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
cd "$ROOT"

DATA_DIR="smoke_test/data/verkko_official"
HIFI="$DATA_DIR/ecoli_hifi_subset24x.fastq.gz"
ONT="$DATA_DIR/ecoli_ont_subset50x.fastq.gz"

HIFI_URL="https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz"
ONT_URL="https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz"

TARGET="results/smoke_test/verkko/ecoli/assembly.fasta"
WORKDIR="results/smoke_test/verkko/ecoli"

mkdir -p "$DATA_DIR"

download_fastq() {
    local url="$1"
    local output="$2"

    if [[ -s "$output" ]] && gzip -t "$output" 2>/dev/null; then
        echo "Reusing valid file: $output"
        return
    fi

    rm -f "$output" "${output}.partial"

    echo "Downloading: $output"

    curl \
        --location \
        --fail \
        --retry 5 \
        --retry-delay 3 \
        --output "${output}.partial" \
        "$url"

    gzip -t "${output}.partial"
    mv "${output}.partial" "$output"

    echo "Download verified: $output"
}

echo
echo "============================================================"
echo "1. DOWNLOADING OFFICIAL VERKKO TEST DATA"
echo "============================================================"

download_fastq "$HIFI_URL" "$HIFI"
download_fastq "$ONT_URL" "$ONT"

echo
echo "============================================================"
echo "2. VERIFYING INPUT FILES"
echo "============================================================"

gzip -t "$HIFI"
gzip -t "$ONT"

ls -lh "$HIFI" "$ONT"

echo
echo "FASTQ headers:"

gzip -cd "$HIFI" | head -n 1 || true
gzip -cd "$ONT" | head -n 1 || true

echo
echo "============================================================"
echo "3. CREATING VERKKO-SPECIFIC SAMPLE SHEET"
echo "============================================================"

cat > smoke_test/samples.verkko.tsv <<EOF
sample	technology	fastq
ecoli	pb	$HIFI
ecoli	ont	$ONT
EOF

cat smoke_test/samples.verkko.tsv

echo
echo "============================================================"
echo "4. CREATING VERKKO-SPECIFIC CONFIG"
echo "============================================================"

cat > smoke_test/config.verkko.yaml <<'EOF'
sample_sheet: smoke_test/samples.verkko.tsv

active_assemblers:
  - verkko

genome_size: 5m
threads: 4
EOF

cat smoke_test/config.verkko.yaml

echo
echo "============================================================"
echo "5. VERIFYING VERKKO CONDA ENVIRONMENT"
echo "============================================================"

VERKKO_ENV=""

for env in .snakemake/conda/*_; do
    if [[ -x "$env/bin/verkko" ]]; then
        VERKKO_ENV="$env"
        break
    fi
done

if [[ -z "$VERKKO_ENV" ]]; then
    echo "ERROR: Verkko environment was not found."
    exit 1
fi

echo "Verkko environment: $VERKKO_ENV"

"$VERKKO_ENV/bin/python" --version
"$VERKKO_ENV/bin/snakemake" --version

"$VERKKO_ENV/bin/python" - <<'PY'
import networkx
import sys

print("Python executable:", sys.executable)
print("NetworkX version:", networkx.__version__)
print("Verkko environment validation: OK")
PY

echo
echo "============================================================"
echo "6. REMOVING ONLY THE PREVIOUS E. COLI TEST"
echo "============================================================"

rm -rf "$WORKDIR"

echo
echo "============================================================"
echo "7. SNAKEMAKE DRY RUN"
echo "============================================================"

snakemake \
    --snakefile Snakefile \
    --configfile smoke_test/config.verkko.yaml \
    --use-conda \
    --cores 4 \
    --dry-run \
    "$TARGET"

echo
echo "============================================================"
echo "8. RUNNING OFFICIAL VERKKO TEST"
echo "============================================================"

snakemake \
    --snakefile Snakefile \
    --configfile smoke_test/config.verkko.yaml \
    --use-conda \
    --cores 4 \
    --printshellcmds \
    --rerun-incomplete \
    "$TARGET"

echo
echo "============================================================"
echo "9. VALIDATING FINAL ASSEMBLY"
echo "============================================================"

if [[ ! -s "$TARGET" ]]; then
    echo "ERROR: Final assembly is missing or empty:"
    echo "$TARGET"
    exit 1
fi

python - "$TARGET" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])

lengths = []
current = 0

with path.open() as handle:
    for line in handle:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current:
                lengths.append(current)
            current = 0
        else:
            current += len(line)

if current:
    lengths.append(current)

if not lengths:
    raise SystemExit("ERROR: No sequences found in assembly.")

lengths.sort(reverse=True)
total = sum(lengths)

running = 0
n50 = 0

for length in lengths:
    running += length

    if running >= total / 2:
        n50 = length
        break

print()
print("VERKKO OFFICIAL TEST SUCCESSFUL")
print("--------------------------------")
print(f"Assembly: {path}")
print(f"Contigs: {len(lengths):,}")
print(f"Total bases: {total:,}")
print(f"Longest contig: {max(lengths):,}")
print(f"N50: {n50:,}")
PY

echo
echo "Final file:"
ls -lh "$TARGET"

echo
echo "============================================================"
echo "VERKKO RULE VALIDATED SUCCESSFULLY"
echo "============================================================"
