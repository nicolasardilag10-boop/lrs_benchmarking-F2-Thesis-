#!/usr/bin/env bash

set -Eeuo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"

CONDA="$HOME/miniconda3/bin/conda"

ENV_PREFIX="$PROJECT/.conda/hg002_extract_htslib124"

EXTRACTION_SCRIPT="$PROJECT/benchmark_chr21_real/scripts/extract_HG002_hifi_chr21_chunked.py"

LOGDIR="$PROJECT/benchmark_chr21_real/results/logs/extraction"

CONSOLE_LOG="$LOGDIR/HG002.hifi.chr21.chunked.htslib124.console.log"

TEST_LOG="$LOGDIR/HG002.hifi.window22.htslib124.test.log"

TEST_COUNT="$LOGDIR/HG002.hifi.window22.htslib124.count.txt"

REMOTE_BAM="https://downloads.pacbcloud.com/public/dataset/HG002-CpG-methylation-202202/HG002.GRCh38.haplotagged.bam"

TEST_REGION="chr21:21000001-22000000"

mkdir -p "$PROJECT/.conda"
mkdir -p "$LOGDIR"

echo "============================================================"
echo "STAGE 1/5: CREATE SAMTOOLS 1.24 ENVIRONMENT"
echo "============================================================"

RECREATE_ENV=0

if [[ ! -x "$ENV_PREFIX/bin/samtools" ]]; then
    RECREATE_ENV=1
else
    INSTALLED_VERSION="$(
        "$ENV_PREFIX/bin/samtools" --version |
        awk 'NR == 1 {print $2}'
    )"

    if [[ "$INSTALLED_VERSION" != "1.24" ]]; then
        RECREATE_ENV=1
    fi
fi

if [[ "$RECREATE_ENV" -eq 1 ]]; then
    rm -rf "$ENV_PREFIX"

    "$CONDA" create \
        --yes \
        --prefix "$ENV_PREFIX" \
        --override-channels \
        --strict-channel-priority \
        --channel conda-forge \
        --channel bioconda \
        python=3.11 \
        samtools=1.24
else
    echo "Existing Samtools 1.24 environment will be reused."
fi

echo
echo "============================================================"
echo "STAGE 2/5: VERIFY EXACT VERSIONS"
echo "============================================================"

"$ENV_PREFIX/bin/python" --version
"$ENV_PREFIX/bin/samtools" --version | head -n 5

SAMTOOLS_VERSION="$(
    "$ENV_PREFIX/bin/samtools" --version |
    awk 'NR == 1 {print $2}'
)"

HTSLIB_VERSION="$(
    "$ENV_PREFIX/bin/samtools" --version |
    awk '/Using htslib/ {print $3; exit}'
)"

if [[ "$SAMTOOLS_VERSION" != "1.24" ]]; then
    echo "ERROR: Expected Samtools 1.24, found $SAMTOOLS_VERSION."
    exit 1
fi

if [[ "$HTSLIB_VERSION" != "1.24" ]]; then
    echo "ERROR: Expected HTSlib 1.24, found $HTSLIB_VERSION."
    exit 1
fi

echo "Samtools validation: OK"
echo "HTSlib validation: OK"

echo
echo "============================================================"
echo "STAGE 3/5: ENABLE MID-STREAM HTTP RECOVERY"
echo "============================================================"

export PATH="$ENV_PREFIX/bin:$PATH"

export HTS_RETRY_MAX=30
export HTS_RETRY_DELAY=2000
export HTS_RETRY_MAX_DELAY=60000
export HTS_LOW_SPEED_LIMIT=1024
export HTS_LOW_SPEED_TIME=120

echo "samtools executable: $(command -v samtools)"
echo "HTS_RETRY_MAX=$HTS_RETRY_MAX"
echo "HTS_RETRY_DELAY=$HTS_RETRY_DELAY milliseconds"
echo "HTS_RETRY_MAX_DELAY=$HTS_RETRY_MAX_DELAY milliseconds"
echo "HTS_LOW_SPEED_LIMIT=$HTS_LOW_SPEED_LIMIT bytes/second"
echo "HTS_LOW_SPEED_TIME=$HTS_LOW_SPEED_TIME seconds"

echo
echo "============================================================"
echo "STAGE 4/5: TEST THE PREVIOUSLY FAILING WINDOW"
echo "============================================================"

rm -f "$TEST_LOG"
rm -f "$TEST_COUNT"

echo "Testing: $TEST_REGION"
echo "This may take several minutes."
echo

if ! samtools view \
    -c \
    "$REMOTE_BAM" \
    "$TEST_REGION" \
    > "$TEST_COUNT" \
    2> >(tee "$TEST_LOG" >&2)
then
    echo
    echo "ERROR: Window 22 still could not be read."
    echo "Test log: $TEST_LOG"
    exit 1
fi

WINDOW_COUNT="$(cat "$TEST_COUNT")"

echo
echo "Window 22 remote-read test passed."
echo "Alignments found: $WINDOW_COUNT"

echo
echo "============================================================"
echo "STAGE 5/5: RESUME SAVED EXTRACTION"
echo "============================================================"

if [[ ! -f "$EXTRACTION_SCRIPT" ]]; then
    echo "ERROR: Extraction script is missing:"
    echo "$EXTRACTION_SCRIPT"
    exit 1
fi

echo "Existing valid windows will be retained."
echo "The extraction should resume from window 22."
echo "Progress and ETA will appear below."
echo

"$ENV_PREFIX/bin/python" \
    -u \
    "$EXTRACTION_SCRIPT" \
    2>&1 |
tee -a "$CONSOLE_LOG"

echo
echo "============================================================"
echo "HG002 HIFI CHROMOSOME-21 EXTRACTION COMPLETED"
echo "============================================================"

echo "Console log:"
echo "$CONSOLE_LOG"
