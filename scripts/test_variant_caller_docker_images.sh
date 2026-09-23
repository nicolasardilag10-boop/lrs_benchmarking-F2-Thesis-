#!/usr/bin/env bash

set -uo pipefail

REPORT_DIR="container_validation/docker_inventory"
SUMMARY="$REPORT_DIR/summary.tsv"
DETAIL_LOG="$REPORT_DIR/details.log"

mkdir -p "$REPORT_DIR"

printf "category\timage\tpull_status\tinspect_status\truntime_status\tnotes\n" \
  > "$SUMMARY"

: > "$DETAIL_LOG"

test_image() {
    local category="$1"
    local image="$2"
    local test_command="$3"
    local notes="$4"

    local pull_status="FAILED"
    local inspect_status="NOT_RUN"
    local runtime_status="NOT_RUN"

    {
        echo
        echo "============================================================"
        echo "CATEGORY: $category"
        echo "IMAGE: $image"
        echo "START: $(date --iso-8601=seconds)"
        echo "============================================================"
    } | tee -a "$DETAIL_LOG"

    if timeout 30m docker pull "$image" 2>&1 | tee -a "$DETAIL_LOG"; then
        pull_status="PASSED"
    else
        if [[ "$image" == storage-node:5000/* ]]; then
            pull_status="BLOCKED_INTERNAL_REGISTRY"
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$category" \
          "$image" \
          "$pull_status" \
          "$inspect_status" \
          "$runtime_status" \
          "$notes" >> "$SUMMARY"

        echo "END: $(date --iso-8601=seconds)" | tee -a "$DETAIL_LOG"
        return
    fi

    if docker image inspect "$image" >/dev/null 2>&1; then
        inspect_status="PASSED"

        docker image inspect "$image" \
          --format 'ID={{.Id}} Created={{.Created}} OS={{.Os}} Architecture={{.Architecture}}' \
          2>&1 | tee -a "$DETAIL_LOG"
    else
        inspect_status="FAILED"
    fi

    if timeout 5m docker run --rm \
        --cpus 2 \
        -m 4g \
        "$image" \
        /bin/bash -lc "$test_command" \
        2>&1 | tee -a "$DETAIL_LOG"
    then
        runtime_status="PASSED"
    else
        runtime_status="FAILED_OR_UNSUPPORTED"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$category" \
      "$image" \
      "$pull_status" \
      "$inspect_status" \
      "$runtime_status" \
      "$notes" >> "$SUMMARY"

    echo "END: $(date --iso-8601=seconds)" | tee -a "$DETAIL_LOG"
}

test_image \
  "SNV Clair3" \
  "hkubal/clair3:v1.1.2" \
  'run_clair3.sh --help || /opt/bin/run_clair3.sh --help || clair3 --help' \
  "Public pinned image"

test_image \
  "SNV DeepVariant CPU" \
  "google/deepvariant:1.9.0" \
  '/opt/deepvariant/bin/run_deepvariant --help || /opt/deepvariant/bin/make_examples --help' \
  "Public pinned CPU image"

test_image \
  "SNV DeepVariant GPU" \
  "google/deepvariant:1.9.0-gpu" \
  '/opt/deepvariant/bin/call_variants --help || /opt/deepvariant/bin/run_deepvariant --help' \
  "Image pull can pass without GPU; full GPU runtime requires NVIDIA support"

test_image \
  "Shared utilities" \
  "storage-node:5000/own/genetic_data_analysis:3.0" \
  'samtools --version && bcftools --version' \
  "Internal registry; may only be reachable from GENMEDBFX or institutional network"

{
    echo
    echo "============================================================"
    echo "NON-DOCKER ENTRIES"
    echo "============================================================"
    echo "cuteSV main caller: NOT_APPLICABLE — Mamba environment"
    echo "pbsv main caller: NOT_APPLICABLE — Mamba environment"
    echo "Sniffles2 main caller: NOT_APPLICABLE — Mamba environment"
    echo "Sawfish main caller: NOT_APPLICABLE — Mamba environment"
    echo "Truvari rules: NOT_APPLICABLE — Mamba environment"
    echo "Native bgzip/tabix/samtools entries: NOT_APPLICABLE — no container reference"
} | tee -a "$DETAIL_LOG"

echo
echo "Validation complete."
echo "Summary: $SUMMARY"
echo "Details: $DETAIL_LOG"
