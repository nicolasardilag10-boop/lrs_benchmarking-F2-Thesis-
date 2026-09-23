#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking"
CONTAINERS="$PROJECT/containers"
CONFIG="$PROJECT/benchmark_chr21_real/config/containers.yaml"

cd "$PROJECT"

declare -A EXPECTED_COMMANDS=(
    [flye2]="flye"
    [goldrush]="goldrush"
    [ntlink]="ntLink"
    [verkko2]="verkko"
)

declare -A EXPECTED_IMAGES=(
    [flye2]="quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1"
    [goldrush]="quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1"
    [ntlink]="quay.io/biocontainers/ntlink:1.3.11--py313h380480a_2"
    [verkko2]="quay.io/biocontainers/verkko:2.3.2--hb0edd9e_0"
)

FAILURES=0


resolve_base_image()
{
    local dockerfile="$1"

    local arg_image
    local from_value

    arg_image="$(
        sed -nE \
            's/^ARG[[:space:]]+BASE_IMAGE=(.+)$/\1/p' \
            "$dockerfile" \
        | head -n 1
    )"

    from_value="$(
        awk '
            /^FROM[[:space:]]+/ {
                print $2
                exit
            }
        ' "$dockerfile"
    )"

    if [[ "$from_value" == '${BASE_IMAGE}' ]]
    then
        printf '%s\n' "$arg_image"
    else
        printf '%s\n' "$from_value"
    fi
}


config_image_for()
{
    local image_name="$1"

    sed -nE \
        "s|^[[:space:]]+$image_name:[[:space:]]+\"([^\"]+)\"[[:space:]]*$|\1|p" \
        "$CONFIG" \
    | head -n 1
}


echo "============================================================"
echo "STATIC ASSEMBLER CONTAINER VALIDATION"
echo "============================================================"

for IMAGE_NAME in flye2 goldrush ntlink verkko2
do
    DOCKERFILE="$CONTAINERS/$IMAGE_NAME/Dockerfile"
    EXPECTED_COMMAND="${EXPECTED_COMMANDS[$IMAGE_NAME]}"
    EXPECTED_IMAGE="${EXPECTED_IMAGES[$IMAGE_NAME]}"

    echo
    echo "------------------------------------------------------------"
    echo "$IMAGE_NAME"
    echo "------------------------------------------------------------"

    if [[ ! -s "$DOCKERFILE" ]]
    then
        echo "FAIL: Dockerfile is missing or empty."
        FAILURES=$((FAILURES + 1))
        continue
    fi

    echo "PASS: Dockerfile exists."

    RESOLVED_IMAGE="$(resolve_base_image "$DOCKERFILE")"

    if [[ -z "$RESOLVED_IMAGE" ]]
    then
        echo "FAIL: base image could not be resolved."
        FAILURES=$((FAILURES + 1))
    else
        echo "Resolved base image:"
        echo "  $RESOLVED_IMAGE"
    fi

    if [[ "$RESOLVED_IMAGE" == "$EXPECTED_IMAGE" ]]
    then
        echo "PASS: exact tested BioContainer image is used."
    else
        echo "FAIL: unexpected base image."
        echo "Expected:"
        echo "  $EXPECTED_IMAGE"
        echo "Observed:"
        echo "  $RESOLVED_IMAGE"
        FAILURES=$((FAILURES + 1))
    fi

    if [[ "$RESOLVED_IMAGE" == *":latest" ]]
    then
        echo "FAIL: latest tag is not allowed."
        FAILURES=$((FAILURES + 1))
    elif [[ "$RESOLVED_IMAGE" == *:* ]]
    then
        echo "PASS: image has a pinned tag."
    else
        echo "FAIL: image tag is missing."
        FAILURES=$((FAILURES + 1))
    fi

    if grep -Fq \
        "command -v $EXPECTED_COMMAND" \
        "$DOCKERFILE"
    then
        echo "PASS: command validation exists: $EXPECTED_COMMAND"
    else
        echo "FAIL: expected command validation is missing: $EXPECTED_COMMAND"
        FAILURES=$((FAILURES + 1))
    fi

    if grep -Eq \
        '^ENTRYPOINT[[:space:]]+\[\][[:space:]]*$' \
        "$DOCKERFILE"
    then
        echo "PASS: inherited entrypoint is cleared."
    else
        echo "FAIL: ENTRYPOINT [] is missing."
        FAILURES=$((FAILURES + 1))
    fi

    if grep -Eq \
        '^WORKDIR[[:space:]]+/work[[:space:]]*$' \
        "$DOCKERFILE"
    then
        echo "PASS: working directory is /work."
    else
        echo "FAIL: WORKDIR /work is missing."
        FAILURES=$((FAILURES + 1))
    fi
done

echo
echo "------------------------------------------------------------"
echo "CONTAINER CONFIGURATION"
echo "------------------------------------------------------------"

if [[ -s "$CONFIG" ]]
then
    echo "PASS: configuration exists:"
    echo "$CONFIG"
else
    echo "FAIL: container configuration is missing."
    FAILURES=$((FAILURES + 1))
fi

for IMAGE_NAME in flye2 goldrush ntlink verkko2
do
    EXPECTED_URI="docker://${EXPECTED_IMAGES[$IMAGE_NAME]}"
    OBSERVED_URI="$(config_image_for "$IMAGE_NAME")"

    echo
    echo "$IMAGE_NAME"

    if [[ "$OBSERVED_URI" == "$EXPECTED_URI" ]]
    then
        echo "PASS: configuration URI matches Dockerfile."
        echo "  $OBSERVED_URI"
    else
        echo "FAIL: configuration URI does not match."
        echo "Expected:"
        echo "  $EXPECTED_URI"
        echo "Observed:"
        echo "  ${OBSERVED_URI:-MISSING}"
        FAILURES=$((FAILURES + 1))
    fi
done

echo
echo "============================================================"

if [[ "$FAILURES" -eq 0 ]]
then
    echo "STATIC VALIDATION PASSED"
    exit 0
else
    echo "STATIC VALIDATION FAILED: $FAILURES problem(s)"
    exit 1
fi
