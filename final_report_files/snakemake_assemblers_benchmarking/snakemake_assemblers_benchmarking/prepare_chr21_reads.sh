#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$PWD"

REFERENCE="$HOME/lrs_benchmarking/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

INPUT_SAMPLES="${PROJECT_DIR}/samples.tsv"
OUTPUT_DIR="${PROJECT_DIR}/data/reads_chr21"
ALIGNMENT_DIR="${PROJECT_DIR}/results/chr21_read_extraction"
OUTPUT_SAMPLES="${PROJECT_DIR}/samples.chr21.tsv"
COVERAGE_TABLE="${PROJECT_DIR}/results/tables/chr21_read_coverage.tsv"

THREADS=4

mkdir -p "$OUTPUT_DIR"
mkdir -p "$ALIGNMENT_DIR"
mkdir -p "$(dirname "$COVERAGE_TABLE")"

# ---------------------------------------------------------------------------
# Validate required programs and files
# ---------------------------------------------------------------------------

for program in minimap2 samtools gzip awk; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $program" >&2
        exit 1
    fi
done

if [[ ! -s "$REFERENCE" ]]; then
    echo "ERROR: Reference FASTA not found:" >&2
    echo "$REFERENCE" >&2
    exit 1
fi

if [[ ! -s "$INPUT_SAMPLES" ]]; then
    echo "ERROR: samples.tsv not found:" >&2
    echo "$INPUT_SAMPLES" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Create/check the FASTA index
# ---------------------------------------------------------------------------

if [[ ! -s "${REFERENCE}.fai" ]]; then
    echo "Creating reference FASTA index..."
    samtools faidx "$REFERENCE"
fi

# Detect whether the reference calls chromosome 21 "chr21" or "21"
if awk '$1 == "chr21" {found=1} END {exit !found}' "${REFERENCE}.fai"; then
    CHR21="chr21"
elif awk '$1 == "21" {found=1} END {exit !found}' "${REFERENCE}.fai"; then
    CHR21="21"
else
    echo "ERROR: Chromosome 21 was not found in ${REFERENCE}.fai" >&2
    exit 1
fi

CHR21_LENGTH=$(
    awk -v chromosome="$CHR21" '
        $1 == chromosome {
            print $2
            exit
        }
    ' "${REFERENCE}.fai"
)

echo "Chromosome name: $CHR21"
echo "Chromosome length: $CHR21_LENGTH bp"

# ---------------------------------------------------------------------------
# Prepare output tables
# ---------------------------------------------------------------------------

printf "sample\ttechnology\tfastq\n" > "$OUTPUT_SAMPLES"

printf \
    "sample\ttechnology\tchromosome\tchromosome_length_bp\tread_count\ttotal_read_bases\testimated_coverage_x\tfastq\n" \
    > "$COVERAGE_TABLE"

# ---------------------------------------------------------------------------
# Map every dataset and extract primary reads assigned to chromosome 21
# ---------------------------------------------------------------------------

while IFS=$'\t' read -r sample technology fastq; do

    [[ "$sample" == "sample" ]] && continue
    [[ -z "${sample:-}" ]] && continue

    if [[ ! -s "$fastq" ]]; then
        echo "ERROR: Input FASTQ missing or empty: $fastq" >&2
        exit 1
    fi

    case "$technology" in
        ont)
            PRESET="map-ont"
            ;;
        pb)
            PRESET="map-hifi"
            ;;
        *)
            echo "ERROR: Unsupported technology: $technology" >&2
            exit 1
            ;;
    esac

    PREFIX="${ALIGNMENT_DIR}/${sample}.${technology}.chr21"
    BAM="${PREFIX}.sorted.bam"
    CHR21_FASTQ="${OUTPUT_DIR}/${sample}.${technology}.chr21.fastq.gz"
    LOG="${PREFIX}.log"

    echo
    echo "============================================================"
    echo "Processing: $sample $technology"
    echo "Preset: $PRESET"
    echo "Input: $fastq"
    echo "============================================================"

    # Map the source reads to GRCh38
    minimap2 \
        -ax "$PRESET" \
        -t "$THREADS" \
        "$REFERENCE" \
        "$fastq" \
        2> "$LOG" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$BAM" \
        -

    samtools index "$BAM"
    samtools quickcheck "$BAM"

    # Extract chromosome-21 primary alignments only.
    # -F 2308 excludes unmapped, secondary, and supplementary alignments.
    samtools view \
        -b \
        -F 2308 \
        "$BAM" \
        "$CHR21" \
    | samtools fastq \
        -n \
        - \
    | gzip -c \
        > "$CHR21_FASTQ"

    READ_COUNT=$(
        gzip -cd "$CHR21_FASTQ" |
        awk 'END {print NR / 4}'
    )

    TOTAL_BASES=$(
        gzip -cd "$CHR21_FASTQ" |
        awk '
            NR % 4 == 2 {
                bases += length($0)
            }
            END {
                print bases + 0
            }
        '
    )

    COVERAGE=$(
        awk \
            -v bases="$TOTAL_BASES" \
            -v chrlen="$CHR21_LENGTH" \
            'BEGIN {
                if (chrlen > 0) {
                    printf "%.3f", bases / chrlen
                } else {
                    print "0.000"
                }
            }'
    )

    printf \
        "%s\t%s\t%s\n" \
        "$sample" \
        "$technology" \
        "$CHR21_FASTQ" \
        >> "$OUTPUT_SAMPLES"

    printf \
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$technology" \
        "$CHR21" \
        "$CHR21_LENGTH" \
        "$READ_COUNT" \
        "$TOTAL_BASES" \
        "$COVERAGE" \
        "$CHR21_FASTQ" \
        >> "$COVERAGE_TABLE"

    echo "Chromosome-21 reads: $READ_COUNT"
    echo "Chromosome-21 bases: $TOTAL_BASES"
    echo "Estimated coverage: ${COVERAGE}x"

done < "$INPUT_SAMPLES"

echo
echo "============================================================"
echo "Chromosome 21 extraction completed"
echo "============================================================"
echo "Sample sheet:"
echo "$OUTPUT_SAMPLES"
echo
echo "Coverage table:"
echo "$COVERAGE_TABLE"
