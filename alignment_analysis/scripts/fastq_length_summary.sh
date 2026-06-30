#!/usr/bin/env bash

set -euo pipefail

SOURCE="/mnt/c/Users/nicol/OneDrive/Documentos/Maestría/Third Semester/F2 Practicals/Material to study/giab(Samples)"

RESULT="$HOME/lrs_benchmarking/alignment_analysis/tables/fastq_length_summary.tsv"

printf "file\tsample\ttechnology\tread_count\ttotal_bases\tmean_length\tmedian_length\tN50\tmin_length\tmax_length\n" > "$RESULT"

for file in "$SOURCE"/*.fastq.gz; do

    temporary_lengths=$(mktemp)

    gzip -cd -- "$file" |
    awk 'NR % 4 == 2 {print length($0)}' \
    > "$temporary_lengths"

    filename=$(basename "$file")

    sample=${filename%%.*}

    remainder=${filename#*.}
    technology=${remainder%%.*}

    read_count=$(
        wc -l < "$temporary_lengths"
    )

    total_bases=$(
        awk '{sum += $1} END {print sum + 0}' \
        "$temporary_lengths"
    )

    mean_length=$(
        awk '
        {
            sum += $1
        }
        END {
            if (NR > 0) {
                printf "%.2f", sum / NR
            }
        }
        ' "$temporary_lengths"
    )

    median_length=$(
        sort -n "$temporary_lengths" |
        awk '
        {
            values[NR] = $1
        }
        END {
            if (NR % 2 == 1) {
                print values[(NR + 1) / 2]
            } else {
                printf "%.2f",
                (values[NR / 2] + values[NR / 2 + 1]) / 2
            }
        }
        '
    )

    n50=$(
        sort -nr "$temporary_lengths" |
        awk '
        {
            values[NR] = $1
            total += $1
        }
        END {
            half = total / 2

            for (i = 1; i <= NR; i++) {
                cumulative += values[i]

                if (cumulative >= half) {
                    print values[i]
                    break
                }
            }
        }
        '
    )

    min_length=$(
        sort -n "$temporary_lengths" |
        head -n 1
    )

    max_length=$(
        sort -nr "$temporary_lengths" |
        head -n 1
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$filename" \
        "$sample" \
        "$technology" \
        "$read_count" \
        "$total_bases" \
        "$mean_length" \
        "$median_length" \
        "$n50" \
        "$min_length" \
        "$max_length" \
        >> "$RESULT"

    rm -f "$temporary_lengths"

done

echo "Saved result: $RESULT"
