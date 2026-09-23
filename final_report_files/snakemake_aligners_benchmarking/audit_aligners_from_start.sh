#!/usr/bin/env bash

set -uo pipefail

REPORT="results/audit"
mkdir -p "$REPORT"

SUMMARY="$REPORT/audit_summary.tsv"
PROBLEMS="$REPORT/audit_problems.txt"

: > "$PROBLEMS"

printf "sample\ttechnology\taligner\tfastq_reads\tbam_primary_reads\tmapped_reads\tmapped_percent\tvalid_cigar\twith_NM\tbases_mapped_cigar\tmismatches\terror_percent\tbam_status\n" \
    > "$SUMMARY"


###############################################################################
# SECTION 1: Verify required workflow files
###############################################################################

echo "=== 1. REQUIRED FILES ==="

for file in \
    config.yaml \
    samples.tsv \
    snakemake_aligners.smk \
    data/reference/genome.fasta
do
    if [[ -e "$file" ]]; then
        echo "OK: $file"
    else
        echo "MISSING: $file" | tee -a "$PROBLEMS"
    fi
done


###############################################################################
# SECTION 2: Verify reference genome
###############################################################################

echo
echo "=== 2. REFERENCE GENOME ==="

if [[ ! -f data/reference/genome.fasta.fai ]]; then
    samtools faidx data/reference/genome.fasta
fi

md5sum data/reference/genome.fasta \
    > "$REPORT/reference.md5"

awk '
{
    total += $2
    contigs++
}
END {
    print "Reference contigs:", contigs
    print "Reference total bases:", total
}' data/reference/genome.fasta.fai


###############################################################################
# SECTION 3: Inspect active aligners and graph configuration
###############################################################################

echo
echo "=== 3. CONFIGURATION ==="

grep -A 10 -n '^active_aligners:' config.yaml \
    > "$REPORT/active_aligners.txt"

grep -nE '^vg_(gbz|min|dist|zipcodes):' config.yaml \
    > "$REPORT/vg_paths.txt" || true

cat "$REPORT/active_aligners.txt"
cat "$REPORT/vg_paths.txt"

if grep -q 'data/graph_test/' "$REPORT/vg_paths.txt"; then
    echo \
      "WARNING: vg Giraffe uses data/graph_test, not a full benchmark graph." \
      | tee -a "$PROBLEMS"
fi


###############################################################################
# SECTION 4: Validate six FASTQ inputs
###############################################################################

echo
echo "=== 4. FASTQ INPUTS ==="

for fastq in data/samples/*.1k.fastq.gz; do
    [[ -e "$fastq" ]] || continue

    name=$(basename "$fastq" .1k.fastq.gz)
    sample=$(echo "$name" | cut -d. -f1)
    technology=$(echo "$name" | cut -d. -f2)

    reads=$(
        gzip -cd "$fastq" |
        awk 'END {print NR/4}'
    )

    bases=$(
        gzip -cd "$fastq" |
        awk 'NR % 4 == 2 {bases += length($0)} END {print bases+0}'
    )

    printf "%-8s %-4s reads=%-5s bases=%s\n" \
        "$sample" "$technology" "$reads" "$bases"

    if [[ "$reads" -ne 1000 ]]; then
        echo \
          "FASTQ_COUNT_FAIL: $fastq contains $reads reads, expected 1000" \
          >> "$PROBLEMS"
    fi
done


###############################################################################
# SECTION 5: Validate every BAM from all four aligners
###############################################################################

echo
echo "=== 5. BAM VALIDATION ==="

for bam in results/bam/*.sorted.bam; do
    [[ -e "$bam" ]] || continue

    filename=$(basename "$bam" .sorted.bam)

    # Skip temporary test or calmd duplicate files
    if [[ "$filename" == *".test"* || "$filename" == *".calmd"* ]]; then
        continue
    fi

    sample=$(echo "$filename" | cut -d. -f1)
    technology=$(echo "$filename" | cut -d. -f2)
    aligner=$(echo "$filename" | cut -d. -f3)

    fastq="data/samples/${sample}.${technology}.1k.fastq.gz"

    if [[ -f "$fastq" ]]; then
        fastq_reads=$(
            gzip -cd "$fastq" |
            awk 'END {print NR/4}'
        )
    else
        fastq_reads="NA"
    fi

    if samtools quickcheck "$bam"; then
        bam_status="OK"
    else
        bam_status="FAILED"
        echo "BAM_QUICKCHECK_FAIL: $bam" >> "$PROBLEMS"
    fi

    sort_order=$(
        samtools view -H "$bam" |
        awk -F'\t' '
        /^@HD/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^SO:/) {
                    sub(/^SO:/, "", $i)
                    print $i
                }
            }
        }'
    )

    if [[ "$sort_order" != "coordinate" ]]; then
        echo \
          "SORT_ORDER_FAIL: $bam has sort order '$sort_order'" \
          >> "$PROBLEMS"
    fi

    primary=$(
        samtools view -c -F 0x900 "$bam"
    )

    mapped=$(
        samtools view -c -F 0x904 "$bam"
    )

    mapped_percent=$(
        awk -v m="$mapped" -v t="$primary" '
        BEGIN {
            if (t > 0) printf "%.6f", 100*m/t
            else print "NA"
        }'
    )

    valid_cigar=$(
        samtools view -F 0x904 "$bam" |
        awk '$6 != "*" {count++} END {print count+0}'
    )

    with_nm=$(
        samtools view -F 0x904 "$bam" |
        awk '
        {
            for (i=12; i<=NF; i++) {
                if ($i ~ /^NM:i:/) {
                    count++
                    break
                }
            }
        }
        END {print count+0}'
    )

    stats="$REPORT/${filename}.stats.txt"
    samtools stats "$bam" > "$stats"

    bases_mapped_cigar=$(
        awk -F'\t' '
        $1=="SN" && $2=="bases mapped (cigar):" {
            gsub(/ .*/, "", $3)
            print $3
        }' "$stats"
    )

    mismatches=$(
        awk -F'\t' '
        $1=="SN" && $2=="mismatches:" {
            gsub(/ .*/, "", $3)
            print $3
        }' "$stats"
    )

    error_rate=$(
        awk -F'\t' '
        $1=="SN" && $2=="error rate:" {
            gsub(/ .*/, "", $3)
            print $3
        }' "$stats"
    )

    error_percent=$(
        awk -v e="$error_rate" '
        BEGIN {
            if (e == "" || e == "NA") print "NA"
            else printf "%.6f", 100*e
        }'
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$technology" \
        "$aligner" \
        "$fastq_reads" \
        "$primary" \
        "$mapped" \
        "$mapped_percent" \
        "$valid_cigar" \
        "$with_nm" \
        "${bases_mapped_cigar:-NA}" \
        "${mismatches:-NA}" \
        "${error_percent:-NA}" \
        "$bam_status" \
        >> "$SUMMARY"

    if [[ "$fastq_reads" != "NA" && "$primary" -ne "$fastq_reads" ]]; then
        echo \
          "READ_COUNT_MISMATCH: $bam FASTQ=$fastq_reads BAM=$primary" \
          >> "$PROBLEMS"
    fi

    if [[ "$mapped" -ne "$valid_cigar" ]]; then
        echo \
          "CIGAR_COUNT_MISMATCH: $bam mapped=$mapped cigar=$valid_cigar" \
          >> "$PROBLEMS"
    fi

    if [[ "$with_nm" -ne 0 && "$with_nm" -ne "$mapped" ]]; then
        echo \
          "PARTIAL_NM_TAGS: $bam mapped=$mapped with_NM=$with_nm" \
          >> "$PROBLEMS"
    fi
done


###############################################################################
# SECTION 6: Check expected number of BAMs per aligner
###############################################################################

echo
echo "=== 6. BAM COUNTS ==="

for aligner in minimap2 pbmm2 vacmap vg_giraffe; do
    count=$(
        awk -F'\t' -v a="$aligner" '
        NR > 1 && $3 == a {count++}
        END {print count+0}
        ' "$SUMMARY"
    )

    echo "$aligner: $count"

    if [[ "$count" -ne 6 ]]; then
        echo \
          "BAM_COUNT_FAIL: $aligner has $count results, expected 6" \
          >> "$PROBLEMS"
    fi
done


###############################################################################
# SECTION 7: Flag known scientific comparability problems
###############################################################################

echo
echo "=== 7. SCIENTIFIC COMPARABILITY ==="

if awk -F'\t' '
NR > 1 && $2=="ont" && $3=="pbmm2" {found=1}
END {exit !found}
' "$SUMMARY"; then
    echo \
      "WARNING: pbmm2 ONT results exist but are not a standard PacBio pbmm2 benchmark." \
      | tee -a "$PROBLEMS"
fi

awk -F'\t' '
NR > 1 && $3=="vacmap" && ($6==0 || $7 < 1) {
    print "VACMAP_SUSPICIOUS:", $1, $2, "mapped_reads=" $6, "mapped_percent=" $7
}
' "$SUMMARY" | tee -a "$PROBLEMS"

if grep -q 'data/graph_test/' "$REPORT/vg_paths.txt"; then
    echo \
      "VG_NOT_DIRECTLY_COMPARABLE: vg Giraffe used a test graph." \
      >> "$PROBLEMS"
fi


###############################################################################
# SECTION 8: Compare regenerated audit values with final summary table
###############################################################################

echo
echo "=== 8. FINAL TABLE PRESENCE ==="

if [[ -f results/tables/alignment_summary.tsv ]]; then
    echo "OK: results/tables/alignment_summary.tsv"

    cp results/tables/alignment_summary.tsv \
       "$REPORT/alignment_summary_workflow.tsv"
else
    echo \
      "MISSING: results/tables/alignment_summary.tsv" \
      | tee -a "$PROBLEMS"
fi


###############################################################################
# FINAL REPORT
###############################################################################

echo
echo "============================================================"
echo "AUDIT COMPLETE"
echo "============================================================"
echo "Summary:  $SUMMARY"
echo "Problems: $PROBLEMS"
echo

if [[ -s "$PROBLEMS" ]]; then
    echo "RESULT: ISSUES FOUND"
    cat "$PROBLEMS"
else
    echo "RESULT: NO STRUCTURAL ISSUES FOUND"
fi
