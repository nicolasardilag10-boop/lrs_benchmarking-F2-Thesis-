#!/usr/bin/env python3

import re
import sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COMPLEMENT = str.maketrans(
    "ACGTNacgtn",
    "TGCANtgcan",
)

EXCLUDED_FLAGS = 0x4 | 0x100 | 0x800


def tag_int(fields, name, default):
    prefix = f"{name}:i:"

    for tag in fields[11:]:
        if tag.startswith(prefix):
            try:
                return int(tag[len(prefix):])
            except ValueError:
                return default

    return default


def aligned_query_bases(cigar):
    if cigar == "*":
        return 0

    return sum(
        int(length)
        for length, operation in CIGAR_RE.findall(cigar)
        if operation in {"M", "I", "=", "X"}
    )


def original_orientation_sequence(fields):
    flag = int(fields[1])
    sequence = fields[9]

    if sequence == "*":
        return sequence

    if flag & 0x10:
        return sequence.translate(COMPLEMENT)[::-1]

    return sequence


def ranking_key(fields):
    mapq = int(fields[4])
    aligned_bases = aligned_query_bases(fields[5])
    alignment_score = tag_int(fields, "AS", -10**15)
    edit_distance = tag_int(fields, "NM", 10**15)
    position = int(fields[3])

    return (
        mapq,
        aligned_bases,
        alignment_score,
        -edit_distance,
        -position,
        "\t".join(fields),
    )


def process_group(records):
    global duplicate_groups
    global removed_records

    candidates = []

    for line in records:
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])

        if not (flag & EXCLUDED_FLAGS):
            candidates.append((line, fields))

    if not candidates:
        sys.stderr.write(
            f"ERROR: no mapped primary record for {records[0].split(chr(9), 1)[0]}\n"
        )
        sys.exit(2)

    sequences = {
        original_orientation_sequence(fields)
        for _, fields in candidates
        if fields[9] != "*"
    }

    if len(sequences) > 1:
        qname = candidates[0][1][0]
        sys.stderr.write(
            f"ERROR: conflicting sequences for read name {qname}\n"
        )
        sys.exit(3)

    if len(candidates) > 1:
        duplicate_groups += 1
        removed_records += len(candidates) - 1

    best_line, _ = max(
        candidates,
        key=lambda item: ranking_key(item[1]),
    )

    sys.stdout.write(best_line)


current_qname = None
group = []

duplicate_groups = 0
removed_records = 0

for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    qname = line.split("\t", 1)[0]

    if current_qname is None:
        current_qname = qname

    if qname != current_qname:
        process_group(group)
        group = []
        current_qname = qname

    group.append(line)

if group:
    process_group(group)

sys.stderr.write(
    f"Duplicated QNAME groups processed: {duplicate_groups}\n"
)
sys.stderr.write(
    f"Primary records removed: {removed_records}\n"
)
