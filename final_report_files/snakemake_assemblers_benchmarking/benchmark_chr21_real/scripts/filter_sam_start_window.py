#!/usr/bin/env python3

from __future__ import annotations

import sys


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


if len(sys.argv) != 6:
    fail(
        "Usage: filter_sam_start_window.py "
        "CHROMOSOME START END EXCLUDED_FLAGS MINIMUM_MAPQ"
    )

chromosome = sys.argv[1]
window_start = int(sys.argv[2])
window_end = int(sys.argv[3])
excluded_flags = int(sys.argv[4])
minimum_mapq = int(sys.argv[5])

input_records = 0
selected_records = 0

for line_number, line in enumerate(sys.stdin, start=1):
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    input_records += 1
    fields = line.rstrip("\n").split("\t")

    if len(fields) < 11:
        fail(f"Malformed SAM record at line {line_number}")

    flag = int(fields[1])
    reference_name = fields[2]
    position = int(fields[3])
    mapq = int(fields[4])

    if reference_name != chromosome:
        continue

    # Assign each read to exactly one window according to alignment start.
    # This prevents duplicate reads at window boundaries.
    if not window_start <= position <= window_end:
        continue

    if flag & excluded_flags:
        continue

    if mapq < minimum_mapq:
        continue

    sys.stdout.write(line)
    selected_records += 1

print(
    (
        f"Window {chromosome}:{window_start}-{window_end}: "
        f"input_records={input_records}, "
        f"selected_records={selected_records}"
    ),
    file=sys.stderr,
)
