from __future__ import annotations

import sys


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: filter_sam_window.py "
        "CHROMOSOME START END EXCLUDED_FLAGS MIN_MAPQ"
    )


chromosome = sys.argv[1]
window_start = int(sys.argv[2])
window_end = int(sys.argv[3])
excluded_flags = int(sys.argv[4])
minimum_mapq = int(sys.argv[5])


for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 11:
        continue

    flag = int(fields[1])
    reference = fields[2]
    position = int(fields[3])
    mapping_quality = int(fields[4])

    # A regional query can return reads beginning before its
    # requested interval. Assign each primary alignment only
    # to the window containing its alignment start position.
    if reference != chromosome:
        continue

    if not window_start <= position <= window_end:
        continue

    if flag & excluded_flags:
        continue

    if mapping_quality < minimum_mapq:
        continue

    sys.stdout.write(line)
