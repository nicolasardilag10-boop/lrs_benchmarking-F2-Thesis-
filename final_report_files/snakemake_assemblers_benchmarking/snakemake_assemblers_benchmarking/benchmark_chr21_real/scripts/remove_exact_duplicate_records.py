#!/usr/bin/env python3

import sys

current_qname = None
seen_records = set()
removed = 0

for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    qname = line.split("\t", 1)[0]

    if qname != current_qname:
        current_qname = qname
        seen_records = set()

    if line in seen_records:
        removed += 1
        continue

    seen_records.add(line)
    sys.stdout.write(line)

sys.stderr.write(f"Exact duplicate records removed: {removed}\n")
