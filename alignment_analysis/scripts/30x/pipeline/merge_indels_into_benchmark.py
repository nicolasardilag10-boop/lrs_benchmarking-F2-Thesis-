#!/usr/bin/env python3
"""Fill legacy insertion/deletion columns from an indel-recovery table."""

from __future__ import annotations

import argparse
import csv
import os
import stat
import tempfile
from pathlib import Path


PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())
TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
DEFAULT_BENCHMARK = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
DEFAULT_INDELS = TABLE_DIR / "derived" / "alignment_benchmark_30x_indel_recovery.tsv"
MISSING = {"", "NA", "N/A", "NAN", "NONE"}


def normalize_technology(value: str) -> str:
    value = value.strip().lower()
    return "pb" if value in {"pb", "pacbio", "pacbio hifi", "hifi"} else value


def row_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("sample", "").strip().upper(),
        normalize_technology(row.get("read_technology", "")),
        row.get("configuration", "").strip().lower(),
    )


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        return list(reader), list(reader.fieldnames)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fill the benchmark's insertions/deletions columns with event "
            "counts from quality_check_aligners_indels.py."
        )
    )
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK)
    parser.add_argument("--indel-table", type=Path, default=DEFAULT_INDELS)
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_BENCHMARK,
        help="Output TSV; defaults to safely updating the canonical table in place.",
    )
    args = parser.parse_args()

    benchmark = args.benchmark.expanduser().resolve()
    indel_table = args.indel_table.expanduser().resolve()
    output = args.out.expanduser().resolve()

    benchmark_rows, fieldnames = read_tsv(benchmark)
    benchmark_mode = stat.S_IMODE(benchmark.stat().st_mode)
    indel_rows, indel_fields = read_tsv(indel_table)

    required_benchmark = {"sample", "read_technology", "configuration", "insertions", "deletions"}
    required_indels = {"sample", "read_technology", "configuration", "insertion_events", "deletion_events"}
    if missing := required_benchmark - set(fieldnames):
        raise ValueError(f"Benchmark lacks columns: {', '.join(sorted(missing))}")
    if missing := required_indels - set(indel_fields):
        raise ValueError(f"Indel table lacks columns: {', '.join(sorted(missing))}")

    recovered: dict[tuple[str, str, str], tuple[str, str]] = {}
    for row in indel_rows:
        insertion = row.get("insertion_events", "").strip()
        deletion = row.get("deletion_events", "").strip()
        if insertion.upper() in MISSING or deletion.upper() in MISSING:
            continue
        key = row_key(row)
        if key in recovered:
            raise ValueError(f"Duplicate recovered indel identity: {key}")
        recovered[key] = (insertion, deletion)

    filled = 0
    unmatched = 0
    for row in benchmark_rows:
        values = recovered.get(row_key(row))
        if values is None:
            unmatched += 1
            continue
        changed = False
        if row.get("insertions", "").strip().upper() in MISSING:
            row["insertions"] = values[0]
            changed = True
        if row.get("deletions", "").strip().upper() in MISSING:
            row["deletions"] = values[1]
            changed = True
        filled += int(changed)

    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(benchmark_rows)
        os.replace(temporary_name, output)
        os.chmod(output, benchmark_mode)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise

    print(f"Benchmark rows: {len(benchmark_rows)}")
    print(f"Rows filled with indel event counts: {filled}")
    print(f"Rows without recovered indel events: {unmatched}")
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
