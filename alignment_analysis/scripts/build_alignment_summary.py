#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


PROJECT = Path.home() / "lrs_benchmarking"
TABLES = PROJECT / "alignment_analysis" / "tables"
OUTPUT = TABLES / "alignment_summary.tsv"

# Accept both standardized ".1k" files and older ".test" files.
PATTERN = re.compile(
    r"^(HG00[234])\.(ont|pb)\.(1k|test)\.samtools_stats\.txt$"
)

EXPECTED = {
    ("HG002", "ont"),
    ("HG002", "pb"),
    ("HG003", "ont"),
    ("HG003", "pb"),
    ("HG004", "ont"),
    ("HG004", "pb"),
}


def parse_sn_file(path: Path) -> dict[str, str]:
    """Read the SN summary section from one samtools stats file."""

    values: dict[str, str] = {}

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                continue

            name = fields[1].strip().rstrip(":")
            value = fields[2].strip()

            values[name] = value

    return values


def get_number(
    statistics: dict[str, str],
    name: str,
    *,
    integer: bool = False,
    default=None,
):
    """Return a numeric SN value, or a default when it is absent."""

    raw_value = statistics.get(name)

    if raw_value is None:
        return default

    number_text = raw_value.split()[0].replace(",", "")

    try:
        number = float(number_text)
    except ValueError:
        return default

    return int(number) if integer else number


def choose_files() -> dict[tuple[str, str], Path]:
    """
    Find the six reports.

    When both '.1k' and '.test' exist for the same dataset,
    prefer the standardized '.1k' file.
    """

    selected: dict[tuple[str, str], Path] = {}

    for path in sorted(TABLES.glob("*.samtools_stats.txt")):
        match = PATTERN.match(path.name)

        if match is None:
            continue

        sample = match.group(1)
        technology = match.group(2)
        suffix = match.group(3)

        dataset = (sample, technology)

        if dataset not in selected:
            selected[dataset] = path
        elif suffix == "1k":
            selected[dataset] = path

    return selected


def build_row(path: Path, sample: str, technology: str) -> dict[str, object]:
    """Create one combined-table row."""

    stats = parse_sn_file(path)

    raw_total = get_number(
        stats,
        "raw total sequences",
        integer=True,
    )

    reads_mapped = get_number(
        stats,
        "reads mapped",
        integer=True,
    )

    if raw_total is None:
        raise ValueError("missing 'raw total sequences'")

    if reads_mapped is None:
        raise ValueError("missing 'reads mapped'")

    reads_unmapped = get_number(
        stats,
        "reads unmapped",
        integer=True,
        default=raw_total - reads_mapped,
    )

    total_length = get_number(
        stats,
        "total length",
        integer=True,
        default=0,
    )

    bases_mapped = get_number(
        stats,
        "bases mapped",
        integer=True,
        default=0,
    )

    bases_mapped_cigar = get_number(
        stats,
        "bases mapped (cigar)",
        integer=True,
        default=bases_mapped,
    )

    mismatches = get_number(
        stats,
        "mismatches",
        integer=True,
        default=0,
    )

    error_rate = get_number(
        stats,
        "error rate",
        default=0.0,
    )

    average_length = get_number(
        stats,
        "average length",
        default=0.0,
    )

    maximum_length = get_number(
        stats,
        "maximum length",
        integer=True,
        default=0,
    )

    non_primary = get_number(
        stats,
        "non-primary alignments",
        integer=True,
        default=0,
    )

    # These fields are optional because some samtools versions
    # do not include them in the SN section.
    insertions = get_number(
        stats,
        "insertions",
        integer=True,
        default="",
    )

    deletions = get_number(
        stats,
        "deletions",
        integer=True,
        default="",
    )

    mapping_percent = (
        reads_mapped / raw_total * 100
        if raw_total > 0
        else 0.0
    )

    mapped_bases_percent = (
        bases_mapped_cigar / total_length * 100
        if total_length > 0
        else 0.0
    )

    technology_label = "ONT" if technology == "ont" else "PacBio"

    return {
        "sample": sample,
        "technology": technology_label,
        "stats_file": path.name,
        "raw_total_sequences": raw_total,
        "reads_mapped": reads_mapped,
        "reads_unmapped": reads_unmapped,
        "mapping_percent": f"{mapping_percent:.4f}",
        "total_length": total_length,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": bases_mapped_cigar,
        "mapped_bases_percent": f"{mapped_bases_percent:.4f}",
        "mismatches": mismatches,
        "error_rate": f"{error_rate:.8f}",
        "error_percent": f"{error_rate * 100:.4f}",
        "insertions": insertions,
        "deletions": deletions,
        "average_length": f"{average_length:.2f}",
        "maximum_length": maximum_length,
        "non_primary_alignments": non_primary,
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)

    selected = choose_files()
    missing = EXPECTED - set(selected)

    if missing:
        missing_text = ", ".join(
            f"{sample}-{technology}"
            for sample, technology in sorted(missing)
        )
        raise FileNotFoundError(
            f"Missing samtools stats files: {missing_text}"
        )

    rows = []

    for sample, technology in sorted(EXPECTED):
        path = selected[(sample, technology)]

        try:
            row = build_row(
                path,
                sample,
                technology,
            )
        except Exception as error:
            raise RuntimeError(
                f"Could not process {path.name}: {error}"
            ) from error

        rows.append(row)

    technology_order = {
        "ONT": 0,
        "PacBio": 1,
    }

    rows.sort(
        key=lambda row: (
            row["sample"],
            technology_order[row["technology"]],
        )
    )

    fieldnames = [
        "sample",
        "technology",
        "stats_file",
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "mapping_percent",
        "total_length",
        "bases_mapped",
        "bases_mapped_cigar",
        "mapped_bases_percent",
        "mismatches",
        "error_rate",
        "error_percent",
        "insertions",
        "deletions",
        "average_length",
        "maximum_length",
        "non_primary_alignments",
    ]

    with OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    print("SUCCESS: alignment summary created")
    print(OUTPUT)
    print()

    for row in rows:
        print(
            f"{row['sample']} "
            f"{row['technology']:<6} "
            f"mapped={row['mapping_percent']}% "
            f"error={row['error_percent']}%"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
