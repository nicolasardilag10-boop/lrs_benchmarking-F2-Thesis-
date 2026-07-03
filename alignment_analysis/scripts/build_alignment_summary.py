#!/usr/bin/env python3

"""
Build one combined alignment summary table from samtools stats reports.

Currently supported alignment configurations:

1. ONT reads aligned with minimap2 map-ont
2. PacBio reads aligned with minimap2 map-hifi
3. PacBio reads aligned with pbmm2 CCS

Input files are read from:

    alignment_analysis/tables/

The combined table is written to:

    alignment_analysis/tables/alignment_summary.tsv
"""

import csv
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

PROJECT = Path.home() / "lrs_benchmarking"
TABLES = PROJECT / "alignment_analysis" / "tables"
OUTPUT = TABLES / "alignment_summary.tsv"


# ---------------------------------------------------------------------
# Supported filename structures
# ---------------------------------------------------------------------

# Examples accepted:
#
# HG002.ont.1k.samtools_stats.txt
# HG002.pb.1k.samtools_stats.txt
# HG002.pb.1k.pbmm2-ccs.samtools_stats.txt
#
# Older ".test" files are also accepted.
PATTERN = re.compile(
    r"^(HG00[234])\."
    r"(ont|pb)\."
    r"(1k|test)"
    r"(?:\.(pbmm2-(?:ccs|subread)))?"
    r"\.samtools_stats\.txt$"
)


# ---------------------------------------------------------------------
# Expected results at the current stage of the project
# ---------------------------------------------------------------------

EXPECTED = {
    ("HG002", "ont", "minimap2"),
    ("HG002", "pb", "minimap2"),
    ("HG002", "ont", "pbmm2-subread"),
    ("HG002", "pb", "pbmm2-subread"),
    ("HG002", "pb", "pbmm2-ccs"),

    ("HG003", "ont", "minimap2"),
    ("HG003", "pb", "minimap2"),
    ("HG003", "ont", "pbmm2-subread"),
    ("HG003", "pb", "pbmm2-subread"),
    ("HG003", "pb", "pbmm2-ccs"),

    ("HG004", "ont", "minimap2"),
    ("HG004", "pb", "minimap2"),
    ("HG004", "ont", "pbmm2-subread"),
    ("HG004", "pb", "pbmm2-subread"),
    ("HG004", "pb", "pbmm2-ccs"),
}


# ---------------------------------------------------------------------
# Read the SN section of a samtools stats file
# ---------------------------------------------------------------------

def parse_sn_file(path: Path) -> dict[str, str]:
    """Read the SN summary fields from one samtools stats report."""

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


# ---------------------------------------------------------------------
# Convert one samtools value to a number
# ---------------------------------------------------------------------

def get_number(
    statistics: dict[str, str],
    name: str,
    *,
    integer: bool = False,
    default=None,
):
    """
    Return one numeric SN value.

    Samtools values sometimes contain explanatory text after the number.
    Only the first field is converted.
    """

    raw_value = statistics.get(name)

    if raw_value is None:
        return default

    number_text = raw_value.split()[0].replace(",", "")

    try:
        number = float(number_text)
    except ValueError:
        return default

    return int(number) if integer else number


# ---------------------------------------------------------------------
# Determine the aligner metadata from the filename
# ---------------------------------------------------------------------

def describe_alignment(
    technology: str,
    method_tag: str | None,
) -> tuple[str, str, str]:
    """
    Return:

        aligner
        preset
        configuration

    The configuration column contains the compact labels intended for
    the final comparison figures.
    """

    if method_tag == "pbmm2-ccs":
        return "pbmm2", "CCS/HIFI", "pbmm2-pb"

    if method_tag == "pbmm2-subread":
        return "pbmm2", "SUBREAD", "pbmm2-ont"

    if technology == "ont":
        return "minimap2", "map-ont", "mm2-ont"

    if technology == "pb":
        return "minimap2", "map-hifi", "mm2-pb"

    raise ValueError(f"Unsupported technology: {technology}")


# ---------------------------------------------------------------------
# Select input files
# ---------------------------------------------------------------------

def choose_files() -> dict[tuple[str, str, str], Path]:
    """
    Find all supported samtools stats reports.

    The selection key is:

        sample, technology, alignment method

    When both ".1k" and ".test" versions exist for the same condition,
    the standardized ".1k" file is preferred.
    """

    selected: dict[tuple[str, str, str], Path] = {}
    selected_suffix: dict[tuple[str, str, str], str] = {}

    for path in sorted(TABLES.glob("*.samtools_stats.txt")):
        match = PATTERN.match(path.name)

        if match is None:
            continue

        sample = match.group(1)
        technology = match.group(2)
        suffix = match.group(3)
        method_tag = match.group(4)

        method = method_tag if method_tag else "minimap2"
        dataset = (sample, technology, method)

        if dataset not in selected:
            selected[dataset] = path
            selected_suffix[dataset] = suffix
            continue

        # Prefer ".1k" over the older ".test" filename.
        if suffix == "1k" and selected_suffix[dataset] != "1k":
            selected[dataset] = path
            selected_suffix[dataset] = suffix

    return selected


# ---------------------------------------------------------------------
# Build one table row
# ---------------------------------------------------------------------

def build_row(
    path: Path,
    sample: str,
    technology: str,
    method: str,
) -> dict[str, object]:
    """Create one combined summary row."""

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

    insertions = get_number(
        stats,
        "insertions",
        integer=True,
        default=None,
    )

    deletions = get_number(
        stats,
        "deletions",
        integer=True,
        default=None,
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

    non_primary_alignments = get_number(
        stats,
        "non-primary alignments",
        integer=True,
        default=0,
    )

    mapped_reads_percent = (
        reads_mapped / raw_total * 100
        if raw_total > 0
        else 0.0
    )

    mapped_bases_percent = (
        bases_mapped_cigar / total_length * 100
        if total_length > 0
        else 0.0
    )

    error_percent = error_rate * 100

    method_tag = None if method == "minimap2" else method

    aligner, preset, configuration = describe_alignment(
        technology,
        method_tag,
    )

    technology_label = {
        "ont": "ONT",
        "pb": "PacBio",
    }[technology]

    return {
        "sample": sample,
        "read_technology": technology_label,
        "aligner": aligner,
        "preset": preset,
        "configuration": configuration,
        "statistics_file": path.name,
        "raw_total_sequences": raw_total,
        "reads_mapped": reads_mapped,
        "reads_unmapped": reads_unmapped,
        "mapped_reads_percent": f"{mapped_reads_percent:.4f}",
        "total_length": total_length,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": bases_mapped_cigar,
        "mapped_bases_percent": f"{mapped_bases_percent:.4f}",
        "mismatches": mismatches,
        "error_rate": f"{error_rate:.8f}",
        "error_percent": f"{error_percent:.4f}",
        "insertions": insertions,
        "deletions": deletions,
        "average_length": f"{average_length:.2f}",
        "maximum_length": maximum_length,
        "non_primary_alignments": non_primary_alignments,
    }


# ---------------------------------------------------------------------
# Sort rows consistently
# ---------------------------------------------------------------------

def sorting_key(row: dict[str, object]) -> tuple[int, int, int]:
    """Define the presentation order in the output table."""

    sample_order = {
        "HG002": 0,
        "HG003": 1,
        "HG004": 2,
    }

    technology_order = {
        "ONT": 0,
        "PacBio": 1,
    }

    configuration_order = {
        "mm2-ont": 0,
        "mm2-pb": 1,
        "pbmm2-ont": 2,
        "pbmm2-pb": 3,
    }

    return (
        sample_order[str(row["sample"])],
        technology_order[str(row["read_technology"])],
        configuration_order[str(row["configuration"])],
    )


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main() -> int:
    """Build and write the combined TSV table."""

    if not TABLES.is_dir():
        print(
            f"ERROR: tables directory not found: {TABLES}",
            file=sys.stderr,
        )
        return 1

    selected = choose_files()

    missing = sorted(EXPECTED - set(selected))

    if missing:
        print(
            "ERROR: the following expected statistics files are missing:",
            file=sys.stderr,
        )

        for sample, technology, method in missing:
            print(
                f"  {sample} | {technology} | {method}",
                file=sys.stderr,
            )

        return 1

    rows: list[dict[str, object]] = []

    for dataset, path in selected.items():
        if dataset not in EXPECTED:
            continue

        sample, technology, method = dataset

        try:
            row = build_row(
                path,
                sample,
                technology,
                method,
            )
        except ValueError as error:
            print(
                f"ERROR: {path.name}: {error}",
                file=sys.stderr,
            )
            return 1

        rows.append(row)

    rows.sort(key=sorting_key)

    fieldnames = [
        "sample",
        "read_technology",
        "aligner",
        "preset",
        "configuration",
        "statistics_file",
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "mapped_reads_percent",
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

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="raise",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Alignment summary written to: {OUTPUT}")
    print(f"Data rows written: {len(rows)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
