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

    alignment_analysis/tables/alignment_summary_30x_samtools.tsv
"""

import csv
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

TABLES = (
    PROJECT
    / "samtools_stats_30x_Christian"
)

VACMAP_VG_TABLES = (
    PROJECT
    / "samtools_stats_30x_Christian"
    / "vacmap_vg_stats"
)

OUTPUT = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv"
)


# ---------------------------------------------------------------------
# Supported filename structures
# ---------------------------------------------------------------------

# Input files are read from:

#samtools_stats_30x_Christian/

#The combined table is written to:

#    samtools_stats_30x_Christian/alignment_summary_30x_samtools.tsv

OLD_PATTERN = re.compile(
    r"^(HG00[234])_"
    r"(ont|pb)_"
    r"30x\.hg38\."
    r"(mm2-(?:ont|pb)|pbmm2-(?:ont|pb))"
    r"\.cram\.stats\.SN\.txt$"
)

NEW_PATTERN = re.compile(
    r"^(HG00[234])\."
    r"(ont|pb)\."
    r"30x\.hg38\."
    r"(vacmap-(?:ont|pb)|vg-(?:ont|pb))"
    r"\.cram\.stats$")


# ---------------------------------------------------------------------
# Expected results at the current stage of the project
# ---------------------------------------------------------------------

VALID_CONFIGURATIONS = (
    ("ont", "mm2-ont"),
    ("ont", "pbmm2-ont"),
    ("ont", "vacmap-ont"),
    ("ont", "vg-ont"),

    ("pb", "mm2-pb"),
    ("pb", "pbmm2-pb"),
    ("pb", "vacmap-pb"),
    ("pb", "vg-pb"),)

EXPECTED = {
    (sample, technology, configuration)
    for sample in ("HG002", "HG003", "HG004")
    for technology, configuration in VALID_CONFIGURATIONS}

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
    default=None,):
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
    method_tag: str,
) -> tuple[str, str, str]:
    """
    Return:

        aligner
        preset
        configuration
    """

    if method_tag == "mm2-ont":
        return "minimap2", "map-ont", "mm2-ont"

    if method_tag == "mm2-pb":
        return "minimap2", "map-hifi", "mm2-pb"

    if method_tag == "pbmm2-ont":
        return "pbmm2", "SUBREAD", "pbmm2-ont"

    if method_tag == "pbmm2-pb":
        return "pbmm2", "CCS", "pbmm2-pb"

    if method_tag == "vacmap-ont":
        return "VACmap", "vacmap-ont", "vacmap-ont"

    if method_tag == "vacmap-pb":
        return "VACmap", "vacmap-pb", "vacmap-pb"

    if method_tag == "vg-ont":
        return "VG Giraffe", "vg-ont", "vg-ont"

    if method_tag == "vg-pb":
        return "VG Giraffe", "vg-pb", "vg-pb"

    raise ValueError(
        f"Unsupported alignment configuration: {method_tag}"
    )

# ---------------------------------------------------------------------
# Select input files
# ---------------------------------------------------------------------

def choose_files():
    selected = {}

    # ---------------------------------------------------------
    # minimap2 / pbmm2
    # ---------------------------------------------------------
    for path in sorted(
        TABLES.glob("*.cram.stats.SN.txt")
    ):
        match = OLD_PATTERN.match(path.name)

        if match is None:
            continue

        sample = match.group(1)
        technology = match.group(2)
        method = match.group(3)

        dataset = (
            sample,
            technology,
            method,
        )

        if dataset in selected:
            raise ValueError(
                f"Duplicate statistics file for {dataset}"
            )

        selected[dataset] = path

    # ---------------------------------------------------------
    # VACmap / VG Giraffe
    # ---------------------------------------------------------
    for path in sorted(
        VACMAP_VG_TABLES.glob("*.cram.stats")
    ):
        match = NEW_PATTERN.match(path.name)

        if match is None:
            continue

        sample = match.group(1)
        technology = match.group(2)
        method = match.group(3)

        dataset = (
            sample,
            technology,
            method,
        )

        if dataset in selected:
            raise ValueError(
                f"Duplicate statistics file for {dataset}"
            )

        selected[dataset] = path

    return selected

# ---------------------------------------------------------------------
# Build one table row
# ---------------------------------------------------------------------

def build_row(
    path: Path,
    sample: str,
    technology: str,
    method: str,) -> dict[str, object]:
    """Create one combined summary row."""

    stats = parse_sn_file(path)

    raw_total = get_number(
        stats,
        "raw total sequences",
        integer=True,)

    reads_mapped = get_number(
        stats,
        "reads mapped",
        integer=True,)

    if raw_total is None:
        raise ValueError("missing 'raw total sequences'")

    if reads_mapped is None:
        raise ValueError("missing 'reads mapped'")

    reads_unmapped = get_number(
        stats,
        "reads unmapped",
        integer=True,
        default=raw_total - reads_mapped,)

    total_length = get_number(
        stats,
        "total length",
        integer=True,
        default=0,)

    bases_mapped = get_number(
        stats,
        "bases mapped",
        integer=True,
        default=0,)

    bases_mapped_cigar = get_number(
        stats,
        "bases mapped (cigar)",
        integer=True,
        default=bases_mapped,)

    mismatches = get_number(
        stats,
        "mismatches",
        integer=True,
        default=0,)

    error_rate = get_number(
        stats,
        "error rate",
        default=0.0,)

    insertions = get_number(
        stats,
        "insertions",
        integer=True,
        default=None,)

    deletions = get_number(
        stats,
        "deletions",
        integer=True,
        default=None,)

    average_length = get_number(
        stats,
        "average length",
        default=0.0,)

    maximum_length = get_number(
        stats,
        "maximum length",
        integer=True,
        default=0,)

    non_primary_alignments = get_number(
        stats,
        "non-primary alignments",
        integer=True,
        default=0,)

    mapped_reads_percent = (
        reads_mapped / raw_total * 100
        if raw_total > 0
        else 0.0)

    mapped_bases_percent = (
        bases_mapped / total_length * 100
        if total_length > 0
        else 0.0)

    error_percent = error_rate * 100

    aligner, preset, configuration = describe_alignment(technology, method,)

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
        "statistics_file": str(path.relative_to(PROJECT)),
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
        "non_primary_alignments": non_primary_alignments,}


# ---------------------------------------------------------------------
# Sort rows consistently
# ---------------------------------------------------------------------

def sorting_key(row: dict[str, object]) -> tuple[int, int, int]:
    """Define the presentation order in the output table."""

    sample_order = {
        "HG002": 0,
        "HG003": 1,
        "HG004": 2,}

    technology_order = {
        "ONT": 0,
        "PacBio": 1,}

    configuration_order = {
        "mm2-ont": 0,
        "pbmm2-ont": 1,
        "vacmap-ont": 2,
        "vg-ont": 3,

        "mm2-pb": 4,
        "pbmm2-pb": 5,
        "vacmap-pb": 6,
        "vg-pb": 7,}

    return (
        sample_order[str(row["sample"])],
        technology_order[str(row["read_technology"])],
        configuration_order[str(row["configuration"])],)


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main() -> int:
    """Build and write the combined TSV table."""

    if not TABLES.is_dir():
        print(f"ERROR: tables directory not found: {TABLES}", file=sys.stderr,)
        return 1

    selected = choose_files()

    missing = sorted(EXPECTED - set(selected))

    if missing:
        print("WARNING: the following expected statistics files are missing:", file=sys.stderr,)

        for sample, technology, method in missing:
            print(f"  {sample} | {technology} | {method}", file=sys.stderr,)

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
                method,)
        except ValueError as error:
            print(
                f"ERROR: {path.name}: {error}",
                file=sys.stderr,)
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
        "non_primary_alignments",]

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,)


    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="raise",)

        writer.writeheader()
        writer.writerows(rows)

    print(f"Alignment summary written to: {OUTPUT}")
    print(f"Data rows written: {len(rows)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
