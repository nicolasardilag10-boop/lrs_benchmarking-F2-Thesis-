#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# Paths
# ============================================================

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

TABLE_DIR = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
)

base_file = (
    TABLE_DIR
    / "source/alignment_summary_30x_combined_metrics_without_ram.tsv"
)

mm2_file = (
    TABLE_DIR
    / "source/mm2.run_metrics.tsv"
)

output_file = (
    TABLE_DIR
    / "final/alignment_benchmark_30x.tsv"
)


# ============================================================
# Load base table
# ============================================================

base = pd.read_csv(
    base_file,
    sep="\t",
)

print("\nBase table:")
print(base_file)
print("Rows:", len(base))


# ============================================================
# Load minimap2 run metrics
#
# IMPORTANT:
# source/mm2.run_metrics.tsv has NO header.
# ============================================================

mm2 = pd.read_csv(
    mm2_file,
    sep="\t",
    header=None,
    names=[
        "sample",
        "technology",
        "coverage",
        "configuration",
        "metrics",
    ],
)

print("\nMinimap2 run metrics:")
print(mm2_file)
print("Rows:", len(mm2))


# ============================================================
# Parse runtime
# ============================================================

mm2["runtime_seconds"] = pd.to_numeric(
    mm2["metrics"].str.extract(
        r"Real time:\s*([0-9.]+)\s*sec",
        expand=False,
    ),
    errors="coerce",
)

mm2["runtime_hours"] = (
    mm2["runtime_seconds"] / 3600.0
)


def seconds_to_hms(seconds):

    if pd.isna(seconds):
        return pd.NA

    total = int(round(seconds))

    hours = total // 3600
    minutes = (total % 3600) // 60
    secs = total % 60

    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


mm2["runtime_hms"] = (
    mm2["runtime_seconds"]
    .apply(seconds_to_hms)
)


# ============================================================
# Normalize source identifiers
# ============================================================

mm2["sample"] = (
    mm2["sample"]
    .astype(str)
    .str.strip()
)

mm2["technology"] = (
    mm2["technology"]
    .astype(str)
    .str.strip()
    .str.lower()
)

mm2["configuration"] = (
    mm2["configuration"]
    .astype(str)
    .str.strip()
    .str.lower()
)

mm2["coverage"] = (
    mm2["coverage"]
    .astype(str)
    .str.strip()
    .str.lower()
)


# ============================================================
# Keep 30x only
# ============================================================

mm2 = mm2[
    mm2["coverage"] == "30x"
].copy()


# ============================================================
# IMPORTANT:
#
# Keep only the technology-matched minimap2 configuration:
#
# ONT     -> mm2-ont
# PacBio  -> mm2-pb
#
# Cross-preset experiments are NOT used.
# ============================================================

matched_mm2 = mm2[
    (
        (mm2["technology"] == "ont")
        & (mm2["configuration"] == "mm2-ont")
    )
    |
    (
        (mm2["technology"] == "pb")
        & (mm2["configuration"] == "mm2-pb")
    )
].copy()


# Match technology naming used by alignment summary
matched_mm2["read_technology"] = (
    matched_mm2["technology"]
    .map({
        "ont": "ONT",
        "pb": "PacBio",
    })
)


# ============================================================
# Show values that will actually be used
# ============================================================

print("\n========================================")
print("MINIMAP2 RUNTIME VALUES TO BE USED")
print("========================================")

print(
    matched_mm2[
        [
            "sample",
            "read_technology",
            "configuration",
            "runtime_seconds",
            "runtime_hours",
            "runtime_hms",
        ]
    ].to_string(index=False)
)


# ============================================================
# Check duplicate sample/technology combinations
# ============================================================

duplicates = matched_mm2.duplicated(
    [
        "sample",
        "read_technology",
    ],
    keep=False,
)

if duplicates.any():

    print("\nDuplicated matched minimap2 rows:")

    print(
        matched_mm2.loc[
            duplicates,
            [
                "sample",
                "read_technology",
                "configuration",
            ]
        ].to_string(index=False)
    )

    raise ValueError(
        "Duplicate minimap2 sample × technology combinations."
    )


# ============================================================
# Build lookup table
# ============================================================

lookup = matched_mm2[
    [
        "sample",
        "read_technology",
        "configuration",
        "runtime_hms",
        "runtime_seconds",
        "runtime_hours",
    ]
].copy()

lookup = lookup.rename(
    columns={
        "configuration": "mm2_runtime_configuration",
        "runtime_hms": "mm2_runtime_hms",
        "runtime_seconds": "mm2_runtime_seconds",
        "runtime_hours": "mm2_runtime_hours",
    }
)


# ============================================================
# Normalize base table identifiers
# ============================================================

base["sample"] = (
    base["sample"]
    .astype(str)
    .str.strip()
)

base["read_technology"] = (
    base["read_technology"]
    .astype(str)
    .str.strip()
    .replace({
        "PacBio HiFi": "PacBio",
        "PB": "PacBio",
        "pb": "PacBio",
        "ont": "ONT",
    })
)

base["aligner"] = (
    base["aligner"]
    .replace({
        "VACMap": "VACmap",
    })
)


# ============================================================
# Merge runtime lookup
# ============================================================

merged = base.merge(
    lookup,
    on=[
        "sample",
        "read_technology",
    ],
    how="left",
    validate="many_to_one",
)


# ============================================================
# Fill ONLY minimap2 rows
#
# Existing VACmap/VG values remain unchanged.
# Existing non-empty values are also preserved.
# ============================================================

minimap2_mask = (
    merged["aligner"]
    .astype(str)
    .str.lower()
    .eq("minimap2")
)


def fill_missing(target_column, source_column):

    mask = (
        minimap2_mask
        & merged[target_column].isna()
        & merged[source_column].notna()
    )

    merged.loc[
        mask,
        target_column,
    ] = merged.loc[
        mask,
        source_column,
    ]


fill_missing(
    "wallclock_runtime_hms",
    "mm2_runtime_hms",
)

fill_missing(
    "wallclock_runtime_seconds",
    "mm2_runtime_seconds",
)

fill_missing(
    "wallclock_runtime_hours",
    "mm2_runtime_hours",
)


# ============================================================
# Fill provenance
# ============================================================

source_mask = (
    minimap2_mask
    & merged["wallclock_runtime_seconds"].notna()
    & merged["wallclock_runtime_source_log"].isna()
)

merged.loc[
    source_mask,
    "wallclock_runtime_source_log",
] = "source/mm2.run_metrics.tsv"


# ============================================================
# Optional configuration provenance
# ============================================================

if "runtime_configuration" not in merged.columns:
    merged["runtime_configuration"] = pd.NA

configuration_mask = (
    minimap2_mask
    & merged["mm2_runtime_configuration"].notna()
)

merged.loc[
    configuration_mask,
    "runtime_configuration",
] = merged.loc[
    configuration_mask,
    "mm2_runtime_configuration",
]


# ============================================================
# Remove temporary merge columns
# ============================================================

merged = merged.drop(
    columns=[
        "mm2_runtime_configuration",
        "mm2_runtime_hms",
        "mm2_runtime_seconds",
        "mm2_runtime_hours",
    ]
)


# ============================================================
# Save NEW table
#
# Do NOT overwrite original until verified.
# ============================================================

merged.to_csv(
    output_file,
    sep="\t",
    index=False,
)


# ============================================================
# Verification
# ============================================================

check = merged[
    merged["aligner"]
    .astype(str)
    .str.lower()
    .eq("minimap2")
][
    [
        "sample",
        "read_technology",
        "aligner",
        "runtime_configuration",
        "wallclock_runtime_hms",
        "wallclock_runtime_seconds",
        "wallclock_runtime_hours",
        "wallclock_runtime_source_log",
    ]
]


print("\n========================================")
print("MINIMAP2 ROWS AFTER MERGE")
print("========================================")

print(
    check.to_string(index=False)
)


print("\n========================================")
print("COMPLETENESS")
print("========================================")

print(
    "Total minimap2 rows:",
    len(check)
)

print(
    "Runtime values filled:",
    check["wallclock_runtime_seconds"]
    .notna()
    .sum()
)

print(
    "Still missing:",
    check["wallclock_runtime_seconds"]
    .isna()
    .sum()
)


print("\nOutput written to:")
print(output_file)
