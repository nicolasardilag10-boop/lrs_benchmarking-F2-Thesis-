#!/usr/bin/env python3

"""
Create a paired error-rate plot comparing ONT and PacBio.

Input:
    alignment_analysis/tables/alignment_summary.tsv

Output:
    alignment_analysis/figures/alignment_error_rate.png
    alignment_analysis/figures/alignment_error_rate.pdf
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


# ============================================================
# 1. Define project paths
# ============================================================

PROJECT_DIR = Path.home() / "lrs_benchmarking"

INPUT_FILE = (
    PROJECT_DIR
    / "alignment_analysis"
    / "tables"
    / "alignment_summary.tsv"
)

FIGURES_DIR = (
    PROJECT_DIR
    / "alignment_analysis"
    / "figures"
)

OUTPUT_PNG = FIGURES_DIR / "alignment_error_rate.png"
OUTPUT_PDF = FIGURES_DIR / "alignment_error_rate.pdf"


# ============================================================
# 2. Define expected samples and technologies
# ============================================================

EXPECTED_SAMPLES = [
    "HG002",
    "HG003",
    "HG004",
]

EXPECTED_TECHNOLOGIES = {
    "ONT",
    "PacBio",
}


# ============================================================
# 3. Read the alignment summary table
# ============================================================

def read_alignment_summary(
    input_file: Path,
) -> dict[str, dict[str, float]]:
    """
    Read error percentages from the alignment summary.

    The returned structure looks like:

        {
            "HG002": {
                "ONT": 5.2,
                "PacBio": 0.4
            }
        }
    """

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input table not found: {input_file}"
        )

    results: dict[str, dict[str, float]] = {}

    with input_file.open(
        mode="r",
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required_columns = {
            "sample",
            "technology",
            "error_percent",
        }

        missing_columns = (
            required_columns
            - set(reader.fieldnames or [])
        )

        if missing_columns:
            raise ValueError(
                "Missing required columns: "
                + ", ".join(sorted(missing_columns))
            )

        for row in reader:

            sample = row["sample"].strip()
            technology = row["technology"].strip()

            if sample not in EXPECTED_SAMPLES:
                continue

            if technology not in EXPECTED_TECHNOLOGIES:
                raise ValueError(
                    f"Unexpected technology: {technology}"
                )

            try:
                error_percent = float(
                    row["error_percent"]
                )
            except ValueError as error:
                raise ValueError(
                    f"Invalid error percentage for "
                    f"{sample}-{technology}: "
                    f"{row['error_percent']}"
                ) from error

            results.setdefault(sample, {})

            if technology in results[sample]:
                raise ValueError(
                    f"Duplicate value for "
                    f"{sample}-{technology}"
                )

            results[sample][technology] = error_percent

    return results


# ============================================================
# 4. Validate paired observations
# ============================================================

def validate_results(
    results: dict[str, dict[str, float]],
) -> None:
    """
    Confirm that every sample contains both technologies.
    """

    for sample in EXPECTED_SAMPLES:

        if sample not in results:
            raise ValueError(
                f"Missing sample: {sample}"
            )

        missing_technologies = (
            EXPECTED_TECHNOLOGIES
            - set(results[sample])
        )

        if missing_technologies:
            raise ValueError(
                f"{sample} is missing: "
                + ", ".join(
                    sorted(missing_technologies)
                )
            )

        for technology, value in results[sample].items():

            if not 0 <= value <= 100:
                raise ValueError(
                    f"Invalid error percentage for "
                    f"{sample}-{technology}: {value}"
                )


# ============================================================
# 5. Create the paired error-rate figure
# ============================================================

def create_plot(
    results: dict[str, dict[str, float]],
) -> None:
    """
    Draw one paired line per sample.

    x = sequencing technology
    y = alignment error percentage
    """

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    technologies = [
        "ONT",
        "PacBio",
    ]

    x_positions = [
        0,
        1,
    ]

    figure, axis = plt.subplots(
        figsize=(7, 5)
    )

    for sample in EXPECTED_SAMPLES:

        values = [
            results[sample]["ONT"],
            results[sample]["PacBio"],
        ]

        # Draw one line connecting the same sample
        # across the two sequencing technologies.
        axis.plot(
            x_positions,
            values,
            marker="o",
            linewidth=1.8,
            markersize=7,
            label=sample,
        )

        # Add the numeric percentage above each point.
        for x_position, value in zip(
            x_positions,
            values,
        ):
            axis.annotate(
                f"{value:.2f}%",
                xy=(x_position, value),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                fontsize=8,
            )

    axis.set_xticks(
        x_positions,
        technologies,
    )

    axis.set_ylabel(
        "Alignment error rate (%)"
    )

    axis.set_xlabel(
        "Sequencing technology"
    )

    axis.set_title(
        "Alignment error rate: ONT versus PacBio"
    )

    # Error percentages cannot be negative.
    axis.set_ylim(
        bottom=0
    )

    # Keep only horizontal reference lines.
    axis.grid(
        axis="y",
        linestyle=":",
        alpha=0.5,
    )

    axis.legend(
        title="Sample",
        frameon=False,
    )

    figure.tight_layout()

    figure.savefig(
        OUTPUT_PNG,
        dpi=300,
        bbox_inches="tight",
    )

    figure.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )

    plt.close(figure)


# ============================================================
# 6. Main program
# ============================================================

def main() -> None:
    results = read_alignment_summary(
        INPUT_FILE
    )

    validate_results(
        results
    )

    print("Error-rate values used:")

    for sample in EXPECTED_SAMPLES:
        print(
            f"{sample}: "
            f"ONT={results[sample]['ONT']:.4f}% | "
            f"PacBio={results[sample]['PacBio']:.4f}%"
        )

    create_plot(
        results
    )

    print()
    print("SUCCESS: error-rate figures created")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


if __name__ == "__main__":

    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)
