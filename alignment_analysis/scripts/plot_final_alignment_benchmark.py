#!/usr/bin/env python3

"""
Create the final alignment error-rate benchmark figure.

Panels:
    1. ONT input reads
    2. PacBio input reads

Each panel contains:
    - mm2-ont
    - mm2-pb
    - pbmm2-ont
    - pbmm2-pb

HG002, HG003, and HG004 are displayed as individual points
connected across configurations.
"""

from pathlib import Path
import statistics
import sys

import matplotlib.pyplot as plt
import pandas as pd


PROJECT = Path.home() / "lrs_benchmarking"

INPUT_TABLE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "alignment_summary.tsv"
)

FIGURE_DIR = (
    PROJECT
    / "alignment_analysis"
    / "figures"
)

PNG_OUTPUT = FIGURE_DIR / "final_alignment_error_rate.png"
PDF_OUTPUT = FIGURE_DIR / "final_alignment_error_rate.pdf"

CONFIGURATIONS = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]

TECHNOLOGIES = [
    "ONT",
    "PacBio",
]

SAMPLES = [
    "HG002",
    "HG003",
    "HG004",
]

SHOW_MEDIAN_LABELS = True


def validate_data(data: pd.DataFrame) -> None:
    """Confirm that the complete 24-condition matrix is present."""

    required_columns = {
        "sample",
        "read_technology",
        "configuration",
        "error_percent",
    }

    missing_columns = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            "Missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, but found {len(data)}"
        )

    duplicate_mask = data.duplicated(
        subset=[
            "sample",
            "read_technology",
            "configuration",
        ],
        keep=False,
    )

    if duplicate_mask.any():
        duplicates = data.loc[
            duplicate_mask,
            [
                "sample",
                "read_technology",
                "configuration",
            ],
        ]

        raise ValueError(
            "Duplicate conditions detected:\n"
            + duplicates.to_string(index=False)
        )

    expected = {
        (sample, technology, configuration)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for configuration in CONFIGURATIONS
    }

    observed = set(
        zip(
            data["sample"],
            data["read_technology"],
            data["configuration"],
        )
    )

    missing_conditions = expected - observed

    if missing_conditions:
        message = "\n".join(
            f"  {sample} | {technology} | {configuration}"
            for sample, technology, configuration
            in sorted(missing_conditions)
        )

        raise ValueError(
            "Missing benchmark conditions:\n" + message
        )


def main() -> int:
    """Read the summary table and generate PNG and PDF figures."""

    if not INPUT_TABLE.is_file():
        print(
            f"ERROR: summary table not found: {INPUT_TABLE}",
            file=sys.stderr,
        )
        return 1

    data = pd.read_csv(
        INPUT_TABLE,
        sep="\t",
    )

    try:
        validate_data(data)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    data["error_percent"] = pd.to_numeric(
        data["error_percent"],
        errors="raise",
    )

    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(14, 7),
        sharey=True,
    )

    x_positions = list(
        range(1, len(CONFIGURATIONS) + 1)
    )

    for axis, technology in zip(
        axes,
        TECHNOLOGIES,
    ):
        panel = data.loc[
            data["read_technology"] == technology
        ].copy()

        boxplot_data = []

        for configuration in CONFIGURATIONS:
            values = (
                panel.loc[
                    panel["configuration"] == configuration,
                    "error_percent",
                ]
                .astype(float)
                .tolist()
            )

            boxplot_data.append(values)

        axis.boxplot(
            boxplot_data,
            positions=x_positions,
            widths=0.55,
            showfliers=False,
            medianprops={'color':'black',"linewidth": 2}
            )

        if SHOW_MEDIAN_LABELS:
            medians = [
                statistics.median(values) if values else float("nan")
                for values in boxplot_data
            ]

            y_min, y_max = axis.get_ylim()
            offset = 0.04 * (y_max - y_min)

            for position, median_value in zip(x_positions, medians):
                axis.text(
                    position,
                    median_value + offset,
                    f"{median_value:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color="black",
                )

        for sample in SAMPLES:
            sample_data = (
                panel.loc[
                    panel["sample"] == sample,
                    [
                        "configuration",
                        "error_percent",
                    ],
                ]
                .set_index("configuration")
                .reindex(CONFIGURATIONS)
            )

            y_values = (
                sample_data["error_percent"]
                .astype(float)
                .tolist()
            )

            axis.plot(
                x_positions,
                y_values,
                marker="o",
                linestyle=":",
                linewidth=1.3,
                markersize=7,
                label=sample,
                zorder=3,
            )

        axis.set_title(
            f"{technology} input reads",
            fontsize=14,
        )


        axis.set_xticks(x_positions)

        axis.set_xticklabels(
            CONFIGURATIONS,
            rotation=25,
            ha="right",
        )

        axis.grid(
            axis="y",
            linestyle=":",
            alpha=False,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[0].set_ylabel(
        "Alignment error rate (%)"
    )

    handles, labels = axes[0].get_legend_handles_labels()

    figure.supxlabel(
    "Alignment",
    fontsize=14,
    y=0.1,
    )

    figure.legend(
        handles,
        labels,
        title="Sample",
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.01),
    )

    figure.suptitle(
        "Alignment error-rate benchmark",
        fontsize=16,
        y=1.07,
    )

    figure.tight_layout(
        rect=[0, 0.06, 1, 0.94]
    
    )

    figure.savefig(
        PNG_OUTPUT,
        dpi=300,
        bbox_inches="tight",
    )

    figure.savefig(
        PDF_OUTPUT,
        bbox_inches="tight",
    )

    plt.close(figure)

    print("Figure successfully created:")
    print(PNG_OUTPUT)
    print(PDF_OUTPUT)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
