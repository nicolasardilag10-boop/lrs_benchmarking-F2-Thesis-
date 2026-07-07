#!/usr/bin/env python3

"""
Plot mapped bases across alignment configurations.

The figure contains one panel per sample:
    HG002
    HG003
    HG004

Within each panel:
    - one line represents ONT reads
    - one line represents PacBio reads
    - points show the four alignment configurations

The plotted value is:

    bases_mapped - minimum(bases_mapped across all results)
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter


# ============================================================
# 1. File paths
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

OUTPUT_PNG = (
    FIGURES_DIR
    / "mapped_bases_by_configuration.png"
)

OUTPUT_PDF = (
    FIGURES_DIR
    / "mapped_bases_by_configuration.pdf"
)


# ============================================================
# 2. Expected benchmark structure
# ============================================================

SAMPLES = [
    "HG002",
    "HG003",
    "HG004",
]

TECHNOLOGIES = [
    "ONT",
    "PacBio",
]

CONFIGURATIONS = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]


# False reproduces the screenshot more closely because each
# sample panel receives its own y-axis range.
#
# Change to True for direct comparison between all samples.
SHARE_Y_AXIS = False


# ============================================================
# 3. Read and validate data
# ============================================================

def read_data() -> pd.DataFrame:
    """Read the alignment summary table."""

    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Input table was not found: {INPUT_FILE}"
        )

    data = pd.read_csv(
        INPUT_FILE,
        sep="\t",
    )

    required_columns = {
        "sample",
        "read_technology",
        "configuration",
        "bases_mapped",
    }

    missing_columns = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    data["bases_mapped"] = pd.to_numeric(
        data["bases_mapped"],
        errors="raise",
    )

    return data


def validate_data(data: pd.DataFrame) -> None:
    """Check that the complete 24-condition benchmark exists."""

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, but found {len(data)}."
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
            "Duplicate benchmark conditions detected:\n"
            + duplicates.to_string(index=False)
        )

    expected_conditions = {
        (sample, technology, configuration)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for configuration in CONFIGURATIONS
    }

    observed_conditions = set(
        zip(
            data["sample"],
            data["read_technology"],
            data["configuration"],
        )
    )

    missing_conditions = (
        expected_conditions - observed_conditions
    )

    if missing_conditions:
        missing_text = "\n".join(
            (
                f"  {sample} | "
                f"{technology} | "
                f"{configuration}"
            )
            for sample, technology, configuration
            in sorted(missing_conditions)
        )

        raise ValueError(
            "Missing benchmark conditions:\n"
            + missing_text
        )


# ============================================================
# 4. Create the figure
# ============================================================

def create_plot(data: pd.DataFrame) -> None:
    """Create the three-panel mapped-bases figure."""

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    global_minimum = data["bases_mapped"].min()

    data = data.copy()

    data["mapped_bases_delta"] = (
        data["bases_mapped"] - global_minimum
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(14, 5.5),
        sharey=SHARE_Y_AXIS,
    )

    x_positions = list(
        range(len(CONFIGURATIONS))
    )

    # Obtain six consistent colors from a Matplotlib palette.
    color_map = plt.get_cmap("tab10")

    series_colors = {
        (sample, technology): color_map(index)
        for index, (sample, technology) in enumerate(
            (
                (sample, technology)
                for sample in SAMPLES
                for technology in TECHNOLOGIES
            )
        )
    }

    legend_handles = []

    for axis, sample in zip(
        axes,
        SAMPLES,
    ):
        sample_data = data.loc[
            data["sample"] == sample
        ].copy()

        for technology in TECHNOLOGIES:
            technology_data = (
                sample_data.loc[
                    sample_data["read_technology"]
                    == technology,
                    [
                        "configuration",
                        "mapped_bases_delta",
                    ],
                ]
                .set_index("configuration")
                .reindex(CONFIGURATIONS)
            )

            values = (
                technology_data["mapped_bases_delta"]
                .astype(float)
                .tolist()
            )

            series_name = (
                f"{sample}_"
                f"{technology.lower()}"
            )

            series_color = series_colors[
                (sample, technology)
            ]

            axis.plot(
                x_positions,
                values,
                marker="o",
                linestyle=":",
                linewidth=1.4,
                markersize=6,
                color=series_color,
                label=series_name,
            )

            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle=":",
                    linewidth=1.4,
                    markersize=6,
                    color=series_color,
                    label=series_name,
                )
            )

        axis.set_title(
            sample,
            fontsize=10,
            pad=7,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.55",
                "boxstyle": "square,pad=0.25",
            },
        )

        axis.set_xticks(
            x_positions
        )

        axis.set_xticklabels(
            CONFIGURATIONS,
            rotation=20,
            ha="right",
            fontsize=8,
        )

        axis.grid(
            axis="y",
            linestyle=":",
            alpha=0.35,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

        # Show large values using scientific notation.
        formatter = ScalarFormatter(
            useMathText=True
        )

        formatter.set_scientific(True)

        formatter.set_powerlimits(
            (0, 0)
        )

        axis.yaxis.set_major_formatter(
            formatter
        )

    figure.suptitle(
        "Mapped bases across alignment configurations",
        fontsize=15,
        y=0.98,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=11,
        y=0.08,
    )

    figure.supylabel(
        "Mapped bases above global minimum",
        fontsize=11,
        x=0.03,
    )

    # Remove duplicate handles caused by plotting one panel at a time.
    unique_handles = {}

    for handle in legend_handles:
        unique_handles[handle.get_label()] = handle

    figure.legend(
        unique_handles.values(),
        unique_handles.keys(),
        title="Input dataset",
        loc="center right",
        bbox_to_anchor=(1.01, 0.55),
        frameon=False,
        fontsize=8,
    )

    figure.text(
        0.05,
        0.015,
        (
            "Plotted value: bases_mapped − "
            f"global minimum ({global_minimum:,.0f} bases)"
        ),
        fontsize=9,
    )

    figure.tight_layout(
        rect=[
            0.04,
            0.10,
            0.88,
            0.93,
        ]
    )

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

    print("Mapped-bases figure created:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


# ============================================================
# 5. Main program
# ============================================================

def main() -> None:
    data = read_data()

    validate_data(
        data
    )

    global_minimum = data["bases_mapped"].min()

    print(
        "Global minimum mapped bases: "
        f"{global_minimum:,.0f}"
    )

    summary = (
        data.groupby(
            "read_technology"
        )["bases_mapped"]
        .mean()
        .sort_values(
            ascending=False
        )
    )

    print()
    print("Mean mapped bases by technology:")

    for technology, value in summary.items():
        print(
            f"{technology}: {value:,.0f}"
        )

    create_plot(
        data
    )


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)