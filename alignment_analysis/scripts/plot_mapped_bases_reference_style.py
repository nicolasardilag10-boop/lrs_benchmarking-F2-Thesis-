#!/usr/bin/env python3

"""
Create a reference-style mapped-bases comparison.

Plotted value:
    bases_mapped - minimum(bases_mapped across all 24 conditions)

Two figures are generated:

1. free_y:
   Closest to the supplied example. Each sample panel has its own
   y-axis limits.

2. shared_y:
   Scientifically safer for direct comparison because all panels
   use the same y-axis limits.

This script uses bases_mapped, not bases_mapped_cigar.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


# ============================================================
# Project paths
# ============================================================

PROJECT = Path.home() / "lrs_benchmarking"

INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "alignment_summary.tsv"
)

FIGURES_DIR = (
    PROJECT
    / "alignment_analysis"
    / "figures"
)


# ============================================================
# Benchmark structure
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

# Compact labels reproduce the reference-image structure.
X_LABELS = {
    "mm2-ont": "mm2-ont",
    "mm2-pb": "mm2-pb",
    "pbmm2-ont": "pbmm2-ont",
    "pbmm2-pb": "pbmm2-pb",
}

TECHNOLOGY_SUFFIX = {
    "ONT": "ont",
    "PacBio": "pb",
}

SERIES_ORDER = [
    (sample, technology)
    for sample in SAMPLES
    for technology in TECHNOLOGIES
]

palette = plt.get_cmap("tab10")

SERIES_COLORS = {
    series: palette(index)
    for index, series in enumerate(SERIES_ORDER)
}


# ============================================================
# Load and validate the table
# ============================================================

def load_data() -> tuple[pd.DataFrame, float]:
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Input table does not exist: {INPUT_FILE}"
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

    missing_columns = sorted(
        required_columns - set(data.columns)
    )

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing_columns)
        )

    data["sample"] = (
        data["sample"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    data["read_technology"] = (
        data["read_technology"]
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({
            "ont": "ONT",
            "pb": "PacBio",
            "pacbio": "PacBio",
        })
    )

    data["configuration"] = (
        data["configuration"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    data["bases_mapped"] = pd.to_numeric(
        data["bases_mapped"],
        errors="raise",
    )

    key_columns = [
        "sample",
        "read_technology",
        "configuration",
    ]

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, found {len(data)}."
        )

    duplicate_mask = data.duplicated(
        key_columns,
        keep=False,
    )

    if duplicate_mask.any():
        raise ValueError(
            "Duplicated benchmark conditions detected:\n"
            + data.loc[
                duplicate_mask,
                key_columns,
            ].to_string(index=False)
        )

    expected_conditions = {
        (sample, technology, configuration)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for configuration in CONFIGURATIONS
    }

    observed_conditions = set(
        data[
            key_columns
        ].itertuples(
            index=False,
            name=None,
        )
    )

    missing_conditions = sorted(
        expected_conditions - observed_conditions
    )

    unexpected_conditions = sorted(
        observed_conditions - expected_conditions
    )

    if missing_conditions:
        raise ValueError(
            "Missing benchmark conditions:\n"
            + "\n".join(
                str(condition)
                for condition in missing_conditions
            )
        )

    if unexpected_conditions:
        raise ValueError(
            "Unexpected benchmark conditions:\n"
            + "\n".join(
                str(condition)
                for condition in unexpected_conditions
            )
        )

    if data["bases_mapped"].isna().any():
        raise ValueError(
            "Missing bases_mapped values detected."
        )

    if (data["bases_mapped"] < 0).any():
        raise ValueError(
            "Negative bases_mapped values detected."
        )

    global_minimum = float(
        data["bases_mapped"].min()
    )

    data["mapped_bases_delta"] = (
        data["bases_mapped"]
        - global_minimum
    )

    return data, global_minimum


# ============================================================
# Figure generation
# ============================================================

def create_figure(
    data: pd.DataFrame,
    global_minimum: float,
    *,
    share_y_axis: bool,
    output_name: str,
) -> None:
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_png = (
        FIGURES_DIR
        / f"{output_name}.png"
    )

    output_pdf = (
        FIGURES_DIR
        / f"{output_name}.pdf"
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(14.5, 5.5),
        sharey=share_y_axis,
    )

    x_positions = np.arange(
        len(CONFIGURATIONS)
    )

    shared_maximum = float(
        data["mapped_bases_delta"].max()
    )

    for axis, sample in zip(
        axes,
        SAMPLES,
    ):
        sample_data = data.loc[
            data["sample"] == sample
        ]

        for technology in TECHNOLOGIES:
            series_data = (
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

            if series_data[
                "mapped_bases_delta"
            ].isna().any():
                raise ValueError(
                    "Missing plotted values after ordering: "
                    f"{sample}, {technology}"
                )

            values = (
                series_data["mapped_bases_delta"]
                .to_numpy(dtype=float)
            )

            axis.plot(
                x_positions,
                values,
                marker="o",
                markersize=6,
                linewidth=1.4,
                linestyle=(0, (1.5, 2.5)),
                color=SERIES_COLORS[
                    (sample, technology)
                ],
                label=(
                    f"{sample}_"
                    f"{TECHNOLOGY_SUFFIX[technology]}"
                ),
            )

        # Facet-style sample strip similar to the example.
        axis.set_title(
            sample,
            fontsize=10,
            pad=8,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.35",
                "boxstyle": "square,pad=0.35",
            },
        )

        axis.set_xticks(
            x_positions
        )

        axis.set_xticklabels(
            [
                X_LABELS[configuration]
                for configuration in CONFIGURATIONS
            ],
            fontsize=8,
        )

        axis.yaxis.set_major_formatter(
            FuncFormatter(
                lambda value, position: (
                    f"{value / 1_000_000:.1f}"
                )
            )
        )

        if share_y_axis:
            axis.set_ylim(
                0,
                shared_maximum * 1.08,
            )
        else:
            panel_values = sample_data[
                "mapped_bases_delta"
            ]

            panel_minimum = float(
                panel_values.min()
            )

            panel_maximum = float(
                panel_values.max()
            )

            panel_range = (
                panel_maximum - panel_minimum
            )

            padding = (
                panel_range * 0.10
                if panel_range > 0
                else max(
                    panel_maximum * 0.05,
                    1,
                )
            )

            axis.set_ylim(
                max(
                    0,
                    panel_minimum - padding,
                ),
                panel_maximum + padding,
            )

        axis.grid(False)

        axis.spines["top"].set_visible(True)
        axis.spines["right"].set_visible(True)

        axis.tick_params(
            axis="both",
            direction="out",
        )

    figure.suptitle(
        "Mapped-base differences across alignment configurations",
        fontsize=14,
        y=0.98,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=11,
        y=0.075,
    )

    figure.supylabel(
        "Mapped bases above global minimum (millions)",
        fontsize=11,
        x=0.02,
    )

    legend_handles = []

    for sample, technology in SERIES_ORDER:
        series_name = (
            f"{sample}_"
            f"{TECHNOLOGY_SUFFIX[technology]}"
        )

        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                markersize=6,
                linewidth=1.4,
                linestyle=(0, (1.5, 2.5)),
                color=SERIES_COLORS[
                    (sample, technology)
                ],
                label=series_name,
            )
        )

    figure.legend(
        handles=legend_handles,
        title="Input dataset",
        loc="center right",
        bbox_to_anchor=(0.995, 0.55),
        frameon=False,
        fontsize=8,
    )

    axis_description = (
        "Shared y-axis across panels."
        if share_y_axis
        else (
            "Panels use independent y-axis ranges "
            "to reproduce the reference style."
        )
    )

    figure.text(
        0.04,
        0.014,
        (
            "Plotted value: bases_mapped − "
            f"global minimum ({global_minimum:,.0f} bases). "
            f"{axis_description}"
        ),
        fontsize=8,
    )

    figure.tight_layout(
        rect=[
            0.04,
            0.12,
            0.86,
            0.93,
        ]
    )

    figure.savefig(
        output_png,
        dpi=300,
        bbox_inches="tight",
    )

    figure.savefig(
        output_pdf,
        bbox_inches="tight",
    )

    plt.close(figure)

    print(f"[CREATED] {output_png}")
    print(f"[CREATED] {output_pdf}")


# ============================================================
# Main
# ============================================================

def main() -> None:
    data, global_minimum = load_data()

    print(
        "Absolute mapped-bases range: "
        f"{data['bases_mapped'].min():,.0f} to "
        f"{data['bases_mapped'].max():,.0f}"
    )

    print(
        "Delta range: "
        f"{data['mapped_bases_delta'].min():,.0f} to "
        f"{data['mapped_bases_delta'].max():,.0f}"
    )

    create_figure(
        data,
        global_minimum,
        share_y_axis=False,
        output_name=(
            "mapped_bases_reference_style_free_y"
        ),
    )

    create_figure(
        data,
        global_minimum,
        share_y_axis=True,
        output_name=(
            "mapped_bases_reference_style_shared_y"
        ),
    )

    print()
    print("STATUS: BOTH REFERENCE-STYLE FIGURES CREATED")


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)
