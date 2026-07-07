#!/usr/bin/env python3

"""
Plot absolute mapped bases for all 24 benchmark conditions.

Structure:
- one panel per sample: HG002, HG003, HG004
- one line for ONT reads
- one line for PacBio reads
- four alignment configurations
- shared y-axis across all panels

This is an absolute mapped-yield figure. It should not be interpreted
as a normalized comparison of alignment efficiency because the ONT and
PacBio inputs may contain different total numbers of bases.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd


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

OUTPUT_PNG = (
    FIGURES_DIR
    / "absolute_mapped_bases_by_sample.png"
)

OUTPUT_PDF = (
    FIGURES_DIR
    / "absolute_mapped_bases_by_sample.pdf"
)


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

DISPLAY_LABELS = {
    "mm2-ont": "minimap2\nmap-ont",
    "mm2-pb": "minimap2\nmap-hifi",
    "pbmm2-ont": "pbmm2\nSUBREAD",
    "pbmm2-pb": "pbmm2\nCCS/HiFi",
}


def read_data() -> pd.DataFrame:
    """Read and standardize the alignment summary table."""

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

    return data


def validate_data(data: pd.DataFrame) -> None:
    """Verify that the complete 24-condition benchmark is present."""

    key_columns = [
        "sample",
        "read_technology",
        "configuration",
    ]

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, but found {len(data)}."
        )

    duplicates = data.duplicated(
        key_columns,
        keep=False,
    )

    if duplicates.any():
        raise ValueError(
            "Duplicated benchmark conditions detected:\n"
            + data.loc[
                duplicates,
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


def create_plot(data: pd.DataFrame) -> None:
    """Create the three-panel absolute mapped-bases figure."""

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(14.5, 5.6),
        sharey=True,
    )

    x_positions = list(
        range(len(CONFIGURATIONS))
    )

    global_minimum = float(
        data["bases_mapped"].min()
    )

    global_maximum = float(
        data["bases_mapped"].max()
    )

    y_range = global_maximum - global_minimum

    y_padding = (
        y_range * 0.08
        if y_range > 0
        else global_maximum * 0.05
    )

    shared_y_minimum = max(
        0,
        global_minimum - y_padding,
    )

    shared_y_maximum = (
        global_maximum + y_padding
    )

    for axis, sample in zip(
        axes,
        SAMPLES,
    ):
        sample_data = data.loc[
            data["sample"] == sample
        ]

        for technology in TECHNOLOGIES:
            technology_data = (
                sample_data.loc[
                    sample_data["read_technology"]
                    == technology,
                    [
                        "configuration",
                        "bases_mapped",
                    ],
                ]
                .set_index("configuration")
                .reindex(CONFIGURATIONS)
            )

            if technology_data["bases_mapped"].isna().any():
                raise ValueError(
                    f"Missing values after reindexing "
                    f"{sample} {technology}."
                )

            values = (
                technology_data["bases_mapped"]
                .astype(float)
                .tolist()
            )

            axis.plot(
                x_positions,
                values,
                marker="o",
                linestyle=":",
                linewidth=1.8,
                markersize=6,
                label=technology,
            )

        axis.set_title(
            sample,
            fontsize=12,
            pad=10,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.55",
                "boxstyle": "square,pad=0.3",
            },
        )

        axis.set_xticks(
            x_positions
        )

        axis.set_xticklabels(
            [
                DISPLAY_LABELS[configuration]
                for configuration in CONFIGURATIONS
            ],
            fontsize=8,
        )

        axis.tick_params(
            axis="both",
            direction="out",
        )

        axis.grid(
            axis="y",
            linestyle=":",
            alpha=0.30,
        )

        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

        axis.set_ylim(
            shared_y_minimum,
            shared_y_maximum,
        )

        axis.yaxis.set_major_formatter(
            FuncFormatter(
                lambda value, position: f"{value / 1_000_000:.1f}"
            )
        )

    figure.suptitle(
        "Absolute mapped bases across alignment configurations",
        fontsize=15,
        y=0.98,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=11,
        y=0.06,
    )

    figure.supylabel(
        "Mapped bases (millions)",
        fontsize=11,
        x=0.025,
    )

    handles, labels = (
        axes[0].get_legend_handles_labels()
    )

    figure.legend(
        handles,
        labels,
        title="Read technology",
        loc="center right",
        bbox_to_anchor=(0.99, 0.55),
        frameon=False,
    )

    figure.text(
        0.045,
        0.012,
        (
            "Absolute mapped-base counts are influenced by total input yield "
            "and should not be interpreted as normalized alignment efficiency."
        ),
        fontsize=8,
    )

    figure.tight_layout(
        rect=[
            0.04,
            0.11,
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

    print("Absolute mapped-bases figure created:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


def main() -> None:
    data = read_data()

    validate_data(
        data
    )

    print(
        "Mapped-bases range: "
        f"{data['bases_mapped'].min():,.0f} to "
        f"{data['bases_mapped'].max():,.0f}"
    )

    print()
    print("Mean mapped bases by read technology:")

    summary = (
        data.groupby(
            "read_technology"
        )["bases_mapped"]
        .mean()
        .reindex(TECHNOLOGIES)
    )

    for technology, value in summary.items():
        print(
            f"  {technology}: {value:,.0f}"
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
