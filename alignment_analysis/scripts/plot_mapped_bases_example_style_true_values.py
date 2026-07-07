#!/usr/bin/env python3

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
    / "mapped_bases_example_style_true_values.png"
)

OUTPUT_PDF = (
    FIGURES_DIR
    / "mapped_bases_example_style_true_values.pdf"
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
    "mm2-ont": "mm2-ont",
    "mm2-pb": "mm2-pb",
    "pbmm2-ont": "pbmm2-ont",
    "pbmm2-pb": "pbmm2-pb",
}


def load_data() -> pd.DataFrame:
    data = pd.read_csv(
        INPUT_FILE,
        sep="\t",
    )

    required = {
        "sample",
        "read_technology",
        "configuration",
        "bases_mapped",
    }

    missing = sorted(
        required - set(data.columns)
    )

    if missing:
        raise ValueError(
            "Missing columns: "
            + ", ".join(missing)
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
    keys = [
        "sample",
        "read_technology",
        "configuration",
    ]

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, found {len(data)}."
        )

    duplicates = data.duplicated(
        keys,
        keep=False,
    )

    if duplicates.any():
        raise ValueError(
            "Duplicated conditions detected."
        )

    expected = {
        (sample, technology, configuration)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for configuration in CONFIGURATIONS
    }

    observed = set(
        data[
            keys
        ].itertuples(
            index=False,
            name=None,
        )
    )

    if expected != observed:
        raise ValueError(
            "The 24 expected conditions do not match the table."
        )


def create_plot(data: pd.DataFrame) -> None:
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(14.5, 5.5),
        sharey=True,
    )

    x_positions = list(
        range(len(CONFIGURATIONS))
    )

    minimum = float(
        data["bases_mapped"].min()
    )

    maximum = float(
        data["bases_mapped"].max()
    )

    padding = (
        maximum - minimum
    ) * 0.08

    shared_minimum = max(
        0,
        minimum - padding,
    )

    shared_maximum = (
        maximum + padding
    )

    for axis, sample in zip(
        axes,
        SAMPLES,
    ):
        sample_data = data.loc[
            data["sample"] == sample
        ]

        for technology in TECHNOLOGIES:
            values = (
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
                ["bases_mapped"]
                .to_numpy(dtype=float)
            )

            axis.plot(
                x_positions,
                values,
                marker="o",
                markersize=6,
                linewidth=1.5,
                linestyle=":",
                label=technology,
            )

        axis.set_title(
            sample,
            fontsize=10,
            pad=8,
            bbox={
                "facecolor": "white",
                "edgecolor": "0.4",
                "boxstyle": "square,pad=0.35",
            },
        )

        axis.set_xticks(
            x_positions
        )

        axis.set_xticklabels(
            [
                DISPLAY_LABELS[item]
                for item in CONFIGURATIONS
            ],
            fontsize=8,
        )

        axis.set_ylim(
            shared_minimum,
            shared_maximum,
        )

        axis.yaxis.set_major_formatter(
            FuncFormatter(
                lambda value, position: (
                    f"{value / 1_000_000:.1f}"
                )
            )
        )

        axis.grid(False)

        axis.spines["top"].set_visible(True)
        axis.spines["right"].set_visible(True)

    figure.suptitle(
        "Absolute mapped bases across alignment configurations",
        fontsize=14,
        y=0.98,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=11,
        y=0.07,
    )

    figure.supylabel(
        "Mapped bases (millions)",
        fontsize=11,
        x=0.02,
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
        0.04,
        0.015,
        (
            "Plotted values are the original bases_mapped counts. "
            "All panels use the same y-axis scale."
        ),
        fontsize=8,
    )

    figure.tight_layout(
        rect=[
            0.04,
            0.12,
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

    print("Figure created:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


def main() -> None:
    data = load_data()
    validate_data(data)

    print(
        "Original bases_mapped range: "
        f"{data['bases_mapped'].min():,.0f} to "
        f"{data['bases_mapped'].max():,.0f}"
    )

    create_plot(data)


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)
