#!/usr/bin/env python3

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np
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
    / "mapped_bases_reference_true_values.png"
)

OUTPUT_PDF = (
    FIGURES_DIR
    / "mapped_bases_reference_true_values.pdf"
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

LEGEND_LABELS = {
    ("HG002", "ONT"): "HG002_ont",
    ("HG002", "PacBio"): "HG002_pb",
    ("HG003", "ONT"): "HG003_ont",
    ("HG003", "PacBio"): "HG003_pb",
    ("HG004", "ONT"): "HG004_ont",
    ("HG004", "PacBio"): "HG004_pb",
}

COLOR_MAP = {
    ("HG002", "ONT"): "#F8766D",
    ("HG002", "PacBio"): "#B79F00",
    ("HG003", "ONT"): "#00BA38",
    ("HG003", "PacBio"): "#00BFC4",
    ("HG004", "ONT"): "#619CFF",
    ("HG004", "PacBio"): "#F564E3",
}


def load_data() -> pd.DataFrame:
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Input table not found: {INPUT_FILE}"
        )

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
            "Missing required columns: "
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
            "Duplicated benchmark conditions detected:\n"
            + data.loc[
                duplicates,
                keys,
            ].to_string(index=False)
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

    missing = sorted(
        expected - observed
    )

    unexpected = sorted(
        observed - expected
    )

    if missing:
        raise ValueError(
            "Missing conditions:\n"
            + "\n".join(
                str(condition)
                for condition in missing
            )
        )

    if unexpected:
        raise ValueError(
            "Unexpected conditions:\n"
            + "\n".join(
                str(condition)
                for condition in unexpected
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
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(16, 5.8),
        sharey=False,
    )

    x_positions = np.arange(
        len(CONFIGURATIONS)
    )

    for axis, sample in zip(
        axes,
        SAMPLES,
    ):
        sample_data = data.loc[
            data["sample"] == sample
        ]

        all_panel_values = []

        for technology in TECHNOLOGIES:
            series = (
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

            if series["bases_mapped"].isna().any():
                raise ValueError(
                    f"Missing values for {sample} {technology}."
                )

            values = (
                series["bases_mapped"]
                .to_numpy(dtype=float)
            )

            all_panel_values.extend(
                values.tolist()
            )

            axis.plot(
                x_positions,
                values,
                marker="o",
                markersize=6,
                linewidth=1.4,
                linestyle=(0, (1.5, 2.5)),
                color=COLOR_MAP[
                    (sample, technology)
                ],
            )

        panel_minimum = min(
            all_panel_values
        )

        panel_maximum = max(
            all_panel_values
        )

        panel_range = (
            panel_maximum - panel_minimum
        )

        padding = (
            panel_range * 0.12
            if panel_range > 0
            else panel_maximum * 0.05
        )

        axis.set_ylim(
            panel_minimum - padding,
            panel_maximum + padding,
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

        axis.yaxis.set_major_formatter(
            FuncFormatter(
                lambda value, position: (
                    f"{value / 1_000_000:.1f}"
                )
            )
        )

        axis.tick_params(
            axis="both",
            direction="out",
        )

        axis.grid(False)

        axis.spines["top"].set_visible(True)
        axis.spines["right"].set_visible(True)
        axis.spines["bottom"].set_visible(True)
        axis.spines["left"].set_visible(True)

        strip = Rectangle(
            (0, 1.0),
            1,
            0.08,
            transform=axis.transAxes,
            clip_on=False,
            facecolor="white",
            edgecolor="black",
            linewidth=1,
        )

        axis.add_patch(
            strip
        )

        axis.text(
            0.5,
            1.04,
            sample,
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )

    figure.suptitle(
        "Absolute mapped bases across alignment configurations",
        fontsize=15,
        y=0.99,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=12,
        y=0.06,
    )

    figure.supylabel(
        "Mapped bases (millions)",
        fontsize=12,
        x=0.015,
    )

    legend_handles = []

    for sample in SAMPLES:
        for technology in TECHNOLOGIES:
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    markersize=6,
                    linewidth=1.4,
                    linestyle=(0, (1.5, 2.5)),
                    color=COLOR_MAP[
                        (sample, technology)
                    ],
                    label=LEGEND_LABELS[
                        (sample, technology)
                    ],
                )
            )

    figure.legend(
        handles=legend_handles,
        title="Input dataset",
        loc="center right",
        bbox_to_anchor=(0.995, 0.55),
        frameon=False,
        fontsize=9,
    )

    figure.text(
        0.04,
        0.012,
        (
            "Original bases_mapped values are plotted without subtraction. "
            "Each sample panel uses its own y-axis range."
        ),
        fontsize=8,
    )

    figure.tight_layout(
        rect=[
            0.035,
            0.11,
            0.87,
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

    print()
    print("Figure created:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


def print_plotted_values(data: pd.DataFrame) -> None:
    table = (
        data.pivot_table(
            index=[
                "sample",
                "read_technology",
            ],
            columns="configuration",
            values="bases_mapped",
            aggfunc="first",
        )
        .reindex(
            columns=CONFIGURATIONS
        )
    )

    print()
    print("Values plotted:")
    print(table.to_string())

    print()
    print(
        "Overall range: "
        f"{data['bases_mapped'].min():,.0f} to "
        f"{data['bases_mapped'].max():,.0f}"
    )


def main() -> None:
    data = load_data()

    validate_data(
        data
    )

    print_plotted_values(
        data
    )

    create_plot(
        data
    )

    print()
    print("STATUS: TRUE-VALUE REFERENCE FIGURE CREATED")


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)
