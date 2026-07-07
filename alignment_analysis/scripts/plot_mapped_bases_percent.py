#!/usr/bin/env python3

from pathlib import Path
import sys

import matplotlib.pyplot as plt
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
    / "mapped_bases_percent_by_sample.png"
)

OUTPUT_PDF = (
    FIGURES_DIR
    / "mapped_bases_percent_by_sample.pdf"
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
    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Input table not found: {INPUT_FILE}"
        )

    data = pd.read_csv(
        INPUT_FILE,
        sep="\t",
    )

    required_columns = {
        "sample",
        "read_technology",
        "configuration",
        "mapped_bases_percent",
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

    data["mapped_bases_percent"] = pd.to_numeric(
        data["mapped_bases_percent"],
        errors="raise",
    )

    return data


def validate_data(data: pd.DataFrame) -> None:
    key_columns = [
        "sample",
        "read_technology",
        "configuration",
    ]

    if len(data) != 24:
        raise ValueError(
            f"Expected 24 rows, found {len(data)}."
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

    if (
        data["mapped_bases_percent"] > 100.0001
    ).any():
        raise ValueError(
            "mapped_bases_percent exceeds 100%."
        )

    if (
        data["mapped_bases_percent"] < 0
    ).any():
        raise ValueError(
            "Negative mapped-bases percentage detected."
        )

    expected = {
        (sample, technology, configuration)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for configuration in CONFIGURATIONS
    }

    observed = set(
        data[
            key_columns
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
            "Missing benchmark conditions:\n"
            + "\n".join(
                str(condition)
                for condition in missing
            )
        )

    if unexpected:
        raise ValueError(
            "Unexpected benchmark conditions:\n"
            + "\n".join(
                str(condition)
                for condition in unexpected
            )
        )


def create_plot(data: pd.DataFrame) -> None:
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(14, 5.5),
        sharey=True,
    )

    x_positions = list(
        range(len(CONFIGURATIONS))
    )

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
                        "mapped_bases_percent",
                    ],
                ]
                .set_index("configuration")
                .reindex(CONFIGURATIONS)
            )

            values = (
                technology_data[
                    "mapped_bases_percent"
                ]
                .astype(float)
                .tolist()
            )

            axis.plot(
                x_positions,
                values,
                marker="o",
                linewidth=1.8,
                markersize=6,
                label=technology,
            )

        axis.set_title(
            sample,
            fontsize=11,
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

        axis.set_ylim(
            0,
            102,
        )

        axis.grid(
            axis="y",
            linestyle=":",
            alpha=0.35,
        )

        axis.spines[
            "top"
        ].set_visible(False)

        axis.spines[
            "right"
        ].set_visible(False)

    figure.suptitle(
        "Mapped input bases across alignment configurations",
        fontsize=15,
        y=0.98,
    )

    figure.supxlabel(
        "Alignment configuration",
        fontsize=11,
        y=0.04,
    )

    figure.supylabel(
        "Mapped input bases (%)",
        fontsize=11,
        x=0.02,
    )

    handles, labels = axes[0].get_legend_handles_labels()

    figure.legend(
        handles,
        labels,
        title="Read technology",
        loc="center right",
        bbox_to_anchor=(0.99, 0.53),
        frameon=False,
    )

    figure.tight_layout(
        rect=[
            0.04,
            0.09,
            0.89,
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

    print("Mapped-bases percentage figure created:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)


def main() -> None:
    data = read_data()

    validate_data(
        data
    )

    print(
        "Percentage range: "
        f"{data['mapped_bases_percent'].min():.4f}%–"
        f"{data['mapped_bases_percent'].max():.4f}%"
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
