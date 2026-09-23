#!/usr/bin/env python3
"""Plot alignment error rate (%) for minimap2, pbmm2, VACmap and VG Giraffe.

Grouped bar chart (one bar per HG002/HG003/HG004 sample within each
aligner) from a true zero baseline.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
except ModuleNotFoundError as error:
    raise SystemExit(
        f"A required plotting package is missing: {error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error

from utils.plot_style import (
    ALIGNER_ORDER,
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    apply_style,
    clean_spines,
    panel_letter,
    sample_legend_handles,
    save_figure,
    subtle_grid,
)

INPUT_FILE = (
    PROJECT / "alignment_analysis" / "tables" / "30x" / "final"
    / "alignment_benchmark_30x.tsv"
)
OUTPUT_PNG = (
    PROJECT / "alignment_analysis" / "figures" / "30x" / "final"
    / "01_alignment_error_rate_30x.png"
)
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

METRIC = "error_percent"

CONFIGURATION_TO_ALIGNER = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
    "vacmap-ont": "VACmap",
    "vacmap-pb": "VACmap",
    "vg-ont": "VG Giraffe",
    "vg-pb": "VG Giraffe",
}

SAMPLE_OFFSETS = {"HG002": -0.23, "HG003": 0.0, "HG004": 0.23}
BAR_WIDTH = 0.19


def load_plot_data() -> pd.DataFrame:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"The input table was not found:\n{INPUT_FILE}")

    data = pd.read_csv(INPUT_FILE, sep="\t")
    required_columns = {"sample", "read_technology", "configuration", METRIC}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"The input table is missing these columns: {sorted(missing_columns)}")

    data[METRIC] = pd.to_numeric(data[METRIC], errors="raise")

    selected = (
        ((data["read_technology"] == "ONT") & data["configuration"].isin(
            ["mm2-ont", "pbmm2-ont", "vacmap-ont", "vg-ont"]))
        | ((data["read_technology"] == "PacBio") & data["configuration"].isin(
            ["mm2-pb", "pbmm2-pb", "vacmap-pb", "vg-pb"]))
    )

    plot_data = data.loc[selected, ["sample", "read_technology", "configuration", METRIC]].copy()
    plot_data["aligner"] = plot_data["configuration"].map(CONFIGURATION_TO_ALIGNER)

    if plot_data["aligner"].isna().any():
        raise ValueError("At least one selected configuration could not be mapped to an aligner.")

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError(
            "Duplicated sample/technology/aligner observations found:\n"
            + plot_data.loc[duplicated].to_string(index=False)
        )

    expected = pd.MultiIndex.from_product(
        [TECHNOLOGY_ORDER, ALIGNER_ORDER, SAMPLE_ORDER],
        names=["read_technology", "aligner", "sample"],
    )
    observed = pd.MultiIndex.from_frame(plot_data[["read_technology", "aligner", "sample"]])
    missing = expected.difference(observed)
    if len(missing):
        raise ValueError(f"Missing sample/aligner/technology combinations: {list(missing)}")

    return plot_data


def main() -> int:
    plot_data = load_plot_data()
    apply_style()

    maximum_error = plot_data[METRIC].max()
    y_max = maximum_error * 1.15

    figure, axes = plt.subplots(
        nrows=1, ncols=2, sharey=True, figsize=(FULL_WIDTH_IN, 3.3))
    aligner_positions = list(range(len(ALIGNER_ORDER)))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for aligner_index, aligner in enumerate(ALIGNER_ORDER):
            values = (
                technology_data[technology_data["aligner"] == aligner]
                .set_index("sample")[METRIC]
                .reindex(SAMPLE_ORDER))

            for sample in SAMPLE_ORDER:
                value = values.loc[sample]
                x_position = aligner_index + SAMPLE_OFFSETS[sample]
                bars = axis.bar(
                    x_position,
                    value,
                    width=BAR_WIDTH,
                    color=SAMPLE_COLORS[sample],
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,)
                axis.annotate(
                    f"{value:.2f}",
                    xy=(bars[0].get_x() + bars[0].get_width() / 2, value),
                    xytext=(0, 2),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                    fontsize=7.5,
                    rotation=45,
                    zorder=4,)

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=6)
        panel = "a" if technology == "ONT" else "b"
        axis.text(
            -0.05,           # X position
            1.02,            # Y position
            panel,
            transform=axis.transAxes,
            fontsize=12,
            fontweight="bold",
            ha="left",
            va="bottom",
            )
        axis.set_ylim(0, y_max)
        axis.set_xlim(-0.5, len(ALIGNER_ORDER) - 0)
        axis.set_xticks(aligner_positions, ALIGNER_ORDER, fontsize=12, rotation=45)
        axis.tick_params(axis="y", labelsize=12)
        subtle_grid(axis, "y")
        clean_spines(axis)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Alignment error rate (%)", fontsize=12)

    figure.legend(
        handles=sample_legend_handles(),
        frameon=False,
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.06),
        columnspacing=1.3,
        handletextpad=0.4,
        fontsize=12,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.12, top=0.85, wspace=0.10)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_FILE}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
