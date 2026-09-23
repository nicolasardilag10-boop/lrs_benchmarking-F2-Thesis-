#!/usr/bin/env python3
"""Plot CIGAR-aligned, mismatch and CIGAR-unaligned bases (%) as a stacked bar.

Three-segment stack (aligned non-mismatch / mismatch / CIGAR-unaligned)
per sample per aligner, from a truncated 88-100% baseline.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch

from utils.plot_style import (
    ALIGNER_ORDER,
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    apply_style,
    panel_letter,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "derived" / "plot_data" / "cigar_mapped_mismatch_unaligned_percent.tsv"
OUTPUT_PNG = FIGURE_DIR / "06_cigar_composition_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.23, "HG003": 0.00, "HG004": 0.23}
BAR_WIDTH = 0.19


def lighten_color(color, amount=0.5):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {"sample", "read_technology", "aligner", "total_length", "bases_mapped_cigar", "mismatches"}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    for column in ["total_length", "bases_mapped_cigar", "mismatches"]:
        data[column] = pd.to_numeric(data[column], errors="raise")

    if (data["total_length"] <= 0).any():
        raise ValueError("At least one total_length value is <= 0.")

    data["cigar_mapped_percent"] = data["bases_mapped_cigar"] / data["total_length"] * 100.0
    data["cigar_unaligned_percent"] = 100.0 - data["cigar_mapped_percent"]
    data["mismatch_total_percent"] = data["mismatches"] / data["total_length"] * 100.0
    data["aligned_non_mismatch_percent"] = (
        (data["bases_mapped_cigar"] - data["mismatches"]) / data["total_length"] * 100.0
    )

    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER),
        [
            "sample", "read_technology", "aligner", "total_length", "bases_mapped_cigar",
            "mismatches", "cigar_mapped_percent", "cigar_unaligned_percent",
            "mismatch_total_percent", "aligned_non_mismatch_percent",
        ],
    ].copy()

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology/aligner rows found:\n" + plot_data.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, found {len(plot_data)}.")

    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    plot_data["aligner"] = pd.Categorical(plot_data["aligner"], ALIGNER_ORDER, ordered=True)
    return plot_data.sort_values(["read_technology", "aligner", "sample"])


def main() -> int:
    plot_data = load_data()

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    apply_style()

    figure, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(FULL_WIDTH_IN, 6.4))
    x_positions = np.arange(len(ALIGNER_ORDER))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for sample in SAMPLE_ORDER:
            sample_data = (
                technology_data[technology_data["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
            )
            aligned_values = sample_data["aligned_non_mismatch_percent"].to_numpy(dtype=float)
            mismatch_values = sample_data["mismatch_total_percent"].to_numpy(dtype=float)
            unaligned_values = sample_data["cigar_unaligned_percent"].to_numpy(dtype=float)
            cigar_mapped_values = sample_data["cigar_mapped_percent"].to_numpy(dtype=float)
            positions = x_positions + SAMPLE_OFFSETS[sample]

            base_color = SAMPLE_COLORS[sample]
            mismatch_color = lighten_color(base_color, amount=0.35)
            unaligned_color = lighten_color(base_color, amount=0.72)

            axis.bar(positions, aligned_values, width=BAR_WIDTH, color=base_color, edgecolor="white", linewidth=0.5, zorder=3)
            mismatch_bars = axis.bar(
                positions, mismatch_values, width=BAR_WIDTH, bottom=aligned_values,
                color=mismatch_color, edgecolor="white", linewidth=0.5, zorder=3,
            )
            axis.bar(
                positions, unaligned_values, width=BAR_WIDTH, bottom=aligned_values + mismatch_values,
                color=unaligned_color, edgecolor="white", linewidth=0.5, zorder=3,
            )

            for mismatch_bar, mismatch_value in zip(mismatch_bars, mismatch_values):
                if np.isnan(mismatch_value):
                    continue
                x_center = mismatch_bar.get_x() + mismatch_bar.get_width() / 2
                y_center = mismatch_bar.get_y() + mismatch_bar.get_height() / 2
                if mismatch_value >= 0.5:
                    axis.text(
                        x_center, y_center, f"{mismatch_value:.2f}",
                        ha="center", va="center", fontsize=7.5, color="black", zorder=7,
                    )

            for position, mapped_value in zip(positions, cigar_mapped_values):
                if np.isnan(mapped_value):
                    continue
                axis.text(
                    position, 100.3, f"{mapped_value:.2f}",
                    ha="center", va="bottom", fontsize=7.5, color=base_color,
                    rotation=45, clip_on=False, zorder=8,
                )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=12, y=1.10)
        panel_letter(axis, "a" if technology == "ONT" else "b", fontsize=12)
        axis.set_xticks(x_positions, ALIGNER_ORDER, fontsize=12, rotation=45)
        axis.set_xlim(-0.6, len(ALIGNER_ORDER) - 0.4)
        axis.set_ylim(88, 100)
        axis.set_yticks([88, 90, 92, 94, 96, 98, 100])
        axis.tick_params(axis="y", labelsize=12)
        axis.grid(axis="y", color="#E5E5E5", linewidth=0.45, zorder=0)
        axis.set_axisbelow(True)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    figure.supylabel("CIGAR-aligned and\nunaligned bases (%)", fontsize=12, x=-0.02, y=0.47)

    sample_legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="white", label=sample) for sample in SAMPLE_ORDER
    ]
    figure.legend(
        handles=sample_legend_handles, frameon=False, ncols=3, loc="upper center",
        bbox_to_anchor=(0.5, 1.00), columnspacing=1.3, handletextpad=0.4, fontsize=12,
    )

    component_legend_handles = [
        Patch(facecolor="#555555", edgecolor="white", label="Aligned, non-mismatch"),
        Patch(facecolor="#999999", edgecolor="white", label="Mismatch bases"),
        Patch(facecolor="#DDDDDD", edgecolor="white", label="CIGAR-unaligned bases"),
    ]
    figure.legend(
        handles=component_legend_handles, frameon=False, ncols=3, loc="upper center",
        bbox_to_anchor=(0.5, 0.975), columnspacing=1.3, handletextpad=0.4, fontsize=12,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.13, right=0.98, bottom=0.11, top=0.84, hspace=0.33)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
