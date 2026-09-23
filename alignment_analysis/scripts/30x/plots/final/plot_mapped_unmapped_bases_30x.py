#!/usr/bin/env python3
"""Plot mapped and unmapped bases (%) for minimap2, pbmm2, VACmap and VG Giraffe.

Stacked bar (mapped bases in the sample color, unmapped bases in a
lightened tint of the same color), from a truncated 90-100% baseline.

Denominator
-----------
Percentages are computed against the verified raw input base count per
sample/technology, not each aligner's own `total_length`. minimap2, pbmm2
and VG Giraffe keep unmapped reads in their output and agree exactly on
`total_length`. VACmap drops unmapped reads (`reads_unmapped == 0`), so its
own `total_length` excludes them and would report 100% mapped. Using the
shared raw input adds VACmap's dropped reads back as unmapped bases, making
it comparable to the other aligners. Same logic as plot_raw_total_count.py.
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
OUTPUT_SUMMARY = TABLE_DIR / "derived" / "plot_data" / "mapped_unmapped_bases_percent_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "12_mapped_unmapped_bases_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.22, "HG003": 0.00, "HG004": 0.22}
BAR_WIDTH = 0.19


def lighten_color(color, amount=0.55):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {"sample", "read_technology", "aligner", "total_length", "bases_mapped", "reads_unmapped"}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER),
        ["sample", "read_technology", "aligner", "total_length", "bases_mapped", "reads_unmapped"],
    ].copy()

    for column in ["total_length", "bases_mapped", "reads_unmapped"]:
        plot_data[column] = pd.to_numeric(plot_data[column], errors="raise")

    # Verified raw input per sample/technology: the total_length that all
    # aligners retaining unmapped reads (reads_unmapped > 0) agree on.
    retaining = plot_data.loc[plot_data["reads_unmapped"] > 0]
    raw_input = retaining.groupby(["sample", "read_technology"])["total_length"].agg(["nunique", "first"])
    if (raw_input["nunique"] != 1).any():
        raise ValueError(
            "Aligners retaining unmapped reads disagree on total_length:\n"
            + raw_input.loc[raw_input["nunique"] != 1].to_string()
        )
    plot_data = plot_data.merge(
        raw_input["first"].rename("raw_input_bases").reset_index(), on=["sample", "read_technology"], how="left",
    )
    if plot_data["raw_input_bases"].isna().any():
        raise ValueError("No aligner retains unmapped reads for at least one sample/technology.")
    if (plot_data["bases_mapped"] > plot_data["raw_input_bases"]).any():
        raise ValueError("bases_mapped exceeds the verified raw input for at least one row.")

    # Bases in reads the aligner dropped from its output (non-zero only for VACmap).
    plot_data["dropped_unmapped_bases"] = plot_data["raw_input_bases"] - plot_data["total_length"]
    plot_data["mapped_bases_percent"] = 100.0 * plot_data["bases_mapped"] / plot_data["raw_input_bases"]
    plot_data["unmapped_bases_percent"] = 100.0 - plot_data["mapped_bases_percent"]

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated observations found:\n" + plot_data.loc[duplicated].to_string(index=False))

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

    x_positions = np.arange(len(ALIGNER_ORDER))
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(FULL_WIDTH_IN, 4.1))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for sample in SAMPLE_ORDER:
            sample_data = (
                technology_data[technology_data["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
            )
            mapped_values = sample_data["mapped_bases_percent"].to_numpy(dtype=float)
            unmapped_values = sample_data["unmapped_bases_percent"].to_numpy(dtype=float)
            positions = x_positions + SAMPLE_OFFSETS[sample]
            base_color = SAMPLE_COLORS[sample]
            light_color = lighten_color(base_color, amount=0.55)

            axis.bar(positions, mapped_values, width=BAR_WIDTH, color=base_color, edgecolor="white", linewidth=0.6, zorder=3)
            axis.bar(
                positions, unmapped_values, width=BAR_WIDTH, bottom=mapped_values,
                color=light_color, edgecolor="white", linewidth=0.6, zorder=3,
            )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=12)
        panel_letter(axis, "a" if technology == "ONT" else "b", x=-0.05, y=1.02, fontsize=12)
        axis.set_xticks(x_positions, ALIGNER_ORDER, fontsize=12, rotation=45)
        axis.set_xlim(-0.6, len(ALIGNER_ORDER) - 0.4)
        axis.set_ylim(90, 100.35)
        axis.set_yticks([90, 92, 94, 96, 98, 100])
        axis.tick_params(axis="y", labelsize=12)
        axis.grid(axis="y", color="#E5E5E5", linewidth=0.45, zorder=0)
        axis.set_axisbelow(True)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Mapped and unmapped bases (%)", fontsize=13)

    sample_legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="white", label=sample) for sample in SAMPLE_ORDER
    ]
    figure.legend(
        handles=sample_legend_handles, frameon=False, ncols=3, loc="upper center",
        bbox_to_anchor=(0.5, 1.03), columnspacing=1.3, handletextpad=0.4, fontsize=12,
    )

    status_legend_handles = [
        Patch(facecolor="#666666", edgecolor="white", label="Mapped bases"),
        Patch(facecolor="#CCCCCC", edgecolor="white", label="Unmapped bases"),
    ]
    figure.legend(
        handles=status_legend_handles, frameon=False, ncols=2, loc="upper center",
        bbox_to_anchor=(0.5, 0.965), columnspacing=1.3, handletextpad=0.4, fontsize=12.0,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.10, top=0.74, wspace=0.10)

    figure.text(
        0.5, -0.12,
        "Note: VACmap removes unmapped reads from its output. Unmapped reads were added back using the\n"
        "raw input counts, to normalize VACmap with the other aligners.",
        ha="center", va="top", fontsize=9, color="#444444", style="italic",
    )

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
