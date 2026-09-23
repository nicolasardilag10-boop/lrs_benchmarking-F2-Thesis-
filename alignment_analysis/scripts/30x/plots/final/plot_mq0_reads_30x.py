#!/usr/bin/env python3
"""Plot MAPQ 0 reads (%) for minimap2, pbmm2, VACmap and VG Giraffe.

Grouped bars per sample per aligner, from a zero baseline.
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

INPUT_TSV = PROJECT / "alignment_analysis" / "tables" / "30x" / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_PNG = PROJECT / "alignment_analysis" / "figures" / "30x" / "final" / "04_mq0_reads_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")
OUTPUT_SUMMARY = (
    PROJECT / "alignment_analysis" / "tables" / "30x" / "derived" / "plot_data" / "mq0_reads_percent.tsv"
)

BAR_WIDTH = 0.18
SAMPLE_OFFSETS = {"HG002": -0.21, "HG003": 0.0, "HG004": 0.21}


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {"sample", "read_technology", "aligner", "reads_mq0", "reads_mapped"}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    for column in ["reads_mq0", "reads_mapped"]:
        data[column] = pd.to_numeric(data[column], errors="raise")

    data["reads_mq0_percent"] = np.where(
        data["reads_mapped"] > 0,
        data["reads_mq0"] / data["reads_mapped"] * 100.0,
        np.nan,
    )

    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER),
        ["sample", "read_technology", "aligner", "reads_mq0", "reads_mapped", "reads_mq0_percent"],
    ].copy()

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError(
            "Duplicated sample/technology/aligner rows found:\n"
            + plot_data.loc[duplicated].to_string(index=False)
        )

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

    y_max = float(plot_data["reads_mq0_percent"].max()) * 1.18
    x_positions = np.arange(len(ALIGNER_ORDER))

    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(FULL_WIDTH_IN, 3.2))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for sample in SAMPLE_ORDER:
            sample_data = (
                technology_data[technology_data["sample"] == sample]
                .set_index("aligner")
                .reindex(ALIGNER_ORDER)
            )
            values = sample_data["reads_mq0_percent"].to_numpy(dtype=float)
            axis.bar(
                x_positions + SAMPLE_OFFSETS[sample],
                values,
                width=BAR_WIDTH,
                color=SAMPLE_COLORS[sample],
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=6)
        panel_letter(axis, "a" if technology == "ONT" else "b", fontsize=12)
        axis.set_xticks(x_positions, ALIGNER_ORDER, fontsize=12, rotation=45)
        axis.set_xlim(-0.6, len(ALIGNER_ORDER) - 0.4)
        axis.set_ylim(0, y_max)
        axis.tick_params(axis="y", labelsize=12)
        subtle_grid(axis, "y")
        clean_spines(axis)
        axis.set_facecolor("white")

    axes[0].set_ylabel("MAPQ 0 reads (%)", fontsize=13)

    figure.legend(
        handles=sample_legend_handles(),
        frameon=False,
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        columnspacing=1.3,
        handletextpad=0.4,
        fontsize=12,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.12, top=0.82, wspace=0.10)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
