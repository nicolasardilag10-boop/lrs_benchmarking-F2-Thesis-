#!/usr/bin/env python3
"""Plot input-normalized CIGAR-aligned yield vs mismatch error.

Promoted from scripts/30x/plots/exploratory/plot_cigar_yield_vs_error.py
(only its scatter panel; that script also drew an intermediate stacked
CIGAR-composition figure to the same output filename, which was always
overwritten by the scatter save immediately after, so it was dead output
and is not reproduced here). No regression line is drawn: the reported
statistic is the Spearman rank correlation (residuals are non-normal per
the Shapiro-Wilk diagnostics in Supplementary Table S1), and an OLS line
would visually imply a linear/Pearson relationship this figure is not
making. Pearson r is reported in Supplementary Table S1 rather than
repeated on this panel.
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
from scipy import stats

from utils.plot_style import (
    ALIGNER_MARKERS,
    ALIGNER_ORDER,
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    aligner_marker_handles,
    apply_style,
    clean_spines,
    panel_letter,
    sample_marker_handles,
    save_figure,
    subtle_grid,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_SUMMARY = TABLE_DIR / "derived" / "plot_data" / "cigar_yield_vs_error_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "13_cigar_yield_vs_error_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {"sample", "read_technology", "aligner", "total_length", "bases_mapped_cigar", "error_percent"}
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    for column in ["total_length", "bases_mapped_cigar", "error_percent"]:
        data[column] = pd.to_numeric(data[column], errors="raise")
    if (data["total_length"] <= 0).any():
        raise ValueError("At least one total_length value is <= 0.")

    selected = data[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER)
    ].copy()

    # Common per-sample/technology input denominator (largest available
    # total_length), matching the exploratory script this is promoted
    # from: minimap2/pbmm2/VG report the same total_length; VACmap's is
    # smaller, so using the max recovers the fullest input representation
    # until the true FASTQ base count is wired in.
    common_input = (
        selected.groupby(["sample", "read_technology"], observed=True)["total_length"]
        .transform("max")
    )
    selected["input_normalized_cigar_yield_percent"] = (
        selected["bases_mapped_cigar"] / common_input * 100.0
    )

    duplicated = selected.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology/aligner rows found:\n" + selected.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER)
    if len(selected) != expected_rows:
        raise ValueError(f"Expected {expected_rows} rows, found {len(selected)}.")

    columns = ["sample", "read_technology", "aligner", "input_normalized_cigar_yield_percent", "error_percent"]
    selected["sample"] = pd.Categorical(selected["sample"], SAMPLE_ORDER, ordered=True)
    selected["read_technology"] = pd.Categorical(selected["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    selected["aligner"] = pd.Categorical(selected["aligner"], ALIGNER_ORDER, ordered=True)
    return selected[columns].sort_values(["read_technology", "aligner", "sample"])


def main() -> int:
    plot_data = load_data()

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    apply_style()

    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(FULL_WIDTH_IN, 3.9))

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for aligner in ALIGNER_ORDER:
            aligner_data = technology_data[technology_data["aligner"] == aligner]
            for _, row in aligner_data.iterrows():
                axis.scatter(
                    row["input_normalized_cigar_yield_percent"], row["error_percent"],
                    marker=ALIGNER_MARKERS[aligner], color=SAMPLE_COLORS[str(row["sample"])],
                    edgecolor="black", linewidth=0.4, s=20, zorder=3,
                )

        spearman = stats.spearmanr(
            technology_data["input_normalized_cigar_yield_percent"], technology_data["error_percent"]
        )
        axis.text(
            0.04, 0.96,
            f"Spearman $\\rho$ = {spearman.statistic:.2f}, $P$ = {spearman.pvalue:.3f}\n$n$ = {len(technology_data)}",
            transform=axis.transAxes, ha="left", va="top", fontsize=10,
        )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=6)
        panel_letter(axis, "a" if technology == "ONT" else "b", fontsize=12)
        axis.tick_params(axis="both", labelsize=12)
        subtle_grid(axis, "both")
        clean_spines(axis)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Mismatch error (%)", fontsize=12)
    figure.supxlabel("Input-normalized CIGAR-aligned yield (%)", y=0.05, fontsize=12)

    # Side-by-side legend blocks: sample legend on the left, aligner legend
    # on the right, each laid out as its own horizontal row so neither
    # sits over the other and neither stacks its entries vertically.
    legend_samples = figure.legend(
        handles=sample_marker_handles(),
        title="GIAB sample",
        fontsize=9,
        title_fontsize=9,
        ncols=3,
        loc="upper left",
        bbox_to_anchor=(0.04, 1.0),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )
    figure.add_artist(legend_samples)

    figure.legend(
        handles=aligner_marker_handles(),
        title="Aligner",
        fontsize=9,
        title_fontsize=9,
        ncols=4,
        bbox_to_anchor=(0.99, 1.0),
        frameon=False,
        columnspacing=1.2,
        handletextpad=0.4,
    )

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.18, top=0.78, wspace=0.10)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
