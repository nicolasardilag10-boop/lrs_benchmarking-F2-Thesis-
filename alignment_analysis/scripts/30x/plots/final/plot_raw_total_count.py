#!/usr/bin/env python3
"""Plot raw input bases as CIGAR-aligned + CIGAR-unaligned bases (Gb).

Absolute stacked bar (CIGAR-aligned bases in the sample color, CIGAR-
unaligned bases in a lightened tint of the same color). Full stack height
equals the verified raw input base count for that sample/technology.

Raw input verification
-----------------------
`total_length` (from samtools stats) only equals the true FASTQ input base
count for an aligner whose output retains unmapped reads. minimap2, pbmm2
and VG Giraffe all report identical `total_length`/`raw_total_sequences`
per sample+technology here (verified below), because all three keep
unmapped reads in their CRAM output. VACmap reports `reads_unmapped == 0`
in every single row: it drops unmapped reads from its own output, so its
own `total_length` undercounts the true input. The verified per-sample/
technology raw input base count is therefore the value the three
retaining aligners agree on, applied to all four aligners -- never a
max/mean/median proxy across aligners.
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
OUTPUT_SUMMARY = TABLE_DIR / "derived" / "plot_data" / "raw_base_counts_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "10_raw_base_counts_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

SAMPLE_OFFSETS = {"HG002": -0.22, "HG003": 0.00, "HG004": 0.22}
BAR_WIDTH = 0.19

ROUNDING_TOLERANCE_BASES = 1.0


def lighten_color(color, amount=0.55):
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


def darken_color(color, amount=0.35):
    rgb = np.array(mcolors.to_rgb(color))
    black = np.array([0.0, 0.0, 0.0])
    return tuple(rgb + (black - rgb) * amount)


def load_and_verify_raw_input() -> pd.DataFrame:
    if not INPUT_TSV.exists():
        raise FileNotFoundError(f"Input TSV not found:\n{INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})

    required_columns = {
        "sample", "read_technology", "aligner",
        "total_length", "bases_mapped_cigar", "reads_unmapped",
    }
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    plot_data = data.loc[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER),
        ["sample", "read_technology", "aligner", "total_length", "bases_mapped_cigar", "reads_unmapped"],
    ].copy()

    for column in ["total_length", "bases_mapped_cigar", "reads_unmapped"]:
        plot_data[column] = pd.to_numeric(plot_data[column], errors="raise")

    duplicated = plot_data.duplicated(subset=["sample", "read_technology", "aligner"], keep=False)
    if duplicated.any():
        raise ValueError("Duplicated sample/technology/aligner rows found:\n" + plot_data.loc[duplicated].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} sample x technology x aligner rows, found {len(plot_data)}.")

    # -----------------------------------------------------------
    # Verify raw input bases per sample/technology.
    #
    # Only aligners that retain unmapped reads in their output
    # (reads_unmapped > 0) report a total_length equal to the true
    # FASTQ base count. An aligner with reads_unmapped == 0 across
    # every row has dropped unmapped reads from its output entirely,
    # so its own total_length cannot be trusted as the denominator.
    # -----------------------------------------------------------

    retaining = plot_data.loc[plot_data["reads_unmapped"] > 0]

    print()
    print("Raw input verification (per sample/technology):")

    verified_raw_input = {}

    for (sample, technology), group in plot_data.groupby(["sample", "read_technology"], observed=True):
        retaining_group = retaining.loc[
            (retaining["sample"] == sample) & (retaining["read_technology"] == technology)
        ]

        if retaining_group.empty:
            raise ValueError(
                f"Cannot verify raw input bases for {sample}/{technology}: "
                "no aligner in this group retains unmapped reads "
                "(reads_unmapped > 0 for none of them). Missing: an "
                "independent FASTQ base count for this sample/technology."
            )

        distinct_totals = retaining_group["total_length"].unique()
        if len(distinct_totals) != 1:
            raise ValueError(
                f"Cannot verify raw input bases for {sample}/{technology}: "
                f"aligners that retain unmapped reads disagree on total_length: "
                + retaining_group[["aligner", "total_length", "reads_unmapped"]].to_string(index=False)
            )

        verified_raw_input[(sample, technology)] = int(distinct_totals[0])

        for _, row in group.sort_values("aligner").iterrows():
            flag = "retains unmapped reads" if row["reads_unmapped"] > 0 else "DROPS unmapped reads (excluded from verification)"
            match = "matches verified total" if row["total_length"] == distinct_totals[0] else "own total_length differs -- overridden by verified total"
            print(
                f"  {sample} {technology:6s} {row['aligner']:11s} "
                f"total_length={int(row['total_length']):>12d}  reads_unmapped={int(row['reads_unmapped']):>7d}  "
                f"({flag}, {match})"
            )

    plot_data["raw_input_bases"] = plot_data.apply(
        lambda row: verified_raw_input[(row["sample"], row["read_technology"])],
        axis=1,
    )

    # -----------------------------------------------------------
    # Compute aligned / unaligned bases against the verified raw input.
    # -----------------------------------------------------------

    plot_data["cigar_aligned_bases"] = plot_data["bases_mapped_cigar"]
    plot_data["cigar_unaligned_bases"] = plot_data["raw_input_bases"] - plot_data["cigar_aligned_bases"]

    if (plot_data["cigar_unaligned_bases"] < 0).any():
        raise ValueError(
            "cigar_unaligned_bases is negative for at least one row -- "
            "bases_mapped_cigar exceeds the verified raw input bases:\n"
            + plot_data.loc[plot_data["cigar_unaligned_bases"] < 0].to_string(index=False)
        )

    reconstructed = plot_data["cigar_aligned_bases"] + plot_data["cigar_unaligned_bases"]
    mismatch = (reconstructed - plot_data["raw_input_bases"]).abs() > ROUNDING_TOLERANCE_BASES
    if mismatch.any():
        raise ValueError(
            "cigar_aligned_bases + cigar_unaligned_bases != raw_input_bases "
            "for at least one row:\n" + plot_data.loc[mismatch].to_string(index=False)
        )

    # Raw input must be identical across all aligners within a sample/technology.
    inconsistent = plot_data.groupby(["sample", "read_technology"], observed=True)["raw_input_bases"].nunique()
    if (inconsistent != 1).any():
        raise ValueError(
            "raw_input_bases is not constant within a sample/technology group:\n"
            + inconsistent.loc[inconsistent != 1].to_string()
        )

    print()
    print("Raw input bases are verified and consistent for every sample/technology group.")

    plot_data["raw_input_gb"] = plot_data["raw_input_bases"] / 1e9
    plot_data["cigar_aligned_gb"] = plot_data["cigar_aligned_bases"] / 1e9
    plot_data["cigar_unaligned_gb"] = plot_data["cigar_unaligned_bases"] / 1e9

    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    plot_data["aligner"] = pd.Categorical(plot_data["aligner"], ALIGNER_ORDER, ordered=True)

    output_columns = [
        "sample", "read_technology", "aligner",
        "raw_input_bases", "raw_input_gb",
        "cigar_aligned_bases", "cigar_aligned_gb",
        "cigar_unaligned_bases", "cigar_unaligned_gb",
    ]

    return plot_data.sort_values(["read_technology", "aligner", "sample"])[output_columns]


def main() -> int:
    plot_data = load_and_verify_raw_input()

    if len(plot_data) != len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER):
        raise ValueError(f"Expected 24 plotted combinations, found {len(plot_data)}.")

    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    print()
    print("Plot data:")
    print(plot_data.to_string(index=False, float_format=lambda value: f"{value:.4f}"))

    apply_style()

    x_positions = np.arange(len(ALIGNER_ORDER))
    figure, axes = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(FULL_WIDTH_IN, 5.0))

    y_max_data = plot_data["raw_input_gb"].max()
    y_min_display = 80.0
    y_max_display = 100.0
    label_font_size = 8.0
    # All raw-input total labels share one baseline just above the tallest
    # bar in the figure, rotated vertically over their own bar, so they
    # read as a tidy row instead of being staggered around each bar top.
    label_baseline = y_max_data + (y_max_display - y_min_display) * 0.015

    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.09, right=0.98, bottom=0.10, top=0.86, wspace=0.10)

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        technology_data = plot_data[plot_data["read_technology"] == technology]

        for sample in SAMPLE_ORDER:
            sample_data = (
                technology_data[technology_data["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
            )
            aligned_values = sample_data["cigar_aligned_gb"].to_numpy(dtype=float)
            unaligned_values = sample_data["cigar_unaligned_gb"].to_numpy(dtype=float)
            raw_values = sample_data["raw_input_gb"].to_numpy(dtype=float)
            positions = x_positions + SAMPLE_OFFSETS[sample]

            base_color = SAMPLE_COLORS[sample]
            light_color = lighten_color(base_color, amount=0.72)

            axis.bar(
                positions, aligned_values, width=BAR_WIDTH,
                color=base_color, edgecolor="white", linewidth=0.6, zorder=3,
            )
            axis.bar(
                positions, unaligned_values, width=BAR_WIDTH, bottom=aligned_values,
                color=light_color, edgecolor="white", linewidth=0.6, zorder=3,
            )

            for position, raw_value in zip(positions, raw_values):
                if np.isnan(raw_value):
                    continue
                label = f"{raw_value:.1f} Gb"
                axis.text(
                    position, label_baseline, label,
                    ha="center", va="bottom", rotation=90, fontsize=label_font_size,
                    color=darken_color(base_color, amount=0.15), zorder=7,
                )

        axis.set_title(TECHNOLOGY_TITLES[technology], pad=6)
        panel_letter(axis, "a" if technology == "ONT" else "b", fontsize=12)
        axis.set_xticks(x_positions, ALIGNER_ORDER, fontsize=12, rotation=45)
        axis.set_xlim(-0.6, len(ALIGNER_ORDER) - 0.4)
        axis.set_ylim(y_min_display, y_max_display)
        axis.tick_params(axis="y", labelsize=12)
        axis.grid(axis="y", color="#E5E5E5", linewidth=0.45, zorder=0)
        axis.set_axisbelow(True)
        for spine_name, spine in axis.spines.items():
            spine.set_visible(spine_name in ("left", "bottom"))
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.set_facecolor("white")

    axes[0].set_ylabel("Base count (Gb)", fontsize=13)

    legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="white", label=sample) for sample in SAMPLE_ORDER
    ] + [
        Patch(facecolor="#555555", edgecolor="white", label="CIGAR-aligned bases"),
        Patch(facecolor="#DDDDDD", edgecolor="white", label="CIGAR-unaligned bases"),
    ]
    figure.legend(
        handles=legend_handles, frameon=False, ncols=5, loc="upper center",
        bbox_to_anchor=(0.5, 1.01), columnspacing=1.3, handletextpad=0.4, fontsize=12,
    )

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print()
    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_SUMMARY}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
