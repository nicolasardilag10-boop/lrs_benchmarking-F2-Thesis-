#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch


# ============================================================
# STEP 1: PATHS
# ============================================================

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

INPUT_TSV = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv"
)

OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "mapped_unmapped_bases_percent_90_100.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "mapped_unmapped_bases_percent_90_100.pdf"
)

OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# STEP 2: HELPER FUNCTION
# ============================================================

def lighten_color(color, amount=0.55):
    """
    Return a lighter version of a color.

    amount=0   -> original color
    amount=1   -> white
    """
    rgb = np.array(mcolors.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    return tuple(rgb + (white - rgb) * amount)


# ============================================================
# STEP 3: COLORS AND PLOTTING ORDER
# ============================================================

sample_colors = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}

sample_order = [
    "HG002",
    "HG003",
    "HG004",
]

aligner_order = [
    "minimap2",
    "pbmm2",
    "VACmap",
    "VG Giraffe",
]

technology_order = [
    "ONT",
    "PacBio",
]

sample_offsets = {
    "HG002": -0.22,
    "HG003": 0.00,
    "HG004": 0.22,
}

bar_width = 0.19


# ============================================================
# STEP 4: LOAD DATA
# ============================================================

if not INPUT_TSV.exists():
    raise FileNotFoundError(
        f"Input table not found:\n{INPUT_TSV}"
    )

data = pd.read_csv(
    INPUT_TSV,
    sep="\t",
)


# ============================================================
# STEP 5: VALIDATE REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "mapped_bases_percent",
}

missing_columns = required_columns.difference(
    data.columns
)

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        f"{sorted(missing_columns)}"
    )

data["mapped_bases_percent"] = pd.to_numeric(
    data["mapped_bases_percent"],
    errors="raise",
)

# Complementary unmapped bases percentage
data["unmapped_bases_percent"] = (
    100.0 - data["mapped_bases_percent"]
)


# ============================================================
# STEP 6: SELECT RELEVANT DATA
# ============================================================

plot_data = data.loc[
    data["sample"].isin(sample_order)
    & data["read_technology"].isin(technology_order)
    & data["aligner"].isin(aligner_order),
    [
        "sample",
        "read_technology",
        "aligner",
        "mapped_bases_percent",
        "unmapped_bases_percent",
    ],
].copy()


# ============================================================
# STEP 7: CHECK FOR DUPLICATED OBSERVATIONS
# ============================================================

duplicated_rows = plot_data.duplicated(
    subset=[
        "sample",
        "read_technology",
        "aligner",
    ],
    keep=False,
)

if duplicated_rows.any():
    raise ValueError(
        "Duplicated observations found:\n"
        + plot_data.loc[
            duplicated_rows
        ].to_string(index=False)
    )


# ============================================================
# STEP 8: DEFINE CATEGORICAL ORDER
# ============================================================

plot_data["sample"] = pd.Categorical(
    plot_data["sample"],
    categories=sample_order,
    ordered=True,
)

plot_data["aligner"] = pd.Categorical(
    plot_data["aligner"],
    categories=aligner_order,
    ordered=True,
)

plot_data["read_technology"] = pd.Categorical(
    plot_data["read_technology"],
    categories=technology_order,
    ordered=True,
)

plot_data = plot_data.sort_values(
    [
        "read_technology",
        "aligner",
        "sample",
    ]
)

print()
print("Data used for plotting:")
print(plot_data.to_string(index=False))


# ============================================================
# STEP 9: CREATE FIGURE
# ============================================================

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(9.1, 5.18),
)

x_positions = np.arange(
    len(aligner_order)
)


# ============================================================
# STEP 10: DRAW STACKED MAPPED / UNMAPPED BARS
# ============================================================

for panel_index, technology in enumerate(technology_order):

    axis = axes[panel_index]

    technology_data = plot_data.loc[
        plot_data["read_technology"] == technology
    ]

    for sample in sample_order:

        sample_data = (
            technology_data.loc[
                technology_data["sample"] == sample
            ]
            .set_index("aligner")
            .reindex(aligner_order)
        )

        mapped_values = (
            sample_data["mapped_bases_percent"]
            .to_numpy(dtype=float)
        )

        unmapped_values = (
            sample_data["unmapped_bases_percent"]
            .to_numpy(dtype=float)
        )

        positions = (
            x_positions + sample_offsets[sample]
        )

        base_color = sample_colors[sample]
        light_color = lighten_color(
            base_color,
            amount=0.55,
        )

        # Mapped bases
        axis.bar(
            positions,
            mapped_values,
            width=bar_width,
            color=base_color,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )

        # Unmapped bases
        axis.bar(
            positions,
            unmapped_values,
            width=bar_width,
            bottom=mapped_values,
            color=light_color,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )

    # ========================================================
    # PANEL TITLES
    # ========================================================

    axis.set_title(
        "ONT" if technology == "ONT" else "PacBio HiFi",
        fontsize=16,
        pad=14,
    )

    axis.text(
        -0.08,
        1.03,
        "a" if technology == "ONT" else "b",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",
    )

    # ========================================================
    # X AXIS
    # ========================================================

    axis.set_xticks(
        x_positions,
        aligner_order,
    )

    # ========================================================
    # Y AXIS
    # ========================================================

    axis.set_ylim(
        90,
        100.35,
    )

    axis.set_yticks(
        [90, 92, 94, 96, 98, 100]
    )

    # ========================================================
    # GRID AND AXIS STYLE
    # ========================================================

    axis.grid(
        axis="y",
        color="#D9D9D9",
        linewidth=0.7,
        linestyle="-",
        alpha=0.75,
        zorder=0,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("black")
    axis.spines["bottom"].set_color("black")

    axis.tick_params(
        axis="both",
        labelsize=11,
        colors="black",
    )

    axis.set_facecolor("white")


# ============================================================
# STEP 11: Y AXIS LABEL
# ============================================================

axes[0].set_ylabel(
    "Mapped and unmapped bases (%)",
    fontsize=13,
)


# ============================================================
# STEP 12: GIAB SAMPLE LEGEND
# ============================================================

sample_legend_handles = [
    Patch(
        facecolor=sample_colors[sample],
        edgecolor="white",
        label=sample,
    )
    for sample in sample_order
]

figure.legend(
    handles=sample_legend_handles,
    title="GIAB sample",
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.99),
    columnspacing=1.0,
    handletextpad=0.4,
    handlelength=1.3,
    fontsize=9,
    title_fontsize=9,
)


# ============================================================
# STEP 13: MAPPED / UNMAPPED LEGEND
# ============================================================

status_legend_handles = [
    Patch(
        facecolor="#666666",
        edgecolor="white",
        label="Mapped bases",
    ),
    Patch(
        facecolor="#CCCCCC",
        edgecolor="white",
        label="Unmapped bases",
    ),
]

figure.legend(
    handles=status_legend_handles,
    frameon=False,
    ncols=2,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.90),
    columnspacing=1.2,
    handletextpad=0.4,
    handlelength=1.3,
    fontsize=9,
)


# ============================================================
# STEP 14: FINAL LAYOUT
# ============================================================

figure.patch.set_facecolor("white")

figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.186,
    top=0.74,
    wspace=0.12,
)


# ============================================================
# STEP 15: SAVE FIGURE
# ============================================================

figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
)

plt.close(figure)


# ============================================================
# STEP 16: CONFIRM SUCCESS
# ============================================================

print()
print("Mapped/unmapped bases figure created successfully.")
print()
print("Y-axis range: 90–100%")
print()
print("PNG:")
print(OUTPUT_PNG)
print()
print("PDF:")
print(OUTPUT_PDF)