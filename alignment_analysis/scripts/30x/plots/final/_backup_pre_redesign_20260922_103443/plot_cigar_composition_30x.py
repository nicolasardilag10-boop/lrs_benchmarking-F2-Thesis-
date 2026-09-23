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
    / "final/alignment_benchmark_30x.tsv"
)

OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/06_cigar_composition_30x.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/06_cigar_composition_30x.pdf"
)

OUTPUT_SUMMARY = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/cigar_mapped_mismatch_unaligned_percent.tsv"
)

OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_SUMMARY.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# STEP 2: HELPER FUNCTION
# ============================================================

def lighten_color(color, amount=0.5):
    """
    Mix a color with white.

    amount = 0.0 -> original color
    amount = 1.0 -> white
    """

    rgb = np.array(
        mcolors.to_rgb(color)
    )

    white = np.array(
        [1.0, 1.0, 1.0]
    )

    return tuple(
        rgb + (white - rgb) * amount
    )


# ============================================================
# STEP 3: ORDER AND SAMPLE COLORS
# ============================================================

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


sample_colors = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}


# Slightly larger spacing between samples.
sample_offsets = {
    "HG002": -0.23,
    "HG003": 0.00,
    "HG004": 0.23,
}

# Wide enough for readable in-bar labels while keeping samples separated.
bar_width = 0.19


# ============================================================
# STEP 4: LOAD TABLE
# ============================================================

if not INPUT_TSV.exists():

    raise FileNotFoundError(
        f"Input TSV not found:\n{INPUT_TSV}"
    )


data = pd.read_csv(
    INPUT_TSV,
    sep="\t",
)


print()
print("Input file:")
print(INPUT_TSV)

print()
print("Input dimensions:")
print(data.shape)


# ============================================================
# STEP 5: REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "total_length",
    "bases_mapped_cigar",
    "mismatches",
}


missing_columns = (
    required_columns
    .difference(
        data.columns
    )
)


if missing_columns:

    raise ValueError(
        "Missing required columns: "
        f"{sorted(missing_columns)}"
    )


# ============================================================
# STEP 6: CONVERT TO NUMERIC
# ============================================================

numeric_columns = [
    "total_length",
    "bases_mapped_cigar",
    "mismatches",
]


for column in numeric_columns:

    data[column] = pd.to_numeric(
        data[column],
        errors="raise",
    )


if (
    data["total_length"] <= 0
).any():

    raise ValueError(
        "At least one total_length value is <= 0."
    )


# ============================================================
# STEP 7: CALCULATE CIGAR METRICS
# ============================================================

# ------------------------------------------------------------
# Total CIGAR-mapped percentage
# ------------------------------------------------------------

data["cigar_mapped_percent"] = (
    data["bases_mapped_cigar"]
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# CIGAR-unaligned percentage
# ------------------------------------------------------------

data["cigar_unaligned_percent"] = (
    (
        data["total_length"]
        - data["bases_mapped_cigar"]
    )
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# Mismatches as percentage of ALL input bases
# ------------------------------------------------------------

data["mismatch_total_percent"] = (
    data["mismatches"]
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# Mismatch percentage within CIGAR-mapped sequence
# ------------------------------------------------------------

data["mismatch_within_cigar_percent"] = np.where(
    data["bases_mapped_cigar"] > 0,

    (
        data["mismatches"]
        / data["bases_mapped_cigar"]
        * 100.0
    ),

    np.nan,
)


# ------------------------------------------------------------
# Aligned non-mismatch percentage
# ------------------------------------------------------------

data["aligned_non_mismatch_percent"] = (
    (
        data["bases_mapped_cigar"]
        - data["mismatches"]
    )
    / data["total_length"]
    * 100.0
)


# ============================================================
# STEP 8: VERIFY STACK SUMS TO 100%
# ============================================================

data["stack_total_percent"] = (
    data["aligned_non_mismatch_percent"]
    + data["mismatch_total_percent"]
    + data["cigar_unaligned_percent"]
)


maximum_stack_error = (
    (
        data["stack_total_percent"]
        - 100.0
    )
    .abs()
    .max()
)


print()
print(
    "Maximum deviation of stacked components from 100%:"
)

print(
    f"{maximum_stack_error:.10f}"
)


if maximum_stack_error > 0.001:

    raise ValueError(
        "Stacked components do not sum to 100%."
    )


# ============================================================
# STEP 9: SELECT BENCHMARK ROWS
# ============================================================

plot_data = data.loc[
    data["sample"].isin(
        sample_order
    )
    & data["read_technology"].isin(
        technology_order
    )
    & data["aligner"].isin(
        aligner_order
    ),
    [
        "sample",
        "read_technology",
        "aligner",
        "total_length",
        "bases_mapped_cigar",
        "mismatches",
        "aligned_non_mismatch_percent",
        "mismatch_total_percent",
        "mismatch_within_cigar_percent",
        "cigar_mapped_percent",
        "cigar_unaligned_percent",
    ],
].copy()


# ============================================================
# STEP 10: CHECK DUPLICATES
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
        "Duplicated sample/technology/aligner rows found:\n"
        + plot_data.loc[
            duplicated_rows
        ].to_string(
            index=False
        )
    )


# ============================================================
# STEP 11: CATEGORICAL ORDER
# ============================================================

plot_data["sample"] = pd.Categorical(
    plot_data["sample"],
    categories=sample_order,
    ordered=True,
)


plot_data["read_technology"] = pd.Categorical(
    plot_data["read_technology"],
    categories=technology_order,
    ordered=True,
)


plot_data["aligner"] = pd.Categorical(
    plot_data["aligner"],
    categories=aligner_order,
    ordered=True,
)


plot_data = plot_data.sort_values(
    [
        "read_technology",
        "aligner",
        "sample",
    ]
)


# ============================================================
# STEP 12: SAVE CALCULATED METRICS
# ============================================================

plot_data.to_csv(
    OUTPUT_SUMMARY,
    sep="\t",
    index=False,
)


print()
print("Calculated values:")

print(
    plot_data.to_string(
        index=False,
        float_format=lambda value: f"{value:.4f}",
    )
)


# ============================================================
# STEP 13: CREATE FIGURE
# ============================================================

figure, axes = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    sharey=True,
    figsize=(
        13.0,
        9.5,
    ),
)


x_positions = np.arange(
    len(
        aligner_order
    )
)


# ============================================================
# STEP 14: DRAW STACKED BARS
# ============================================================

for panel_index, technology in enumerate(
    technology_order
):

    axis = axes[
        panel_index
    ]


    technology_data = plot_data.loc[
        plot_data[
            "read_technology"
        ]
        == technology
    ]


    for sample in sample_order:

        sample_data = (
            technology_data.loc[
                technology_data[
                    "sample"
                ]
                == sample
            ]
            .set_index(
                "aligner"
            )
            .reindex(
                aligner_order
            )
        )


        # ====================================================
        # VALUES
        # ====================================================

        aligned_values = (
            sample_data[
                "aligned_non_mismatch_percent"
            ]
            .to_numpy(
                dtype=float
            )
        )


        mismatch_values = (
            sample_data[
                "mismatch_total_percent"
            ]
            .to_numpy(
                dtype=float
            )
        )


        unaligned_values = (
            sample_data[
                "cigar_unaligned_percent"
            ]
            .to_numpy(
                dtype=float
            )
        )


        cigar_mapped_values = (
            sample_data[
                "cigar_mapped_percent"
            ]
            .to_numpy(
                dtype=float
            )
        )


        positions = (
            x_positions
            + sample_offsets[
                sample
            ]
        )


        base_color = (
            sample_colors[
                sample
            ]
        )


        mismatch_color = lighten_color(
            base_color,
            amount=0.35,
        )


        unaligned_color = lighten_color(
            base_color,
            amount=0.72,
        )


        # ====================================================
        # ALIGNED NON-MISMATCH SEGMENT
        # ====================================================

        axis.bar(
            positions,
            aligned_values,
            width=bar_width,
            color=base_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )


        # ====================================================
        # MISMATCH SEGMENT
        # ====================================================

        mismatch_bars = axis.bar(
            positions,
            mismatch_values,
            width=bar_width,
            bottom=aligned_values,
            color=mismatch_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )


        # ====================================================
        # CIGAR-UNALIGNED SEGMENT
        # ====================================================

        axis.bar(
            positions,
            unaligned_values,
            width=bar_width,
            bottom=(
                aligned_values
                + mismatch_values
            ),
            color=unaligned_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )


        # ====================================================
        # MISMATCH LABELS
        #
        # Put each mismatch percentage horizontally in the
        # centre of its own mismatch segment.  This keeps labels
        # associated with the correct bar and avoids the diagonal
        # overlap seen in the previous version.
        # ====================================================

        for (
            mismatch_bar,
            mismatch_value,
        ) in zip(
            mismatch_bars,
            mismatch_values,
        ):

            if np.isnan(
                mismatch_value
            ):
                continue

            x_center = (
                mismatch_bar.get_x()
                + mismatch_bar.get_width()
                / 2
            )

            y_center = (
                mismatch_bar.get_y()
                + mismatch_bar.get_height()
                / 2
            )

            # Use a slightly smaller font only for the thinnest
            # mismatch bands (e.g. ~0.6-0.8% in PacBio HiFi).
            mismatch_fontsize = (
                9.5
                if mismatch_value < 0.85
                else 10.5
            )

            axis.text(
                x_center,
                y_center,
                f"{mismatch_value:.2f}%",
                ha="center",
                va="center",
                fontsize=mismatch_fontsize,
                fontweight="normal",
                color="black",
                rotation=0,
                clip_on=True,
                zorder=7,
            )


        # ====================================================
        # CIGAR-MAPPED PERCENTAGE ABOVE BARS
        # ====================================================

        # Use three compact horizontal label tiers above the
        # complete 100% stacks.  The tiering prevents neighbouring
        # HG002/HG003/HG004 values from colliding.
        mapped_label_y = {
            "HG002": 100.18,
            "HG003": 100.55,
            "HG004": 100.92,
        }[sample]


        for (
            position,
            mapped_value,
        ) in zip(
            positions,
            cigar_mapped_values,
        ):

            if np.isnan(
                mapped_value
            ):
                continue

            axis.text(
                position,
                mapped_label_y,
                f"{mapped_value:.2f}%",
                ha="center",
                va="bottom",
                fontsize=9.5,
                fontweight="bold",
                color=base_color,
                rotation=0,
                clip_on=False,
                zorder=8,
            )


    # ========================================================
    # PANEL TITLE
    # ========================================================

    if technology == "ONT":

        panel_title = "ONT"

    else:

        panel_title = "PacBio HiFi"


    axis.set_title(
        panel_title,
        fontsize=16,
        pad=14,
    )


    # ========================================================
    # X AXIS
    # ========================================================

    axis.set_xticks(
        x_positions,
        aligner_order,
    )


    if technology != technology_order[-1]:

        axis.tick_params(
            axis="x",
            labelbottom=False,
        )


    # ========================================================
    # Y AXIS
    # ========================================================

    axis.set_ylim(
        88,
        102.1,
    )


    axis.set_yticks(
        [
            88,
            90,
            92,
            94,
            96,
            98,
            100,
        ]
    )


    # ========================================================
    # GRID / STYLE
    # ========================================================

    axis.grid(
        axis="y",
        color="#D9D9D9",
        linewidth=0.7,
        linestyle="-",
        alpha=0.75,
        zorder=0,
    )


    axis.spines[
        "top"
    ].set_visible(
        False
    )


    axis.spines[
        "right"
    ].set_visible(
        False
    )


    axis.spines[
        "left"
    ].set_color(
        "black"
    )


    axis.spines[
        "bottom"
    ].set_color(
        "black"
    )


    axis.tick_params(
        axis="both",
        labelsize=11,
        colors="black",
    )


    # This figure's canvas (figsize=(13.0, 9.5)) is wider than the
    # other benchmark figures (~9.1 in), so once every figure is
    # scaled to the same \linewidth in the report, an 11pt aligner
    # label here renders visibly smaller than the matching 11pt
    # labels elsewhere. Bump just the x-tick (aligner name) size to
    # compensate for that extra scale-down.
    axis.tick_params(
        axis="x",
        labelsize=14.7,
        pad=2,
    )


    axis.set_facecolor(
        "white"
    )



# ============================================================
# STEP 16: GIAB SAMPLE LEGEND
# ============================================================

sample_legend_handles = [

    Patch(
        facecolor=sample_colors[
            sample
        ],
        edgecolor="white",
        label=sample,
    )

    for sample in sample_order
]


figure.legend(
    handles=sample_legend_handles,
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.992,
    ),
    columnspacing=1.8,
    fontsize=17,
)


# ============================================================
# STEP 17: COMPONENT LEGEND
# ============================================================

component_legend_handles = [

    Patch(
        facecolor="#555555",
        edgecolor="white",
        label="Aligned, non-mismatch",
    ),

    Patch(
        facecolor="#999999",
        edgecolor="white",
        label="Mismatch bases",
    ),

    Patch(
        facecolor="#DDDDDD",
        edgecolor="white",
        label="CIGAR-unaligned bases",
    ),
]


figure.legend(
    handles=component_legend_handles,
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.865,
    ),
    columnspacing=1.8,
    fontsize=16,
)


# ============================================================
# STEP 18: NOTE
# ============================================================

figure.text(
    0.5,
    0.815,

    (
        "Mismatch percentages are centered inside mismatch segments; "
        "labels above bars show total CIGAR-mapped bases"
    ),

    ha="center",
    va="center",

    fontsize=9,

    color="#555555",
)


# ============================================================
# STEP 19: FINAL LAYOUT
# ============================================================

figure.patch.set_facecolor(
    "white"
)


figure.subplots_adjust(
    left=0.135,
    right=0.98,
    bottom=0.077,
    top=0.755,
    hspace=0.18,
)


# ============================================================
# STEP 19b: SHARED Y-AXIS LABEL AND PANEL LETTERS
#
# The shared y-axis label sits close to the axis (immediately
# left of the tick labels). The "a"/"b" panel letters sit inside
# each panel's own top-left corner. The y-axis label position is
# measured from the actual rendered tick-label extents so it
# never overlaps regardless of font/DPI.
# ============================================================

panel_letters = {
    "ONT": "a",
    "PacBio": "b",
}

renderer = figure.canvas.get_renderer()

inches_to_x_fraction = 1.0 / figure.get_figwidth()
label_gap = 0.05 * inches_to_x_fraction

tick_label_left_edges = []

for axis in axes:

    for tick_label in axis.get_yticklabels():

        bbox = tick_label.get_window_extent(
            renderer=renderer
        )

        tick_label_left_edges.append(
            bbox.x0
        )

tick_labels_left_display = min(
    tick_label_left_edges
)

tick_labels_left_fraction = figure.transFigure.inverted().transform(
    (tick_labels_left_display, 0)
)[0]

target_ylabel_right_fraction = (
    tick_labels_left_fraction
    - label_gap
)

# supylabel defaults to vertically centring on the full figure, but
# the legends above the axes push the axes block itself off-centre.
# Anchor to the axes block's own vertical centre instead, so a larger
# font does not push the label up into the legend area.
axes_block_top = axes[0].get_position().y1
axes_block_bottom = axes[-1].get_position().y0
axes_block_center = (axes_block_top + axes_block_bottom) / 2

# supylabel's `x` sets the horizontal centre of the (rotated) text
# bounding box, not its right edge, so place it once, measure the
# rendered right edge, then correct by the measured offset.
ylabel_artist = figure.supylabel(
    "CIGAR-aligned and unaligned bases (%)",
    fontsize=17,
    x=target_ylabel_right_fraction,
    y=axes_block_center,
)

figure.canvas.draw()
renderer = figure.canvas.get_renderer()

ylabel_right_fraction = figure.transFigure.inverted().transform(
    (
        ylabel_artist.get_window_extent(
            renderer=renderer
        ).x1,
        0,
    )
)[0]

ylabel_artist.set_x(
    target_ylabel_right_fraction
    - (
        ylabel_right_fraction
        - target_ylabel_right_fraction
    )
)

figure.canvas.draw()
renderer = figure.canvas.get_renderer()

for panel_index, technology in enumerate(
    technology_order
):

    axes[panel_index].text(
        -0.08,
        1.02,
        panel_letters[technology],
        transform=axes[panel_index].transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",
    )


# ============================================================
# STEP 20: SAVE
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


plt.close(
    figure
)


# ============================================================
# STEP 21: CONFIRM
# ============================================================

print()

print(
    "CIGAR stacked figure created successfully."
)

print()

print(
    "PNG:"
)
print(
    OUTPUT_PNG
)

print()

print(
    "PDF:"
)
print(
    OUTPUT_PDF
)

print()

print(
    "Calculated metrics TSV:"
)
print(
    OUTPUT_SUMMARY
)
