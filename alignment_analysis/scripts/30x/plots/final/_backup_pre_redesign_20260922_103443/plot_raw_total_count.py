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

INPUT_XLSX = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "final/alignment_benchmark_30x.xlsx"
)

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
    / "final/10_raw_base_counts_30x.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/10_raw_base_counts_30x.pdf"
)

OUTPUT_SUMMARY = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/raw_base_counts_30x.tsv"
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
    "VACMap",
    "VG Giraffe",
]

technology_order = [
    "ONT",
    "PacBio",
]

configuration_to_aligner = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
    "vacmap-ont": "VACMap",
    "vacmap-pb": "VACMap",
    "vg-ont": "VG Giraffe",
    "vg-pb": "VG Giraffe",
}

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

if INPUT_XLSX.exists():
    print()
    print("Input file:")
    print(INPUT_XLSX)
    data = pd.read_excel(INPUT_XLSX)

elif INPUT_TSV.exists():
    print()
    print("Input file:")
    print(INPUT_TSV)
    data = pd.read_csv(INPUT_TSV, sep="\t")

else:
    raise FileNotFoundError(
        "Neither the requested Excel file nor the equivalent TSV was found:\n"
        f"{INPUT_XLSX}\n{INPUT_TSV}"
    )

print()
print("Input dimensions:")
print(data.shape)


# ============================================================
# STEP 5: REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "configuration",
    "total_length",
    "bases_mapped_cigar",
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
# STEP 6: SELECT BENCHMARK ROWS
# ============================================================

selected_conditions = (
    (
        (data["read_technology"] == "ONT")
        & data["configuration"].isin([
            "mm2-ont",
            "pbmm2-ont",
            "vacmap-ont",
            "vg-ont",
        ])
    )
    |
    (
        (data["read_technology"] == "PacBio")
        & data["configuration"].isin([
            "mm2-pb",
            "pbmm2-pb",
            "vacmap-pb",
            "vg-pb",
        ])
    )
)

plot_data = data.loc[
    selected_conditions,
    [
        "sample",
        "read_technology",
        "configuration",
        "total_length",
        "bases_mapped_cigar",
    ] + (
        ["input_fastq_total_bases"]
        if "input_fastq_total_bases" in data.columns
        else []
    ),
].copy()

plot_data["aligner"] = plot_data["configuration"].map(
    configuration_to_aligner
)


# ============================================================
# STEP 7: CONVERT TO NUMERIC
# ============================================================

numeric_columns = [
    "total_length",
    "bases_mapped_cigar",
]

for column in numeric_columns:

    plot_data[column] = pd.to_numeric(
        plot_data[column],
        errors="raise",
    )

if (
    plot_data["total_length"] <= 0
).any():

    raise ValueError(
        "At least one total_length value is <= 0."
    )


# ============================================================
# STEP 8: INPUT-BASES DENOMINATOR
# ============================================================

uses_fastq_bases = "input_fastq_total_bases" in plot_data.columns

if uses_fastq_bases:

    plot_data["input_fastq_total_bases"] = pd.to_numeric(
        plot_data["input_fastq_total_bases"],
        errors="raise",
    )

    if plot_data["input_fastq_total_bases"].isna().any():
        raise ValueError(
            "input_fastq_total_bases contains missing values."
        )

    input_bases = plot_data["input_fastq_total_bases"]

    print()
    print("Using actual input_fastq_total_bases as denominator.")

else:

    # Temporary diagnostic proxy only. Replace this with the actual
    # original FASTQ base count for the final publication calculation.
    input_bases = (
        plot_data
        .groupby(["sample", "read_technology"])["total_length"]
        .transform("max")
    )

    print()
    print(
        "WARNING: input_fastq_total_bases is unavailable. Using the "
        "maximum total_length across aligners as a temporary proxy."
    )

plot_data["input_bases_denominator"] = input_bases


# ============================================================
# STEP 9: CALCULATE RAW BASE COUNTS (Gb)
# ============================================================

plot_data["aligned_bases_gb"] = (
    plot_data["bases_mapped_cigar"]
    / 1e9
)

plot_data["unaligned_bases_gb"] = (
    (
        plot_data["input_bases_denominator"]
        - plot_data["bases_mapped_cigar"]
    )
    / 1e9
)

plot_data["total_bases_gb"] = (
    plot_data["aligned_bases_gb"]
    + plot_data["unaligned_bases_gb"]
)


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

expected_rows = (
    len(sample_order)
    * len(technology_order)
    * len(aligner_order)
)

if len(plot_data) != expected_rows:

    raise ValueError(
        f"Expected {expected_rows} technology-matched rows, "
        f"but found {len(plot_data)}."
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

panel_letters = {
    "ONT": "a",
    "PacBio": "b",
}

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

    y_maximum = (
        technology_data["total_bases_gb"].max()
        * 1.14
    )

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
                "aligned_bases_gb"
            ]
            .to_numpy(
                dtype=float
            )
        )

        unaligned_values = (
            sample_data[
                "unaligned_bases_gb"
            ]
            .to_numpy(
                dtype=float
            )
        )

        total_values = (
            sample_data[
                "total_bases_gb"
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

        unaligned_color = lighten_color(
            base_color,
            amount=0.72,
        )

        # ====================================================
        # ALIGNED SEGMENT
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
        # UNALIGNED SEGMENT
        # ====================================================

        unaligned_bars = axis.bar(
            positions,
            unaligned_values,
            width=bar_width,
            bottom=aligned_values,
            color=unaligned_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )

        # ====================================================
        # UNALIGNED LABELS
        #
        # Centered inside each unaligned segment, mirroring how
        # the mismatch label sits inside its own thin segment in
        # the CIGAR-composition figure. This keeps the value
        # anchored to its bar instead of colliding with the
        # total-bases label placed above the stack.
        # ====================================================

        for (
            unaligned_bar,
            unaligned_value,
        ) in zip(
            unaligned_bars,
            unaligned_values,
        ):

            if np.isnan(
                unaligned_value
            ):
                continue

            x_center = (
                unaligned_bar.get_x()
                + unaligned_bar.get_width()
                / 2
            )

            y_top = (
                unaligned_bar.get_y()
                + unaligned_bar.get_height()
            )

            axis.text(
                x_center,
                y_top + y_maximum * 0.010,
                f"{unaligned_value:.2f} Gb",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="normal",
                color="#444444",
                rotation=0,
                clip_on=True,
                zorder=7,
            )

        # ====================================================
        # TOTAL BASES LABEL ABOVE BARS
        # ====================================================

        # Anchored to each bar's own stack height (not a shared
        # panel-wide tier) so it always clears the unaligned-Gb
        # label directly beneath it; a per-sample step still keeps
        # neighbouring HG002/HG003/HG004 labels from colliding.
        total_tier_gap = (
            y_maximum
            * (
                0.060
                + 0.040
                * sample_order.index(sample)
            )
        )

        for (
            position,
            total_value,
        ) in zip(
            positions,
            total_values,
        ):

            if np.isnan(
                total_value
            ):
                continue

            axis.text(
                position,
                total_value + total_tier_gap,
                f"{total_value:.1f} Gb",
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
        0,
        y_maximum * 1.12,
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
# STEP 15: GIAB SAMPLE LEGEND
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
# STEP 16: COMPONENT LEGEND
# ============================================================

component_legend_handles = [

    Patch(
        facecolor="#555555",
        edgecolor="white",
        label="CIGAR-aligned bases",
    ),

    Patch(
        facecolor="#DDDDDD",
        edgecolor="white",
        label="Unaligned bases",
    ),
]

figure.legend(
    handles=component_legend_handles,
    frameon=False,
    ncols=2,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.865,
    ),
    columnspacing=1.8,
    fontsize=16,
)


# ============================================================
# STEP 17: NOTE
# ============================================================

denominator_note = (
    "Denominator: actual FASTQ bases"
    if uses_fastq_bases
    else "TEMPORARY DIAGNOSTIC PROXY: denominator = max(total_length) per sample/technology"
)

figure.text(
    0.5,
    0.815,

    (
        "Unaligned-segment labels (Gb) sit above their bar; "
        "labels above each stack show total input bases. "
        f"{denominator_note}."
    ),

    ha="center",
    va="center",

    fontsize=9,

    color="#555555",
)


# ============================================================
# STEP 18: FINAL LAYOUT
# ============================================================

figure.patch.set_facecolor(
    "white"
)

figure.subplots_adjust(
    left=0.09,
    right=0.98,
    bottom=0.077,
    top=0.755,
    hspace=0.18,
)


# ============================================================
# STEP 18b: SHARED Y-AXIS LABEL AND PANEL LETTERS
#
# The shared y-axis label sits close to the axis (immediately
# left of the tick labels). The "a"/"b" panel letters sit inside
# each panel's own top-left corner. The y-axis label position is
# measured from the actual rendered tick-label extents so it
# never overlaps regardless of font/DPI.
# ============================================================

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
    "Raw base count (Gb)",
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
# STEP 19: SAVE
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
# STEP 20: CONFIRM
# ============================================================

print()

print(
    "Raw base count stacked figure created successfully."
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
