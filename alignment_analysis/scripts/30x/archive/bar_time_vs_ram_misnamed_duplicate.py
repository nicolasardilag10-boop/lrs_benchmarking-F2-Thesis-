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
    / "cigar_mapped_mismatch_unaligned_percent.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "cigar_mapped_mismatch_unaligned_percent.pdf"
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
    "HG002": -0.24,
    "HG003": 0.00,
    "HG004": 0.24,
}

# Slightly narrower bars.
bar_width = 0.18


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
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(
        13.0,
        7.8,
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
        # Rotate labels 45 degrees to reduce horizontal space.
        # ====================================================

        mismatch_label_y_offsets = {
            "HG002": 0.00,
            "HG003": 0.03,
            "HG004": 0.06,
        }


        mismatch_label_x_offsets = {
            "HG002": -0.006,
            "HG003":  0.000,
            "HG004":  0.006,
        }


        for (
            mismatch_bar,
            mismatch_value,
            aligned_value,
        ) in zip(
            mismatch_bars,
            mismatch_values,
            aligned_values,
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


            mismatch_bottom = (
                aligned_value
            )


            mismatch_top = (
                aligned_value
                + mismatch_value
            )


            mismatch_middle = (
                mismatch_bottom
                + mismatch_value / 2
            )


            # -----------------------------------------------
            # LARGE / MEDIUM mismatch segments
            # -----------------------------------------------

            if mismatch_value >= 0.55:

                axis.text(
                    x_center
                    + mismatch_label_x_offsets[
                        sample
                    ],

                    mismatch_middle
                    + mismatch_label_y_offsets[
                        sample
                    ],

                    f"{mismatch_value:.2f}%",

                    ha="center",
                    va="center",

                    fontsize=7.0,

                    color="black",

                    rotation=45,
                    rotation_mode="anchor",

                    clip_on=False,

                    zorder=7,
                )


            # -----------------------------------------------
            # SMALL mismatch segments
            # -----------------------------------------------

            elif mismatch_value >= 0.20:

                axis.text(
                    x_center
                    + mismatch_label_x_offsets[
                        sample
                    ],

                    mismatch_top
                    + 0.05
                    + mismatch_label_y_offsets[
                        sample
                    ],

                    f"{mismatch_value:.2f}%",

                    ha="center",
                    va="bottom",

                    fontsize=6.7,

                    color="black",

                    rotation=45,
                    rotation_mode="anchor",

                    clip_on=False,

                    zorder=7,
                )


        # ====================================================
        # CIGAR-MAPPED PERCENTAGE ABOVE BARS
        # ====================================================

        mapped_label_offsets = {
            "HG002": 0.00,
            "HG003": 0.30,
            "HG004": 0.60,
        }


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

                100.20
                + mapped_label_offsets[
                    sample
                ],

                f"{mapped_value:.2f}%",

                ha="center",
                va="bottom",

                fontsize=7.4,

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
        panel_letter = "a"

    else:

        panel_title = "PacBio HiFi"
        panel_letter = "b"


    axis.set_title(
        panel_title,
        fontsize=16,
        pad=14,
    )


    axis.text(
        -0.08,
        1.035,
        panel_letter,
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


    axis.set_xlabel(
        "Aligner",
        fontsize=13,
        labelpad=7,
    )


    # ========================================================
    # Y AXIS
    # ========================================================

    axis.set_ylim(
        88,
        102.00,
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


    axis.set_facecolor(
        "white"
    )


# ============================================================
# STEP 15: Y-AXIS LABEL
# ============================================================

axes[0].set_ylabel(
    "CIGAR-aligned and unaligned bases (%)",
    fontsize=13,
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
    title="GIAB sample",
    title_fontsize=14,
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.99,
    ),
    columnspacing=1.8,
    fontsize=11,
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
        0.905,
    ),
    columnspacing=1.8,
    fontsize=10.5,
)


# ============================================================
# STEP 18: NOTE
# ============================================================

figure.text(
    0.5,
    0.855,

    (
        "Mismatch percentages are displayed diagonally; "
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
    left=0.08,
    right=0.98,
    bottom=0.13,
    top=0.78,
    wspace=0.12,
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