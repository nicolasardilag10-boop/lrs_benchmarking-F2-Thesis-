#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
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
    / "final/04_mq0_reads_30x.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/04_mq0_reads_30x.pdf"
)

OUTPUT_SUMMARY = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/mq0_reads_percent.tsv"
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
# STEP 2: ORDERS
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


# ============================================================
# STEP 3: SAMPLE COLORS
# ============================================================

sample_colors = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}


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
# STEP 5: STANDARDIZE ALIGNER NAMES
# ============================================================

data["aligner"] = (
    data["aligner"]
    .replace(
        {
            "VACMap": "VACmap",
        }
    )
)


# ============================================================
# STEP 6: REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "reads_mq0",
    "reads_mapped",
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
# STEP 7: NUMERIC CONVERSION
# ============================================================

numeric_columns = [
    "reads_mq0",
    "reads_mapped",
]


for column in numeric_columns:

    data[column] = pd.to_numeric(
        data[column],
        errors="raise",
    )


# ============================================================
# STEP 8: CALCULATE MQ0 PERCENT
# ============================================================

data["reads_mq0_percent"] = np.where(

    data["reads_mapped"] > 0,

    (
        data["reads_mq0"]
        / data["reads_mapped"]
        * 100.0
    ),

    np.nan,
)


# ============================================================
# STEP 9: SELECT BENCHMARK ROWS
# ============================================================

plot_data = data.loc[
    data["sample"].isin(sample_order)
    & data["read_technology"].isin(technology_order)
    & data["aligner"].isin(aligner_order),
    [
        "sample",
        "read_technology",
        "aligner",
        "reads_mq0",
        "reads_mapped",
        "reads_mq0_percent",
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
# STEP 11: CHECK COMPLETENESS
# ============================================================

expected_rows = (
    len(sample_order)
    * len(technology_order)
    * len(aligner_order)
)


print()
print(
    f"Expected benchmark rows: {expected_rows}"
)

print(
    f"Observed benchmark rows: {len(plot_data)}"
)


if len(plot_data) != expected_rows:

    print()
    print(
        "WARNING: dataset is not complete."
    )


# ============================================================
# STEP 12: CATEGORICAL ORDER
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
# STEP 13: SAVE CALCULATED VALUES
# ============================================================

plot_data.to_csv(
    OUTPUT_SUMMARY,
    sep="\t",
    index=False,
)


print()
print("MQ0 values:")

print(
    plot_data.to_string(
        index=False,
        float_format=lambda value: f"{value:.4f}",
    )
)


# ============================================================
# STEP 14: FIGURE SETTINGS
# ============================================================

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(
        11.2,
        5.46,
    ),
)


x_positions = np.arange(
    len(aligner_order))

# Narrower bars
bar_width = 0.18

# Slightly larger separation between the bars inside each aligner group
sample_spacing = 0.21

sample_offsets = {
    "HG002": -sample_spacing,
    "HG003": 0.0,
    "HG004": sample_spacing,}


# ============================================================
# STEP 15: SHARED Y LIMIT
# ============================================================

maximum_value = (
    plot_data[
        "reads_mq0_percent"
    ]
    .max())


# Round upward to a clean quarter-percentage interval.
y_upper = (
    np.ceil(
        maximum_value
        / 0.25
    )
    * 0.25
    + 0.25
)


print()
print(
    f"Maximum MQ0 value: {maximum_value:.4f}%"
)

print(
    f"Y-axis upper limit: {y_upper:.2f}%"
)


# ============================================================
# STEP 16: DRAW PANELS
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


    # ========================================================
    # DRAW SAMPLE BARS
    # ========================================================

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


        mq0_values = (
            sample_data[
                "reads_mq0_percent"
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


        bars = axis.bar(
            positions,
            mq0_values,

            width=bar_width,

            color=sample_colors[
                sample
            ],

            edgecolor="white",

            linewidth=1.2,

            zorder=3,
        )


        # ====================================================
        # VALUE LABELS
        # ====================================================

        #for bar, value in zip(
        #    bars,
        #    mq0_values,
        #):

        #    if np.isnan(
        #        value
        #    ):
        #        continue


        #    axis.text(
        #        bar.get_x()
        #        + bar.get_width()
        #        / 2,

        #        value
        #        + y_upper * 0.018,

        #        f"{value:.2f}%",

        #        ha="center",
        #        va="bottom",

        #        fontsize=10,

        #        rotation=45,

        #        rotation_mode="anchor",

        #        color="black",

        #        clip_on=False,

        #        zorder=5,
        #    )


    # ========================================================
    # PANEL TITLES
    # ========================================================

    if technology == "ONT":

        panel_title = "ONT"
        panel_letter = "a"

    else:

        panel_title = "PacBio HiFi"
        panel_letter = "b"


    axis.set_title(
        panel_title,

        fontsize=20,

        pad=20,
    )


    axis.text(
        -0.08,
        1.035,

        panel_letter,

        transform=axis.transAxes,

        ha="left",
        va="bottom",

        fontsize=20,

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
        0,
        y_upper,
    )


    # ========================================================
    # GRID
    # ========================================================

    axis.grid(
        axis="y",

        color="#D9D9D9",

        linewidth=0.7,

        linestyle="-",

        alpha=0.75,

        zorder=0,
    )


    # ========================================================
    # SPINES
    # ========================================================

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


    # ========================================================
    # TICKS
    # ========================================================

    axis.tick_params(
        axis="x",

        labelsize=13,

        colors="black",
    )


    axis.tick_params(
        axis="y",

        labelsize=12,

        colors="black",
    )


    axis.set_facecolor(
        "white"
    )


# ============================================================
# STEP 17: Y LABEL
# ============================================================

axes[0].set_ylabel(
    "MQ0 reads (%)",
    fontsize=13,)


# ============================================================
# STEP 18: LEGEND
# ============================================================

sample_legend_handles = [

    Patch(
        facecolor=sample_colors[
            sample],

        edgecolor="white",

        label=sample,
    )

    for sample in sample_order]


figure.legend(
    handles=sample_legend_handles,

    title="GIAB sample",

    frameon=False,

    ncols=3,

    loc="upper center",

    bbox_to_anchor=(
        0.5,
        0.99,),

    columnspacing=1.0,
    handletextpad=0.4,

    fontsize=9,
    title_fontsize=9,)


# ============================================================
# STEP 19: FINAL LAYOUT
# ============================================================

figure.patch.set_facecolor(
    "white")


figure.subplots_adjust(
    left=0.075,
    right=0.98,

    bottom=0.186,
    top=0.80,

    wspace=0.12,)


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
    "MQ0 grouped figure created successfully."
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
    "Summary:"
)

print(
    OUTPUT_SUMMARY
)