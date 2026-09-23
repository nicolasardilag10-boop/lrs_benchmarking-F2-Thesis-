#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats


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


# ------------------------------------------------------------
# Existing stacked CIGAR figure
# ------------------------------------------------------------

OUTPUT_STACKED_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "cigar_mismatch_vs_error_rate.png"
)

OUTPUT_STACKED_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "cigar_mismatch_vs_error_rate.pdf"
)


# ------------------------------------------------------------
# NEW: normalized yield vs error figure
# ------------------------------------------------------------

OUTPUT_SCATTER_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "cigar_mismatch_vs_error_rate.png"
)

OUTPUT_SCATTER_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "cigar_mismatch_vs_error_rate.pdf"
)


# ------------------------------------------------------------
# NEW: normalized metrics table
# ------------------------------------------------------------

OUTPUT_SUMMARY = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/cigar_mismatch_vs_error_rate.tsv"
)


OUTPUT_STACKED_PNG.parent.mkdir(
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
# STEP 3: ORDER / COLORS
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


sample_offsets = {
    "HG002": -0.24,
    "HG003": 0.00,
    "HG004": 0.24,
}

bar_width = 0.18


# Different marker for each aligner in the scatter plot.
aligner_markers = {
    "minimap2": "o",
    "pbmm2": "s",
    "VACmap": "^",
    "VG Giraffe": "D",
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
# STEP 5: REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "total_length",
    "bases_mapped_cigar",
    "mismatches",
    "error_percent",
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
# STEP 6: CONVERT NUMERIC COLUMNS
# ============================================================

numeric_columns = [
    "total_length",
    "bases_mapped_cigar",
    "mismatches",
    "error_percent",
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
# STEP 7: STANDARD CIGAR METRICS
# ============================================================

# ------------------------------------------------------------
# CIGAR mapped percentage using each aligner's own denominator
# ------------------------------------------------------------

data["cigar_mapped_percent_internal"] = (
    data["bases_mapped_cigar"]
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# CIGAR unaligned percentage using each aligner's denominator
# ------------------------------------------------------------

data["cigar_unaligned_percent_internal"] = (
    (
        data["total_length"]
        - data["bases_mapped_cigar"]
    )
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# Mismatches as percentage of total bases
# ------------------------------------------------------------

data["mismatch_total_percent"] = (
    data["mismatches"]
    / data["total_length"]
    * 100.0
)


# ------------------------------------------------------------
# Mismatch rate within CIGAR-aligned bases
# ------------------------------------------------------------

data["mismatch_error_calculated_percent"] = np.where(
    data["bases_mapped_cigar"] > 0,

    (
        data["mismatches"]
        / data["bases_mapped_cigar"]
        * 100.0
    ),

    np.nan,
)


# ============================================================
# STEP 8: DEFINE COMMON INPUT DENOMINATOR
# ============================================================

# IMPORTANT:
#
# For each sample × technology, use one common denominator.
#
# In the current table:
#
# minimap2 / pbmm2 / VG generally report the same total_length,
# while VACmap has a smaller total_length.
#
# Using the maximum total_length therefore recovers the largest
# available representation of the original input.
#
# Best future version:
# replace this with the actual base count from the original FASTQ.

common_input_table = (
    data[
        data["sample"].isin(sample_order)
        & data["read_technology"].isin(technology_order)
        & data["aligner"].isin(aligner_order)
    ]
    .groupby(
        [
            "sample",
            "read_technology",
        ],
        observed=True,
    )["total_length"]
    .max()
    .rename("common_input_bases")
    .reset_index()
)


data = data.merge(
    common_input_table,
    on=[
        "sample",
        "read_technology",
    ],
    how="left",
)


print()
print("Common input denominators:")
print(
    common_input_table.to_string(
        index=False
    )
)


# ============================================================
# STEP 9: INPUT-NORMALIZED CIGAR YIELD
# ============================================================

# Equation:
#
# input-normalized CIGAR yield (%) =
#
# bases_mapped_cigar
# ------------------ × 100
# common_input_bases

data["input_normalized_cigar_yield_percent"] = (
    data["bases_mapped_cigar"]
    / data["common_input_bases"]
    * 100.0
)


# ============================================================
# STEP 10: INPUT-NORMALIZED UNALIGNED FRACTION
# ============================================================

data["input_normalized_cigar_unaligned_percent"] = (
    (
        data["common_input_bases"]
        - data["bases_mapped_cigar"]
    )
    / data["common_input_bases"]
    * 100.0
)


# ============================================================
# STEP 11: ERROR-ADJUSTED CIGAR YIELD
# ============================================================

# Equation:
#
# error-adjusted yield (%) =
#
# bases_mapped_cigar - mismatches
# ------------------------------- × 100
# common_input_bases

data["error_adjusted_cigar_yield_percent"] = (
    (
        data["bases_mapped_cigar"]
        - data["mismatches"]
    )
    / data["common_input_bases"]
    * 100.0
)


# ============================================================
# STEP 12: COMPARE CALCULATED MISMATCH RATE WITH ERROR_PERCENT
# ============================================================

data["error_percent_difference"] = (
    data["mismatch_error_calculated_percent"]
    - data["error_percent"]
)


data["absolute_error_percent_difference"] = (
    data["error_percent_difference"]
    .abs()
)


print()
print(
    "Maximum difference between calculated mismatch error "
    "and TSV error_percent:"
)

print(
    f"{data['absolute_error_percent_difference'].max():.6f} "
    "percentage points"
)


# ============================================================
# STEP 13: SELECT BENCHMARK ROWS
# ============================================================

plot_data = data.loc[
    data["sample"].isin(sample_order)
    & data["read_technology"].isin(technology_order)
    & data["aligner"].isin(aligner_order),
    [
        "sample",
        "read_technology",
        "aligner",

        "total_length",
        "common_input_bases",

        "bases_mapped_cigar",
        "mismatches",

        "cigar_mapped_percent_internal",
        "cigar_unaligned_percent_internal",

        "mismatch_total_percent",

        "mismatch_error_calculated_percent",
        "error_percent",
        "error_percent_difference",

        "input_normalized_cigar_yield_percent",
        "input_normalized_cigar_unaligned_percent",
        "error_adjusted_cigar_yield_percent",
    ],
].copy()


# ============================================================
# STEP 14: CHECK DUPLICATES
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
# STEP 15: CATEGORICAL ORDER
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
# STEP 16: SAVE NORMALIZED METRICS
# ============================================================

plot_data.to_csv(
    OUTPUT_SUMMARY,
    sep="\t",
    index=False,
)


print()
print("Normalized metrics:")
print(
    plot_data.to_string(
        index=False,
        float_format=lambda value: f"{value:.4f}",
    )
)


# ============================================================
# STEP 17: CREATE STACKED CIGAR FIGURE
# ============================================================

# For the 100% composition plot, retain the internal denominator
# because these three components must sum to exactly 100%.

plot_data["aligned_non_mismatch_percent"] = (
    (
        plot_data["bases_mapped_cigar"]
        - plot_data["mismatches"]
    )
    / plot_data["total_length"]
    * 100.0
)


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
    len(aligner_order)
)


for panel_index, technology in enumerate(
    technology_order
):

    axis = axes[
        panel_index
    ]


    technology_data = plot_data.loc[
        plot_data["read_technology"]
        == technology
    ]


    for sample in sample_order:

        sample_data = (
            technology_data.loc[
                technology_data["sample"]
                == sample
            ]
            .set_index("aligner")
            .reindex(aligner_order)
        )


        aligned_values = (
            sample_data[
                "aligned_non_mismatch_percent"
            ]
            .to_numpy(dtype=float)
        )


        mismatch_values = (
            sample_data[
                "mismatch_total_percent"
            ]
            .to_numpy(dtype=float)
        )


        unaligned_values = (
            sample_data[
                "cigar_unaligned_percent_internal"
            ]
            .to_numpy(dtype=float)
        )


        cigar_mapped_values = (
            sample_data[
                "cigar_mapped_percent_internal"
            ]
            .to_numpy(dtype=float)
        )


        positions = (
            x_positions
            + sample_offsets[sample]
        )


        base_color = (
            sample_colors[sample]
        )


        mismatch_color = lighten_color(
            base_color,
            amount=0.35,
        )


        unaligned_color = lighten_color(
            base_color,
            amount=0.72,
        )


        # ----------------------------------------------------
        # aligned non-mismatch
        # ----------------------------------------------------

        axis.bar(
            positions,
            aligned_values,
            width=bar_width,
            color=base_color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )


        # ----------------------------------------------------
        # mismatches
        # ----------------------------------------------------

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


        # ----------------------------------------------------
        # CIGAR unaligned
        # ----------------------------------------------------

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


        # ----------------------------------------------------
        # diagonal mismatch labels
        # ----------------------------------------------------

        for (
            mismatch_bar,
            mismatch_value,
        ) in zip(
            mismatch_bars,
            mismatch_values,
        ):

            if np.isnan(mismatch_value):
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


            if mismatch_value >= 0.20:

                axis.text(
                    x_center,
                    y_center,
                    f"{mismatch_value:.2f}%",
                    ha="center",
                    va="center",
                    fontsize=6.8,
                    color="black",
                    rotation=45,
                    rotation_mode="anchor",
                    clip_on=False,
                    zorder=7,
                )


        # ----------------------------------------------------
        # CIGAR mapped percentage above bars
        # ----------------------------------------------------

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

            if np.isnan(mapped_value):
                continue


            axis.text(
                position,
                100.20
                + mapped_label_offsets[sample],
                f"{mapped_value:.2f}%",
                ha="center",
                va="bottom",
                fontsize=7.4,
                color=base_color,
                clip_on=False,
                zorder=8,
            )


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


    axis.set_xticks(
        x_positions,
        aligner_order,
    )


    axis.set_xlabel(
        "Aligner",
        fontsize=13,
    )


    axis.set_ylim(
        88,
        102,
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


    axis.grid(
        axis="y",
        color="#D9D9D9",
        linewidth=0.7,
        alpha=0.75,
        zorder=0,
    )


    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    axis.tick_params(
        axis="both",
        labelsize=11,
    )


axes[0].set_ylabel(
    "CIGAR-aligned and unaligned bases (%)",
    fontsize=13,
)


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
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.99),
    fontsize=11,
)


component_legend_handles = [
    Patch(
        facecolor="#555555",
        label="Aligned, non-mismatch",
    ),
    Patch(
        facecolor="#999999",
        label="Mismatch bases",
    ),
    Patch(
        facecolor="#DDDDDD",
        label="CIGAR-unaligned bases",
    ),
]


figure.legend(
    handles=component_legend_handles,
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.905),
    fontsize=10.5,
)


figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.13,
    top=0.78,
    wspace=0.12,
)


figure.savefig(
    OUTPUT_STACKED_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)


figure.savefig(
    OUTPUT_STACKED_PDF,
    bbox_inches="tight",
    facecolor="white",
)


plt.close(figure)


# ============================================================
# STEP 18: SCATTER PLOT
#
# Input-normalized CIGAR yield vs mismatch error
# ============================================================

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(
        9.1,
        4.76,
    ),
    sharey=True,
)


for panel_index, technology in enumerate(
    technology_order
):

    axis = axes[
        panel_index
    ]


    technology_data = plot_data.loc[
        plot_data["read_technology"]
        == technology
    ]


    for aligner in aligner_order:

        aligner_data = technology_data.loc[
            technology_data["aligner"]
            == aligner
        ]


        for _, row in aligner_data.iterrows():

            sample = str(
                row["sample"]
            )


            axis.scatter(
                row[
                    "input_normalized_cigar_yield_percent"
                ],

                row[
                    "error_percent"
                ],

                s=95,

                marker=aligner_markers[
                    aligner
                ],

                color=sample_colors[
                    sample
                ],

                edgecolor="black",

                linewidth=0.7,

                alpha=0.90,

                zorder=4,
            )


    # --------------------------------------------------------
    # Correlation
    # --------------------------------------------------------

    valid_data = technology_data[
        [
            "input_normalized_cigar_yield_percent",
            "error_percent",
        ]
    ].dropna()


    if len(valid_data) >= 2:

        pearson_result = stats.pearsonr(
            valid_data["input_normalized_cigar_yield_percent"],
            valid_data["error_percent"],
        )

        pearson_r = pearson_result.statistic
        pearson_p = pearson_result.pvalue

        spearman_result = stats.spearmanr(
            valid_data["input_normalized_cigar_yield_percent"],
            valid_data["error_percent"],
        )

        spearman_rho = spearman_result.statistic
        spearman_p = spearman_result.pvalue

    else:

        pearson_r = np.nan
        pearson_p = np.nan
        spearman_rho = np.nan
        spearman_p = np.nan


    # Spearman is the primary statistic reported: the regression
    # residuals for this relationship are non-normal (Shapiro-Wilk,
    # see the supplementary correlation-diagnostics figures/table),
    # so the rank-based, outlier-robust Spearman correlation is the
    # statistically valid one here. Pearson is shown underneath only
    # as a diagnostic reference.

    def _stars(p_value):

        if np.isnan(p_value):
            return ""
        if p_value < 0.001:
            return "***"
        if p_value < 0.01:
            return "**"
        if p_value < 0.05:
            return "*"

        return "ns"


    axis.text(
        0.04,
        0.96,

        f"Spearman $\\rho$ = {spearman_rho:.2f} ({_stars(spearman_p)})\n"
        f"Pearson r = {pearson_r:.2f} ({_stars(pearson_p)})",

        transform=axis.transAxes,

        ha="left",
        va="top",

        fontsize=10,
    )


    if technology == "ONT":

        title = "ONT"
        panel_letter = "a"

    else:

        title = "PacBio HiFi"
        panel_letter = "b"


    axis.set_title(
        title,
        fontsize=16,
        pad=14,
    )


    axis.text(
        -0.08,
        1.03,
        panel_letter,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",
    )


    axis.grid(
        color="#D9D9D9",
        linewidth=0.7,
        alpha=0.75,
        zorder=0,
    )


    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


axes[0].set_ylabel(
    "Mismatch error (%)",
    fontsize=12,
)


figure.supxlabel(
    "Input-normalized CIGAR-aligned yield (%)",
    fontsize=12,
    y=0.10,
)


# ============================================================
# STEP 19: SCATTER LEGENDS
# ============================================================

sample_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor=sample_colors[sample],
        markeredgecolor="black",
        markersize=6,
        label=sample,
    )
    for sample in sample_order
]


aligner_handles = [
    Line2D(
        [0],
        [0],
        marker=aligner_markers[aligner],
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=6,
        label=aligner,
    )
    for aligner in aligner_order
]


legend_style = dict(
    frameon=False,
    fontsize=9,
    title_fontsize=9,
    columnspacing=0.9,
    handletextpad=0.35,
    borderaxespad=0.0,
    handlelength=1.2,
)


legend_samples = figure.legend(
    handles=sample_handles,
    title="GIAB sample",
    ncols=3,
    loc="upper right",
    bbox_to_anchor=(
        0.44,
        1.01,
    ),
    **legend_style,
)


figure.add_artist(
    legend_samples
)


figure.legend(
    handles=aligner_handles,
    title="Aligner",
    ncols=4,
    loc="upper left",
    bbox_to_anchor=(
        0.46,
        1.01,
    ),
    **legend_style,
)


figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.214,
    top=0.757,
    wspace=0.12,
)


figure.savefig(
    OUTPUT_SCATTER_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)


figure.savefig(
    OUTPUT_SCATTER_PDF,
    bbox_inches="tight",
    facecolor="white",
)


plt.close(figure)


# ============================================================
# STEP 20: SUMMARY BY ALIGNER / TECHNOLOGY
# ============================================================

summary = (
    plot_data
    .groupby(
        [
            "read_technology",
            "aligner",
        ],
        observed=True,
    )
    .agg(
        mean_input_normalized_cigar_yield=(
            "input_normalized_cigar_yield_percent",
            "mean",
        ),

        sd_input_normalized_cigar_yield=(
            "input_normalized_cigar_yield_percent",
            "std",
        ),

        mean_error_percent=(
            "error_percent",
            "mean",
        ),

        sd_error_percent=(
            "error_percent",
            "std",
        ),

        mean_error_adjusted_cigar_yield=(
            "error_adjusted_cigar_yield_percent",
            "mean",
        ),

        sd_error_adjusted_cigar_yield=(
            "error_adjusted_cigar_yield_percent",
            "std",
        ),
    )
    .reset_index()
)


print()
print("============================================================")
print("SUMMARY BY TECHNOLOGY / ALIGNER")
print("============================================================")

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    )
)


# ============================================================
# STEP 21: FINAL OUTPUTS
# ============================================================

print()
print("Analysis completed successfully.")

print()
print("Stacked CIGAR figure:")
print(OUTPUT_STACKED_PDF)

print()
print("Normalized yield vs error figure:")
print(OUTPUT_SCATTER_PDF)

print()
print("Normalized metrics TSV:")
print(OUTPUT_SUMMARY)