#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scipy import stats

sys.path.insert(0, str(next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
) / "alignment_analysis" / "scripts"))
from utils.plot_style import detect_font_family, mirror_to_report

# Match the font used across the main (final/) figures so this
# supplementary diagnostic figure is typographically consistent
# with the rest of the manuscript.
plt.rcParams.update({
    "font.family": detect_font_family(),
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# STEP 1: PROJECT PATHS
# ============================================================

# Script location:
#
# lrs_benchmarking/
# └── alignment_analysis/
#     └── scripts/
#         └── 30x/
#             └── Normalization_homogenization_graphs.py
#

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())


INPUT_TSV = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv"
)


OUTPUT_DIR = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "correlation_diagnostics"
)


OUTPUT_TABLE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/cigar_yield_error_correlation_diagnostics.tsv"
)


OUTPUT_NORMALIZED_DATA = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/cigar_yield_error_normalized_data.tsv"
)


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_TABLE.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# STEP 2: ORDERS AND COLORS
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


technology_labels = {
    "ONT": "ONT",
    "PacBio": "PacBio HiFi",
}


sample_colors = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}


aligner_markers = {
    "minimap2": "o",
    "pbmm2": "s",
    "VACmap": "^",
    "VG Giraffe": "D",
}


# ============================================================
# STEP 3: LOAD ORIGINAL ALIGNMENT TABLE
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
# STEP 4: VALIDATE REQUIRED COLUMNS
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
# STEP 5: NUMERIC CONVERSION
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


# ============================================================
# STEP 6: KEEP ONLY BENCHMARK CONDITIONS
# ============================================================

plot_data = data.loc[
    data["sample"].isin(sample_order)
    & data["read_technology"].isin(technology_order)
    & data["aligner"].isin(aligner_order)
].copy()


# ============================================================
# STEP 7: COMMON INPUT DENOMINATOR
# ============================================================

# For every HG sample × sequencing technology:
#
# use the same denominator for all four aligners.
#
# In the current summary, minimap2/pbmm2/VG generally retain
# the original total_length while VACmap reports a reduced
# total_length.
#
# Therefore we use the maximum available total_length within
# each sample × technology as the common denominator.
#
# Ideally this should later be replaced with the exact number
# of bases calculated directly from the original FASTQ.

common_input = (
    plot_data
    .groupby(
        [
            "sample",
            "read_technology",
        ],
        observed=True,
    )["total_length"]
    .max()
    .rename(
        "common_input_bases"
    )
    .reset_index()
)


plot_data = plot_data.merge(
    common_input,
    on=[
        "sample",
        "read_technology",
    ],
    how="left",
)


print()
print("Common input denominators:")
print(
    common_input.to_string(
        index=False
    )
)


# ============================================================
# STEP 8: NORMALIZED CIGAR YIELD
# ============================================================

# Input-normalized CIGAR-aligned yield (%)
#
#              CIGAR-aligned bases
# 100 × --------------------------------
#              common input bases

plot_data[
    "input_normalized_cigar_yield_percent"
] = (
    plot_data["bases_mapped_cigar"]
    / plot_data["common_input_bases"]
    * 100.0
)


# ============================================================
# STEP 9: CALCULATED MISMATCH ERROR
# ============================================================

plot_data[
    "calculated_mismatch_error_percent"
] = (
    plot_data["mismatches"]
    / plot_data["bases_mapped_cigar"]
    * 100.0
)


# Compare against existing error_percent.
plot_data[
    "error_percent_difference"
] = (
    plot_data[
        "calculated_mismatch_error_percent"
    ]
    - plot_data[
        "error_percent"
    ]
)


print()
print(
    "Maximum difference between calculated mismatch error "
    "and existing error_percent:"
)

print(
    f"{plot_data['error_percent_difference'].abs().max():.6f}"
)


# ============================================================
# STEP 10: ERROR-ADJUSTED CIGAR YIELD
# ============================================================

plot_data[
    "error_adjusted_cigar_yield_percent"
] = (
    (
        plot_data["bases_mapped_cigar"]
        - plot_data["mismatches"]
    )
    / plot_data["common_input_bases"]
    * 100.0
)


# ============================================================
# STEP 11: SAVE NORMALIZED VALUES
# ============================================================

normalized_columns = [
    "sample",
    "read_technology",
    "aligner",
    "total_length",
    "common_input_bases",
    "bases_mapped_cigar",
    "mismatches",
    "input_normalized_cigar_yield_percent",
    "error_percent",
    "calculated_mismatch_error_percent",
    "error_percent_difference",
    "error_adjusted_cigar_yield_percent",
]


plot_data[
    normalized_columns
].to_csv(
    OUTPUT_NORMALIZED_DATA,
    sep="\t",
    index=False,
)


print()
print("Normalized data saved:")
print(OUTPUT_NORMALIZED_DATA)


# ============================================================
# STEP 12: STORAGE FOR STATISTICAL RESULTS
# ============================================================

diagnostic_results = []
technology_stats = {}


# ============================================================
# STEP 13: TECHNOLOGY-SPECIFIC ANALYSIS
# ============================================================

# Only the statistics are computed in this loop. The combined
# figure (STEP 19) is built afterwards so that ONT and PacBio
# HiFi can be drawn side by side instead of as two separate
# full-page figures.

for technology in technology_order:

    technology_data = (
        plot_data.loc[
            plot_data[
                "read_technology"
            ]
            == technology,
            [
                "sample",
                "aligner",
                "input_normalized_cigar_yield_percent",
                "error_percent",
            ],
        ]
        .dropna()
        .copy()
        .reset_index(
            drop=True
        )
    )


    x = (
        technology_data[
            "input_normalized_cigar_yield_percent"
        ]
        .to_numpy(
            dtype=float
        )
    )


    y = (
        technology_data[
            "error_percent"
        ]
        .to_numpy(
            dtype=float
        )
    )


    n = len(x)


    if n < 4:

        print()
        print(
            f"Not enough observations for {technology}."
        )

        continue


    # ========================================================
    # STEP 14: LINEAR REGRESSION
    # ========================================================

    regression = stats.linregress(
        x,
        y,
    )


    slope = (
        regression.slope
    )

    intercept = (
        regression.intercept
    )


    fitted = (
        intercept
        + slope * x
    )


    residuals = (
        y
        - fitted
    )


    r_squared = (
        regression.rvalue ** 2
    )


    # ========================================================
    # STEP 15: PEARSON
    # ========================================================

    pearson_result = stats.pearsonr(
        x,
        y,
    )


    pearson_r = (
        pearson_result.statistic
    )


    pearson_p = (
        pearson_result.pvalue
    )


    # ========================================================
    # STEP 16: SPEARMAN
    # ========================================================

    spearman_result = stats.spearmanr(
        x,
        y,
    )


    spearman_rho = (
        spearman_result.statistic
    )


    spearman_p = (
        spearman_result.pvalue
    )


    # ========================================================
    # STEP 17: SHAPIRO-WILK OF REGRESSION RESIDUALS
    # ========================================================

    shapiro_result = stats.shapiro(
        residuals
    )


    shapiro_W = (
        shapiro_result.statistic
    )


    shapiro_p = (
        shapiro_result.pvalue
    )


    # ========================================================
    # STEP 18: STORE STATISTICS
    # ========================================================

    diagnostic_results.append(
        {
            "technology": technology,

            "n": n,

            "pearson_r": pearson_r,
            "pearson_p": pearson_p,

            "spearman_rho": spearman_rho,
            "spearman_p": spearman_p,

            "shapiro_residual_W": shapiro_W,
            "shapiro_residual_p": shapiro_p,

            "r_squared": r_squared,

            "regression_slope": slope,

            "pearson_spearman_difference": abs(
                pearson_r
                - spearman_rho
            ),
        }
    )


    (
        theoretical_quantiles,
        ordered_residuals,
    ), (
        qq_slope,
        qq_intercept,
        qq_r,
    ) = stats.probplot(
        residuals,
        dist="norm",
    )


    technology_stats[technology] = {
        "technology_data": technology_data,
        "x": x,
        "y": y,
        "n": n,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "shapiro_W": shapiro_W,
        "shapiro_p": shapiro_p,
        "theoretical_quantiles": theoretical_quantiles,
        "ordered_residuals": ordered_residuals,
        "qq_slope": qq_slope,
        "qq_intercept": qq_intercept,
    }


# ============================================================
# STEP 19: COMBINED SIX-PANEL FIGURE
# ============================================================

# Layout: 3 rows (graph type) x 2 columns (technology), so
# ONT and PacBio HiFi can be compared directly within each
# row instead of scanning two separate full-page figures.
#
# Column 0 = ONT, column 1 = PacBio HiFi.
# Row 0 = yield-error relationship
# Row 1 = residual Q-Q diagnostic
# Row 2 = correlation consistency
#
# Panel letters follow row-major reading order:
# a = ONT row 0, b = PacBio row 0,
# c = ONT row 1, d = PacBio row 1,
# e = ONT row 2, f = PacBio row 2.

panel_letters = ["a", "b", "c", "d", "e", "f"]

row_titles = [
    "Yield–error relationship",
    "Residual Q-Q diagnostic",
    "Correlation consistency",
]

figure, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(
        11.5,
        13.2,
    ),
)


for col, technology in enumerate(technology_order):

    stats_for_technology = technology_stats[technology]
    technology_label = technology_labels[technology]

    x = stats_for_technology["x"]
    slope = stats_for_technology["slope"]
    intercept = stats_for_technology["intercept"]

    # --------------------------------------------------------
    # ROW 0: NORMALIZED CIGAR YIELD VS ERROR
    # --------------------------------------------------------

    axis = axes[0, col]

    for _, row in stats_for_technology["technology_data"].iterrows():

        sample = str(row["sample"])
        aligner = str(row["aligner"])

        axis.scatter(
            row["input_normalized_cigar_yield_percent"],
            row["error_percent"],
            s=100,
            marker=aligner_markers[aligner],
            color=sample_colors[sample],
            edgecolor="black",
            linewidth=0.7,
            alpha=0.95,
            zorder=4,
        )

    x_line = np.linspace(x.min(), x.max(), 200)
    y_line = intercept + slope * x_line

    axis.plot(
        x_line,
        y_line,
        color="black",
        linewidth=1.4,
        zorder=3,
    )

    axis.set_xlabel(
        "Input-normalized CIGAR-aligned yield (%)",
        fontsize=11,
    )

    axis.set_ylabel(
        "Mismatch error (%)",
        fontsize=11,
    )

    axis.set_title(
        row_titles[0],
        fontsize=12,
        fontweight="bold",
    )

    axis.text(
        0.04,
        0.96,
        (
            f"$R^2$ = {stats_for_technology['r_squared']:.2f}\n"
            f"n = {stats_for_technology['n']}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    # --------------------------------------------------------
    # ROW 1: RESIDUAL Q-Q PLOT
    # --------------------------------------------------------

    axis = axes[1, col]

    theoretical_quantiles = stats_for_technology["theoretical_quantiles"]
    ordered_residuals = stats_for_technology["ordered_residuals"]
    qq_slope = stats_for_technology["qq_slope"]
    qq_intercept = stats_for_technology["qq_intercept"]

    axis.scatter(
        theoretical_quantiles,
        ordered_residuals,
        s=55,
        facecolor="white",
        edgecolor="black",
        linewidth=0.9,
        zorder=4,
    )

    qq_x = np.asarray(
        [
            theoretical_quantiles.min(),
            theoretical_quantiles.max(),
        ]
    )

    qq_y = qq_intercept + qq_slope * qq_x

    axis.plot(
        qq_x,
        qq_y,
        color="black",
        linewidth=1.2,
        zorder=3,
    )

    axis.set_xlabel(
        "Theoretical normal quantiles",
        fontsize=11,
    )

    axis.set_ylabel(
        "Ordered residuals",
        fontsize=11,
    )

    axis.set_title(
        row_titles[1],
        fontsize=12,
        fontweight="bold",
    )

    axis.text(
        0.04,
        0.96,
        (
            "Shapiro-Wilk\n"
            f"W = {stats_for_technology['shapiro_W']:.3f}\n"
            f"p = {stats_for_technology['shapiro_p']:.3f}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    # --------------------------------------------------------
    # ROW 2: PEARSON VS SPEARMAN
    # --------------------------------------------------------

    axis = axes[2, col]

    correlation_values = [
        stats_for_technology["pearson_r"],
        stats_for_technology["spearman_rho"],
    ]

    correlation_labels = [
        "Pearson r",
        "Spearman ρ",
    ]

    y_positions = np.array([1, 0])

    bars = axis.barh(
        y_positions,
        correlation_values,
        height=0.52,
        color="white",
        edgecolor="black",
        linewidth=0.9,
        zorder=3,
    )

    axis.axvline(
        0,
        color="black",
        linewidth=1.0,
        zorder=2,
    )

    for bar, value in zip(bars, correlation_values):

        if value >= 0:
            x_label = value + 0.04
            ha = "left"
        else:
            x_label = value - 0.04
            ha = "right"

        axis.text(
            x_label,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}",
            ha=ha,
            va="center",
            fontsize=11,
            fontweight="bold",
            color="black",
        )

    axis.set_yticks(y_positions, correlation_labels)
    axis.set_xlim(-1.15, 1.15)
    axis.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])

    axis.set_xlabel(
        "Correlation coefficient",
        fontsize=11,
    )

    axis.set_title(
        row_titles[2],
        fontsize=12,
        fontweight="bold",
    )

    axis.text(
        0.04,
        0.08,
        (
            f"Pearson p = {stats_for_technology['pearson_p']:.3g}\n"
            f"Spearman p = {stats_for_technology['spearman_p']:.3g}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )


# ============================================================
# PANEL LETTERS (row-major: a-b top row, c-d middle, e-f bottom)
# ============================================================

for row in range(3):

    for col in range(2):

        axis = axes[row, col]

        axis.text(
            -0.10,
            1.10,
            panel_letters[row * 2 + col],
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )


# ============================================================
# COMMON AXIS STYLE
# ============================================================

for axis in axes.flat:

    axis.grid(
        alpha=0.20,
        linewidth=0.6,
        zorder=0,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=10)


# ============================================================
# SHARED LEGENDS
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
    fontsize=7.5,
    title_fontsize=7.5,
    columnspacing=0.8,
    handletextpad=0.3,
    borderaxespad=0.0,
    handlelength=1.1,
)

sample_legend = figure.legend(
    handles=sample_handles,
    title="GIAB sample",
    ncols=3,
    loc="upper right",
    bbox_to_anchor=(0.47, 0.975),
    **legend_style,
)

figure.add_artist(sample_legend)

figure.legend(
    handles=aligner_handles,
    title="Aligner",
    ncols=4,
    loc="upper left",
    bbox_to_anchor=(0.49, 0.975),
    **legend_style,
)


# ============================================================
# TITLE AND LAYOUT
# ============================================================

figure.subplots_adjust(
    left=0.07,
    right=0.98,
    bottom=0.045,
    top=0.85,
    wspace=0.30,
    hspace=0.42,
)


# ============================================================
# COLUMN HEADERS (technology, above row 0)
# ============================================================

figure.canvas.draw()

column_centers = [
    (
        axes[0, col].get_position().x0
        + axes[0, col].get_position().x1
    )
    / 2
    for col in range(2)
]

for col, technology in enumerate(technology_order):

    figure.text(
        column_centers[col],
        0.905,
        technology_labels[technology],
        ha="center",
        va="bottom",
        fontsize=17,
        fontweight="bold",
    )


# ============================================================
# SAVE FIGURE
# ============================================================

output_png = (
    OUTPUT_DIR
    / "cigar_yield_error_diagnostics_combined.png"
)

output_pdf = (
    OUTPUT_DIR
    / "cigar_yield_error_diagnostics_combined.pdf"
)

figure.savefig(
    output_png,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

figure.savefig(
    output_pdf,
    bbox_inches="tight",
    facecolor="white",
)

plt.close(figure)
mirror_to_report(output_pdf)

print()
print("Combined ONT / PacBio HiFi diagnostic figure:")
print(output_pdf)


# ============================================================
# STEP 20: SAVE STATISTICAL SUMMARY
# ============================================================

results_table = pd.DataFrame(
    diagnostic_results
)


results_table.to_csv(
    OUTPUT_TABLE,

    sep="\t",

    index=False,
)

mirror_to_report(OUTPUT_TABLE)


print()
print(
    "============================================================"
)

print(
    "CORRELATION DIAGNOSTIC SUMMARY"
)

print(
    "============================================================"
)


print(
    results_table.to_string(
        index=False,

        float_format=lambda value: f"{value:.4f}",
    )
)


print()
print(
    "Statistical summary:"
)

print(
    OUTPUT_TABLE
)


print()
print(
    "Analysis completed successfully."
)
