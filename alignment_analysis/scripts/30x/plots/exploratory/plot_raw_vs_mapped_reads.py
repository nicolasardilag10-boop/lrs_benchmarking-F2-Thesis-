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
    / "final"
    / "alignment_benchmark_30x.tsv"
)

OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "raw_vs_mapped_reads_30x.png"
)

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "raw_vs_mapped_reads_30x.pdf"
)

OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,)


# ============================================================
# STEP 2: COLORS AND PLOTTING ORDER
# ============================================================

# Same GIAB colors used in the other benchmarking figures.
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

# "Raw reads" is not an aligner: it is the shared sequencing input
# that all four aligners were given, drawn once per sample before
# any tool-specific bar.
category_order = [
    "Raw reads",
    *aligner_order,
]

technology_order = [
    "ONT",
    "PacBio",
]


# ============================================================
# STEP 3: LOAD DATA
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
# STEP 4: VALIDATE REQUIRED COLUMNS
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "raw_total_sequences",
    "reads_mapped",
}

missing_columns = required_columns.difference(
    data.columns
)

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        f"{sorted(missing_columns)}"
    )


data["raw_total_sequences"] = pd.to_numeric(
    data["raw_total_sequences"],
    errors="raise",
)

data["reads_mapped"] = pd.to_numeric(
    data["reads_mapped"],
    errors="raise",
)


# ============================================================
# STEP 5: SELECT RELEVANT DATA
# ============================================================

plot_data = data.loc[
    data["sample"].isin(sample_order)
    & data["read_technology"].isin(technology_order)
    & data["aligner"].isin(aligner_order),
    [
        "sample",
        "read_technology",
        "aligner",
        "raw_total_sequences",
        "reads_mapped",
    ],
].copy()


# ============================================================
# STEP 6: CHECK FOR DUPLICATED OBSERVATIONS
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
        ].to_string(
            index=False
        )
    )


# ============================================================
# STEP 7: RESOLVE A SINGLE SHARED "RAW READS" VALUE
# ============================================================
#
# raw_total_sequences is reported per aligner because it comes from
# each aligner's own samtools stats output. minimap2, pbmm2 and VG
# Giraffe always agree with each other (they all report every input
# read, mapped or not). VACmap's own stats file does not contain
# unmapped-read records, so its raw_total_sequences under-counts the
# true sequencing input by exactly its unmapped-read count.
#
# Using each aligner's own value as "raw reads" would therefore make
# VACmap look like it started from fewer reads than the others, which
# is not true. Instead, the true shared raw-read count is taken from
# the aligners whose raw_total_sequences agree with each other, and
# that single value is used for every aligner's "Raw reads" bar.

reference_aligners = [
    "minimap2",
    "pbmm2",
    "VG Giraffe",
]

raw_reads_reference = plot_data.loc[
    plot_data["aligner"].isin(reference_aligners)
]

raw_reads_by_sample = raw_reads_reference.groupby(
    [
        "sample",
        "read_technology",
    ]
)["raw_total_sequences"].nunique()

if (raw_reads_by_sample > 1).any():
    disagreeing = raw_reads_by_sample.loc[
        raw_reads_by_sample > 1
    ]
    raise ValueError(
        "minimap2, pbmm2 and VG Giraffe disagree on raw_total_sequences "
        f"for:\n{disagreeing.to_string()}"
    )

raw_reads = (
    raw_reads_reference.groupby(
        [
            "sample",
            "read_technology",
        ]
    )["raw_total_sequences"]
    .first()
)


# ============================================================
# STEP 8: CREATE FIGURE
# ============================================================

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(18.0, 10.2),
)

x_positions = np.arange(
    len(category_order)
)


# ============================================================
# STEP 9: DRAW STACKED RAW / MAPPED BARS
# ============================================================

for panel_index, technology in enumerate(
    technology_order
):

    axis = axes[panel_index]

    technology_data = plot_data.loc[
        plot_data["read_technology"]
        == technology
    ]

    bottoms = np.zeros(
        len(category_order)
    )

    for sample in sample_order:

        base_color = sample_colors[sample]

        sample_data = (
            technology_data.loc[
                technology_data["sample"]
                == sample
            ]
            .set_index("aligner")
            .reindex(aligner_order)
        )

        mapped_values = (
            sample_data["reads_mapped"]
            .to_numpy(dtype=float)
        )

        raw_value = raw_reads.loc[
            (sample, technology)
        ]

        # First category ("Raw reads") + one value per aligner.
        values = np.concatenate(
            [
                [raw_value],
                mapped_values,
            ]
        ) / 1e6

        axis.bar(
            x_positions,
            values,
            width=0.6,
            bottom=bottoms,
            color=base_color,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )

        bottoms += values


    # ========================================================
    # PANEL TITLES
    # ========================================================

    axis.set_title(
        (
            "ONT"
            if technology == "ONT"
            else "PacBio HiFi"
        ),
        fontsize=21,
        pad=16,
    )

    axis.text(
        -0.08,
        1.035,
        (
            "a"
            if technology == "ONT"
            else "b"
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=23,
        fontweight="bold",
    )


    # ========================================================
    # X AXIS
    # ========================================================

    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        ALIGNER_ORDER,
        rotation=45,
        ha="right",
        fontsize=10,
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
        labelsize=15,
        colors="black",
    )

    axis.set_facecolor("white")


# ============================================================
# STEP 10: Y AXIS LABEL
# ============================================================

axes[0].set_ylabel(
    "Number of reads (millions)",
    fontsize=17,
)


# ============================================================
# STEP 11: GIAB SAMPLE LEGEND
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
    title_fontsize=20,
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.98,
    ),
    columnspacing=1.8,
    handlelength=1.8,
    fontsize=15,
)


# ============================================================
# STEP 12: FINAL LAYOUT
# ============================================================

figure.patch.set_facecolor("white")

figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.22,
    top=0.82,
    wspace=0.12,
)


# ============================================================
# STEP 13: SAVE FIGURE
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
# STEP 14: CONFIRM SUCCESS
# ============================================================

print()
print("Raw vs mapped reads figure created successfully.")
print()
print("PNG:")
print(OUTPUT_PNG)
print()
print("PDF:")
print(OUTPUT_PDF)
