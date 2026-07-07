#!/usr/bin/env python3

from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator, ScalarFormatter


# ============================================================
# File paths
# ============================================================

INPUT_FILE = Path(
    "alignment_analysis/tables/alignment_summary.tsv"
)

OUTPUT_PREFIX = Path(
    "alignment_analysis/figures/"
    "mapped_bases_cigar_reference_style"
)


# ============================================================
# Expected order
# ============================================================

ALIGNMENT_ORDER = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]

SAMPLE_ORDER = [
    "HG002",
    "HG003",
    "HG004",
]


# ============================================================
# Colors similar to the reference graph
# ============================================================

COLORS = {
    "HG002_ont": "#F8766D",
    "HG002_pb": "#B79F00",
    "HG003_ont": "#00BA38",
    "HG003_pb": "#00BFC4",
    "HG004_ont": "#619CFF",
    "HG004_pb": "#F564E3",
}


# ============================================================
# Helper functions
# ============================================================

def clean_text(value):
    """
    Normalize spaces, underscores, slashes and capitalization.
    """

    text = str(value).strip().lower()

    text = re.sub(
        r"[\s_/]+",
        "-",
        text,
    )

    text = re.sub(
        r"-+",
        "-",
        text,
    )

    return text.strip("-")


def normalize_technology(value):
    """
    Convert different technology names into:
    - ont
    - pb
    """

    text = clean_text(value)

    if (
        "ont" in text
        or "nanopore" in text
    ):
        return "ont"

    if (
        text == "pb"
        or "pacbio" in text
        or "pac-bio" in text
        or "hifi" in text
        or "ccs" in text
    ):
        return "pb"

    return None


def normalize_alignment(value):
    """
    Convert different configuration labels into:
    - mm2-ont
    - mm2-pb
    - pbmm2-ont
    - pbmm2-pb
    """

    text = clean_text(value)

    # Detect mapper
    if "pbmm2" in text:
        mapper = "pbmm2"

    elif (
        "minimap2" in text
        or text.startswith("mm2")
    ):
        mapper = "mm2"

    else:
        mapper = None

    # Detect alignment preset
    if "ont" in text:
        preset = "ont"

    elif (
        "map-pb" in text
        or text.endswith("-pb")
        or "pacbio" in text
        or "hifi" in text
        or "ccs" in text
    ):
        preset = "pb"

    else:
        preset = None

    if mapper is None or preset is None:
        return None

    return f"{mapper}-{preset}"


# ============================================================
# Check input file
# ============================================================

if not INPUT_FILE.exists():
    sys.exit(
        "\nERROR: The input file was not found:\n"
        f"  {INPUT_FILE.resolve()}\n\n"
        "Check your available files with:\n"
        "  find alignment_analysis -maxdepth 3 -type f | sort\n"
    )


# ============================================================
# Read original table
# ============================================================

data = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)

data.columns = [
    str(column).strip().lower()
    for column in data.columns
]


# ============================================================
# Validate columns
# ============================================================

required_columns = {
    "sample",
    "read_technology",
    "configuration",
    "bases_mapped_cigar",
}

missing_columns = (
    required_columns
    - set(data.columns)
)

if missing_columns:
    sys.exit(
        "\nERROR: Missing required columns:\n  "
        + ", ".join(sorted(missing_columns))
        + "\n\nColumns found:\n  "
        + ", ".join(data.columns)
    )


data = data[
    [
        "sample",
        "read_technology",
        "configuration",
        "bases_mapped_cigar",
    ]
].copy()


# ============================================================
# Clean and normalize columns
# ============================================================

data["sample"] = (
    data["sample"]
    .astype(str)
    .str.strip()
    .str.upper()
)

data["technology"] = (
    data["read_technology"]
    .map(normalize_technology)
)

data["alignment"] = (
    data["configuration"]
    .map(normalize_alignment)
)

data["bases_mapped_cigar"] = pd.to_numeric(
    data["bases_mapped_cigar"]
    .astype(str)
    .str.replace(",", "", regex=False)
    .str.strip(),
    errors="coerce",
)


# ============================================================
# Detect invalid rows
# ============================================================

invalid_rows = data[
    data[
        [
            "technology",
            "alignment",
            "bases_mapped_cigar",
        ]
    ]
    .isna()
    .any(axis=1)
]

if not invalid_rows.empty:
    sys.exit(
        "\nERROR: Some rows could not be interpreted.\n"
        "Check read_technology, configuration and "
        "bases_mapped_cigar:\n\n"
        + invalid_rows.to_string(index=False)
    )


# ============================================================
# Detect duplicated runs
# ============================================================

duplicate_rows = data[
    data.duplicated(
        [
            "sample",
            "technology",
            "alignment",
        ],
        keep=False,
    )
]

if not duplicate_rows.empty:
    duplicate_rows = duplicate_rows.sort_values(
        [
            "sample",
            "technology",
            "alignment",
        ]
    )

    sys.exit(
        "\nERROR: Duplicated runs were found.\n"
        "There must be exactly one row for every:\n"
        "sample × technology × alignment\n\n"
        + duplicate_rows.to_string(index=False)
    )


# ============================================================
# Create original-file label
# ============================================================

data["orig_file"] = (
    data["sample"]
    + "_"
    + data["technology"]
)


# ============================================================
# Apply category order
# ============================================================

data["alignment"] = pd.Categorical(
    data["alignment"],
    categories=ALIGNMENT_ORDER,
    ordered=True,
)

data = data.sort_values(
    [
        "sample",
        "technology",
        "alignment",
    ]
).reset_index(drop=True)


# ============================================================
# Determine sample order
# ============================================================

samples_present = list(
    data["sample"].unique()
)

samples = [
    sample
    for sample in SAMPLE_ORDER
    if sample in samples_present
]

samples.extend(
    sample
    for sample in samples_present
    if sample not in samples
)

if not samples:
    sys.exit(
        "\nERROR: No valid samples were found.\n"
    )


# ============================================================
# Check missing combinations
# ============================================================

expected_combinations = pd.MultiIndex.from_product(
    [
        samples,
        ["ont", "pb"],
        ALIGNMENT_ORDER,
    ],
    names=[
        "sample",
        "technology",
        "alignment",
    ],
)

observed_combinations = pd.MultiIndex.from_frame(
    data[
        [
            "sample",
            "technology",
            "alignment",
        ]
    ].astype(str)
)

missing_combinations = (
    expected_combinations
    .difference(observed_combinations)
)


# ============================================================
# Print exact original values
# ============================================================

print()
print("=" * 95)
print("EXACT ORIGINAL VALUES USED IN THE GRAPH")
print("=" * 95)

print(
    data[
        [
            "sample",
            "technology",
            "alignment",
            "bases_mapped_cigar",
        ]
    ].to_string(index=False)
)

print()
print(
    "Minimum bases_mapped_cigar: "
    f"{data['bases_mapped_cigar'].min():,.0f}"
)

print(
    "Maximum bases_mapped_cigar: "
    f"{data['bases_mapped_cigar'].max():,.0f}"
)

if len(missing_combinations) == 0:
    print()
    print(
        "All expected sample × technology × alignment "
        "combinations are present."
    )

else:
    print()
    print("WARNING: Missing combinations:")

    for combination in missing_combinations:
        print(
            "  - "
            + " | ".join(combination)
        )

print("=" * 95)
print()


# ============================================================
# Create figure
# ============================================================

figure, axes = plt.subplots(
    1,
    len(samples),
    figsize=(
        4.25 * len(samples),
        4.7,
    ),
    sharex=True,
    sharey=False,
    squeeze=False,
)

axes = axes.ravel()

x_positions = np.arange(
    len(ALIGNMENT_ORDER)
)


# ============================================================
# Draw each sample panel
# ============================================================

for axis, sample in zip(
    axes,
    samples,
):

    sample_data = data[
        data["sample"] == sample
    ]

    for technology in [
        "ont",
        "pb",
    ]:

        line_data = sample_data[
            sample_data["technology"]
            == technology
        ]

        if line_data.empty:
            continue

        series = (
            line_data
            .set_index("alignment")[
                "bases_mapped_cigar"
            ]
            .reindex(ALIGNMENT_ORDER)
        )

        y_values = series.to_numpy(
            dtype=float
        )

        legend_label = (
            f"{sample}_{technology}"
        )

        # Dotted gray connection line
        axis.plot(
            x_positions,
            y_values,
            linestyle=(
                0,
                (
                    1.5,
                    2.5,
                ),
            ),
            linewidth=1.0,
            color="#BFBFBF",
            zorder=1,
        )

        # Colored observations
        axis.scatter(
            x_positions,
            y_values,
            s=32,
            color=COLORS.get(
                legend_label,
                "#444444",
            ),
            edgecolor="white",
            linewidth=0.55,
            zorder=2,
        )

    # X-axis
    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        ALIGNMENT_ORDER,
        fontsize=8,
    )

    # Y-axis tick count
    axis.yaxis.set_major_locator(
        MaxNLocator(
            nbins=5
        )
    )

    # Scientific notation:
    # tick labels such as 1.0, 1.5, 2.0
    # with ×10^9 shown beside the axis.
    scientific_formatter = ScalarFormatter(
        useMathText=True
    )

    scientific_formatter.set_scientific(
        True
    )

    scientific_formatter.set_powerlimits(
        (0, 0)
    )

    scientific_formatter.set_useOffset(
        False
    )

    axis.yaxis.set_major_formatter(
        scientific_formatter
    )

    axis.yaxis.get_offset_text().set_fontsize(
        8
    )

    axis.tick_params(
        axis="y",
        labelsize=8,
    )

    axis.grid(
        False
    )

    axis.margins(
        x=0.08,
        y=0.10,
    )

    # Black panel borders
    for spine in axis.spines.values():
        spine.set_linewidth(
            0.8
        )

        spine.set_color(
            "black"
        )

    # Facet header
    facet_strip = Rectangle(
        (
            0,
            1.0,
        ),
        1,
        0.075,
        transform=axis.transAxes,
        facecolor="#F7F7F7",
        edgecolor="black",
        linewidth=0.8,
        clip_on=False,
    )

    axis.add_patch(
        facet_strip
    )

    axis.text(
        0.5,
        1.037,
        sample,
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=9,
    )


# ============================================================
# Figure labels
# ============================================================

figure.suptitle(
    "Mapped bases across alignment configurations",
    x=0.055,
    y=0.985,
    ha="left",
    fontsize=12,
)

figure.supxlabel(
    "Alignment configuration",
    y=0.055,
    fontsize=10,
)

figure.supylabel(
    "Bases mapped from CIGAR",
    x=0.012,
    fontsize=10,
)


# ============================================================
# Legend
# ============================================================

available_labels = set(
    data["orig_file"]
)

legend_labels = [
    f"{sample}_{technology}"
    for sample in samples
    for technology in [
        "ont",
        "pb",
    ]
    if f"{sample}_{technology}"
    in available_labels
]

legend_handles = [
    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="none",
        markerfacecolor=COLORS.get(
            label,
            "#444444",
        ),
        markeredgecolor="white",
        markersize=6,
        label=label,
    )
    for label in legend_labels
]

figure.legend(
    handles=legend_handles,
    title="orig_file",
    loc="center left",
    bbox_to_anchor=(
        0.885,
        0.52,
    ),
    frameon=False,
    fontsize=8,
    title_fontsize=9,
)


# ============================================================
# Layout
# ============================================================

figure.subplots_adjust(
    left=0.075,
    right=0.875,
    bottom=0.16,
    top=0.865,
    wspace=0.24,
)


# ============================================================
# Save outputs
# ============================================================

OUTPUT_PREFIX.parent.mkdir(
    parents=True,
    exist_ok=True,
)

png_output = (
    OUTPUT_PREFIX
    .with_suffix(".png")
)

pdf_output = (
    OUTPUT_PREFIX
    .with_suffix(".pdf")
)

values_output = (
    OUTPUT_PREFIX
    .with_name(
        OUTPUT_PREFIX.name
        + "_values.tsv"
    )
)

figure.savefig(
    png_output,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

figure.savefig(
    pdf_output,
    bbox_inches="tight",
    facecolor="white",
)

plt.close(
    figure
)

data[
    [
        "sample",
        "technology",
        "alignment",
        "bases_mapped_cigar",
    ]
].to_csv(
    values_output,
    sep="\t",
    index=False,
)


# ============================================================
# Final report
# ============================================================

print("Created successfully:")

print(
    f"  PNG: {png_output.resolve()}"
)

print(
    f"  PDF: {pdf_output.resolve()}"
)

print(
    f"  Values: {values_output.resolve()}"
)
