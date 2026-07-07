from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


# ------------------------------------------------------------
# Input and output files
# ------------------------------------------------------------
input_path = Path(
    "alignment_analysis/tables/alignment_summary.tsv"
)

output_png = Path(
    "alignment_analysis/figures/"
    "cigar_mapped_percent_by_technology.png"
)

output_pdf = Path(
    "alignment_analysis/figures/"
    "cigar_mapped_percent_by_technology.pdf"
)


# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
data = pd.read_csv(input_path, sep="\t")


required_columns = [
    "sample",
    "read_technology",
    "configuration",
    "total_length",
    "bases_mapped_cigar",
]

missing_columns = [
    column
    for column in required_columns
    if column not in data.columns
]

if missing_columns:
    raise ValueError(
        "Missing required columns: "
        + ", ".join(missing_columns)
    )


# ------------------------------------------------------------
# Convert numeric columns
# ------------------------------------------------------------
data["total_length"] = pd.to_numeric(
    data["total_length"],
    errors="raise",
)

data["bases_mapped_cigar"] = pd.to_numeric(
    data["bases_mapped_cigar"],
    errors="raise",
)


# ------------------------------------------------------------
# Normalize CIGAR-mapped bases
# ------------------------------------------------------------
data["cigar_mapped_percent"] = (
    data["bases_mapped_cigar"]
    / data["total_length"]
    * 100
)


# ------------------------------------------------------------
# Order of categories
# ------------------------------------------------------------
technology_order = [
    "ONT",
    "PacBio",
]

configuration_order = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]

data["read_technology"] = pd.Categorical(
    data["read_technology"],
    categories=technology_order,
    ordered=True,
)

data["configuration"] = pd.Categorical(
    data["configuration"],
    categories=configuration_order,
    ordered=True,
)


data = data.sort_values(
    [
        "read_technology",
        "configuration",
        "sample",
    ]
)


# ------------------------------------------------------------
# Print values used
# ------------------------------------------------------------
print("\nIndividual values plotted:\n")

print(
    data[
        [
            "sample",
            "read_technology",
            "configuration",
            "bases_mapped_cigar",
            "total_length",
            "cigar_mapped_percent",
        ]
    ].to_string(
        index=False,
        formatters={
            "cigar_mapped_percent": lambda value: f"{value:.4f}%"
        },
    )
)


summary = (
    data.groupby(
        [
            "read_technology",
            "configuration",
        ],
        observed=True,
    )["cigar_mapped_percent"]
    .agg(
        mean="mean",
        minimum="min",
        maximum="max",
    )
    .reset_index()
)


print("\nMean values by technology and configuration:\n")

print(
    summary.to_string(
        index=False,
        formatters={
            "mean": lambda value: f"{value:.4f}%",
            "minimum": lambda value: f"{value:.4f}%",
            "maximum": lambda value: f"{value:.4f}%",
        },
    )
)


# ------------------------------------------------------------
# Create figure
# ------------------------------------------------------------
fig, ax = plt.subplots(
    figsize=(10, 7)
)


technology_positions = np.arange(
    len(technology_order)
)

configuration_offsets = np.linspace(
    -0.27,
    0.27,
    len(configuration_order),
)


default_colors = plt.rcParams[
    "axes.prop_cycle"
].by_key()["color"]


# ------------------------------------------------------------
# Plot each alignment configuration
# ------------------------------------------------------------
for configuration_index, configuration in enumerate(
    configuration_order
):
    color = default_colors[
        configuration_index % len(default_colors)
    ]

    mean_x_values = []
    mean_y_values = []

    for technology_index, technology in enumerate(
        technology_order
    ):
        group = data[
            (
                data["read_technology"] == technology
            )
            &
            (
                data["configuration"] == configuration
            )
        ].copy()

        if group.empty:
            continue

        base_x = (
            technology_positions[technology_index]
            + configuration_offsets[configuration_index]
        )

        # Small deterministic horizontal separation between samples
        sample_offsets = np.linspace(
            -0.025,
            0.025,
            len(group),
        )

        sample_x = base_x + sample_offsets

        # Individual sample values
        ax.scatter(
            sample_x,
            group["cigar_mapped_percent"],
            s=55,
            alpha=0.75,
            color=color,
            zorder=3,
        )

        mean_value = group[
            "cigar_mapped_percent"
        ].mean()

        minimum_value = group[
            "cigar_mapped_percent"
        ].min()

        maximum_value = group[
            "cigar_mapped_percent"
        ].max()

        # Range across HG002, HG003 and HG004
        ax.vlines(
            base_x,
            minimum_value,
            maximum_value,
            color=color,
            linewidth=2,
            alpha=0.8,
            zorder=2,
        )

        # Mean value
        ax.scatter(
            base_x,
            mean_value,
            marker="D",
            s=95,
            color=color,
            edgecolor="black",
            linewidth=0.7,
            zorder=4,
        )

        mean_x_values.append(base_x)
        mean_y_values.append(mean_value)

    # Connect ONT and PacBio means for the same configuration
    if len(mean_x_values) == 2:
        ax.plot(
            mean_x_values,
            mean_y_values,
            color=color,
            linewidth=1.8,
            alpha=0.75,
            label=configuration,
            zorder=1,
        )


# ------------------------------------------------------------
# Reference line at 100%
# ------------------------------------------------------------
ax.axhline(
    100,
    linestyle="--",
    linewidth=1.2,
    alpha=0.65,
)


# ------------------------------------------------------------
# Axis limits
# ------------------------------------------------------------
minimum_plot_value = data[
    "cigar_mapped_percent"
].min()

maximum_plot_value = data[
    "cigar_mapped_percent"
].max()

ax.set_ylim(
    minimum_plot_value - 0.7,
    maximum_plot_value + 0.7,
)


# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------
ax.set_xticks(
    technology_positions
)

ax.set_xticklabels(
    technology_order,
    fontsize=12,
)

ax.set_xlabel(
    "Read technology",
    fontsize=12,
)

ax.set_ylabel(
    "CIGAR-mapped bases as percentage of total read bases",
    fontsize=12,
)

ax.set_title(
    "CIGAR-mapped bases compared by sequencing technology",
    fontsize=15,
    fontweight="bold",
    pad=16,
)


ax.yaxis.set_major_formatter(
    PercentFormatter(
        xmax=100,
        decimals=1,
    )
)


# ------------------------------------------------------------
# Appearance
# ------------------------------------------------------------
ax.grid(
    axis="y",
    linestyle="--",
    alpha=0.25,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


ax.legend(
    title="Alignment configuration",
    frameon=False,
    loc="lower right",
)


fig.text(
    0.5,
    0.015,
    (
        "Small circles represent HG002, HG003 and HG004. "
        "Diamonds represent the mean. Vertical lines show the sample range."
    ),
    ha="center",
    fontsize=9,
)


fig.tight_layout(
    rect=[0, 0.05, 1, 1]
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
output_png.parent.mkdir(
    parents=True,
    exist_ok=True,
)

fig.savefig(
    output_png,
    dpi=300,
    bbox_inches="tight",
)

fig.savefig(
    output_pdf,
    bbox_inches="tight",
)

plt.close(fig)


print("\nFigure created:")
print(output_png.resolve())
print(output_pdf.resolve())
