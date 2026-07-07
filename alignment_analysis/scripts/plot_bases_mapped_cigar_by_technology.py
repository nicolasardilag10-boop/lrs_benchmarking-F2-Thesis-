from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MultipleLocator


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
input_path = Path(
    "alignment_analysis/tables/alignment_summary.tsv"
)

output_png = Path(
    "alignment_analysis/figures/"
    "bases_mapped_cigar_by_technology.png"
)

output_pdf = Path(
    "alignment_analysis/figures/"
    "bases_mapped_cigar_by_technology.pdf"
)


# ------------------------------------------------------------
# Read table
# ------------------------------------------------------------
data = pd.read_csv(
    input_path,
    sep="\t",
)


required_columns = [
    "sample",
    "read_technology",
    "configuration",
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
# Convert the metric to numeric
# No percentages and no normalization are performed here
# ------------------------------------------------------------
data["bases_mapped_cigar"] = pd.to_numeric(
    data["bases_mapped_cigar"],
    errors="raise",
)


# ------------------------------------------------------------
# Category order
# ------------------------------------------------------------
technology_order = [
    "ONT",
    "PacBio",
]

preferred_configuration_order = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]

configuration_order = [
    configuration
    for configuration in preferred_configuration_order
    if configuration in data["configuration"].unique()
]


unexpected_configurations = [
    configuration
    for configuration in data["configuration"].unique()
    if configuration not in configuration_order
]

configuration_order.extend(
    sorted(unexpected_configurations)
)


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
# Verification
# ------------------------------------------------------------
minimum_value = data["bases_mapped_cigar"].min()
maximum_value = data["bases_mapped_cigar"].max()


print("\nExact raw values used in the graph:\n")

print(
    data[
        [
            "sample",
            "read_technology",
            "configuration",
            "bases_mapped_cigar",
        ]
    ].to_string(
        index=False,
        formatters={
            "bases_mapped_cigar": lambda value: f"{value:,.0f}"
        },
    )
)


print("\nRaw range:")
print(f"{minimum_value:,.0f} to {maximum_value:,.0f} bases")

print("\nSame range in units of 10^7:")
print(
    f"{minimum_value / 1e7:.4f} to "
    f"{maximum_value / 1e7:.4f}"
)


# Sanity check for the requested graphical scale
if minimum_value < 1e7 or maximum_value > 2.5e7:
    raise ValueError(
        "Some values fall outside the requested plotting range "
        "of 1.0 × 10^7 to 2.5 × 10^7 bases."
    )


# ------------------------------------------------------------
# Create grouped boxplot
# ------------------------------------------------------------
fig, ax = plt.subplots(
    figsize=(11, 7),
)


technology_positions = np.array(
    [0, 1],
    dtype=float,
)


number_of_configurations = len(
    configuration_order
)


group_width = 0.68

configuration_spacing = (
    group_width
    / number_of_configurations
)

configuration_offsets = (
    np.arange(number_of_configurations)
    - (number_of_configurations - 1) / 2
) * configuration_spacing


colors = plt.rcParams[
    "axes.prop_cycle"
].by_key()["color"]


legend_handles = []


for configuration_index, configuration in enumerate(
    configuration_order
):
    color = colors[
        configuration_index % len(colors)
    ]

    legend_handles.append(
        Patch(
            facecolor=color,
            edgecolor="black",
            label=configuration,
            alpha=0.7,
        )
    )

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

        values = group[
            "bases_mapped_cigar"
        ].to_numpy()

        position = (
            technology_positions[technology_index]
            + configuration_offsets[configuration_index]
        )


        # Boxplot of HG002, HG003 and HG004
        boxplot = ax.boxplot(
            values,
            positions=[position],
            widths=configuration_spacing * 0.68,
            patch_artist=True,
            showfliers=False,
            medianprops={
                "color": "black",
                "linewidth": 1.5,
            },
            boxprops={
                "edgecolor": "black",
                "linewidth": 1.0,
            },
            whiskerprops={
                "color": "black",
                "linewidth": 1.0,
            },
            capprops={
                "color": "black",
                "linewidth": 1.0,
            },
        )

        boxplot["boxes"][0].set_facecolor(
            color
        )

        boxplot["boxes"][0].set_alpha(
            0.65
        )


        # Small fixed offsets so the three samples are visible
        point_offsets = np.linspace(
            -0.025,
            0.025,
            len(group),
        )

        ax.scatter(
            position + point_offsets,
            values,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            s=48,
            zorder=3,
        )


# ------------------------------------------------------------
# Y-axis
#
# The underlying values remain raw integers.
# Only the labels are displayed in units of 10^7.
# ------------------------------------------------------------
ax.set_ylim(
    1.0e7,
    2.5e7,
)

ax.yaxis.set_major_locator(
    MultipleLocator(
        0.25e7
    )
)

ax.yaxis.set_major_formatter(
    FuncFormatter(
        lambda value, position: f"{value / 1e7:.2f}".rstrip("0").rstrip(".")
    )
)


# ------------------------------------------------------------
# X-axis
# ------------------------------------------------------------
ax.set_xticks(
    technology_positions
)

ax.set_xticklabels(
    technology_order,
    fontsize=12,
)


# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------
ax.set_xlabel(
    "Read technology",
    fontsize=12,
)

ax.set_ylabel(
    "CIGAR-mapped bases (×10⁷ bases)",
    fontsize=12,
)

ax.set_title(
    "CIGAR-mapped bases by sequencing technology",
    fontsize=15,
    fontweight="bold",
    pad=16,
)


# ------------------------------------------------------------
# Appearance
# ------------------------------------------------------------
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    axis="both",
    labelsize=11,
)

ax.legend(
    handles=legend_handles,
    title="Alignment configuration",
    frameon=False,
    loc="upper right",
)


fig.text(
    0.5,
    0.015,
    (
        "Boxes and points represent HG002, HG003 and HG004. "
        "The y-axis shows raw CIGAR-mapped bases in units of 10⁷."
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
