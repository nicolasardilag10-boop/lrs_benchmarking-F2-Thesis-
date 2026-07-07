from pathlib import Path
import math

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import PercentFormatter


# ------------------------------------------------------------
# File locations
# ------------------------------------------------------------
input_path = Path(
    "alignment_analysis/tables/alignment_summary.tsv"
)

output_png = Path(
    "alignment_analysis/figures/"
    "cigar_mapped_percent_diagnostic.png"
)

output_pdf = Path(
    "alignment_analysis/figures/"
    "cigar_mapped_percent_diagnostic.pdf"
)


# ------------------------------------------------------------
# Read and validate the table
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
# Calculate CIGAR-mapped bases as percentage of total bases
# ------------------------------------------------------------
data["cigar_mapped_percent"] = (
    data["bases_mapped_cigar"]
    / data["total_length"]
    * 100
)


# ------------------------------------------------------------
# Configuration order
# Keep the names exactly as they currently appear in the table
# ------------------------------------------------------------
configuration_order = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]

data["configuration"] = pd.Categorical(
    data["configuration"],
    categories=configuration_order,
    ordered=True,
)

data = data.sort_values(
    [
        "sample",
        "read_technology",
        "configuration",
    ]
)


# ------------------------------------------------------------
# Print the values used in the graph
# ------------------------------------------------------------
print("\nValues plotted:\n")

print(
    data[
        [
            "sample",
            "read_technology",
            "configuration",
            "total_length",
            "bases_mapped_cigar",
            "cigar_mapped_percent",
        ]
    ].to_string(
        index=False,
        formatters={
            "cigar_mapped_percent": lambda value: f"{value:.4f}%"
        },
    )
)


# ------------------------------------------------------------
# Create figure
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 7))

markers = {
    "ONT": "o",
    "PacBio": "s",
}

for (sample, technology), group in data.groupby(
    ["sample", "read_technology"],
    observed=True,
):
    group = group.sort_values("configuration")

    ax.plot(
        group["configuration"].astype(str),
        group["cigar_mapped_percent"],
        marker=markers.get(technology, "o"),
        linewidth=2,
        markersize=7,
        label=f"{sample} · {technology}",
    )


# Reference line
ax.axhline(
    100,
    linestyle="--",
    linewidth=1.2,
    alpha=0.7,
)


# ------------------------------------------------------------
# Dynamic zoomed y-axis
# Appropriate for a point/line diagnostic plot
# ------------------------------------------------------------
minimum_value = data["cigar_mapped_percent"].min()
maximum_value = data["cigar_mapped_percent"].max()

lower_limit = math.floor((minimum_value - 0.5) * 2) / 2
upper_limit = math.ceil((maximum_value + 0.5) * 2) / 2

ax.set_ylim(lower_limit, upper_limit)


# ------------------------------------------------------------
# Labels and appearance
# ------------------------------------------------------------
ax.set_title(
    "CIGAR-mapped bases by alignment configuration",
    fontsize=16,
    fontweight="bold",
    pad=16,
)

ax.set_xlabel(
    "Alignment configuration",
    fontsize=12,
)

ax.set_ylabel(
    "CIGAR-mapped bases / total read bases",
    fontsize=12,
)

ax.yaxis.set_major_formatter(
    PercentFormatter(
        xmax=100,
        decimals=1,
    )
)

ax.grid(
    axis="y",
    linestyle="--",
    alpha=0.3,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(
    title="Sample · technology",
    frameon=False,
    ncol=2,
    fontsize=9,
)

fig.text(
    0.5,
    0.01,
    (
        "Diagnostic plot using current samtools statistics. "
        "Values above 100% can indicate supplementary or secondary "
        "alignment contributions."
    ),
    ha="center",
    fontsize=9,
)

fig.tight_layout(rect=[0, 0.05, 1, 1])


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
