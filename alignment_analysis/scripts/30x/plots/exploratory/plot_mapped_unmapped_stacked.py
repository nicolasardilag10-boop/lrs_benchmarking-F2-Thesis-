#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------
PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

input_tsv = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv")

output_png = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "mapped_unmapped_bases_stacked.png")

output_png.parent.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def lighten_color(color, amount=0.50):
    """
    Make a lighter version of a color by mixing it with white.
    amount=0   -> original color
    amount=1   -> white
    """
    c = np.array(mcolors.to_rgb(color))
    white = np.array([1, 1, 1])
    return tuple(c + (white - c) * amount)


# --------------------------------------------------
# Match the scientific palette used in the alignment error
# rate figures for a consistent professional look.
# --------------------------------------------------
sample_colors = {
    "HG002": "#0072B2",  # blue
    "HG003": "#D55E00",  # orange
    "HG004": "#009E73",  # green
}

aligner_order = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
sample_order = ["HG002", "HG003", "HG004"]
technology_order = ["ONT", "PacBio"]


# --------------------------------------------------
# Load data
# --------------------------------------------------
df = pd.read_csv(input_tsv, sep="\t")

# Unmapped bases
df["unmapped_bases"] = (df["total_length"] - df["bases_mapped"]).clip(lower=0)

# Convert to Gbp
df["mapped_gbp"] = df["bases_mapped"] / 1e9
df["unmapped_gbp"] = df["unmapped_bases"] / 1e9

# Enforce plotting order
df["aligner"] = pd.Categorical(df["aligner"], categories=aligner_order, ordered=True)
df["sample"] = pd.Categorical(df["sample"], categories=sample_order, ordered=True)
df["read_technology"] = pd.Categorical(df["read_technology"], categories=technology_order, ordered=True)

df = df.sort_values(["read_technology", "aligner", "sample"])


# --------------------------------------------------
# Plot
# --------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
bar_width = 0.22
x = np.arange(len(aligner_order))

for ax, tech in zip(axes, technology_order):
    sub = df[df["read_technology"] == tech].copy()

    for i, sample in enumerate(sample_order):
        sample_df = sub[sub["sample"] == sample].set_index("aligner").reindex(aligner_order)

        xpos = x + (i - 1) * bar_width

        # Use SAME sample colors as first figure
        base_color = sample_colors[sample]
        light_color = lighten_color(base_color, amount=0.50)

        mapped_vals = sample_df["mapped_gbp"].values
        unmapped_vals = sample_df["unmapped_gbp"].values

        # mapped bases
        ax.bar(
            xpos,
            mapped_vals,
            width=bar_width,
            color=base_color,
            edgecolor="black",
            linewidth=0.6
        )

        # unmapped bases
        ax.bar(
            xpos,
            unmapped_vals,
            width=bar_width,
            bottom=mapped_vals,
            color=light_color,
            edgecolor="black",
            linewidth=0.6
        )

    ax.set_title("PacBio HiFi" if tech == "PacBio" else "ONT", fontsize=20, pad=16)
    ax.set_xticks(x)
    ax.set_xticklabels(aligner_order, fontsize=13)
    ax.set_xlabel("Aligner", fontsize=15)
    ax.tick_params(axis="y", labelsize=12)

axes[0].set_ylabel("Bases (Gbp)", fontsize=16)

# Panel labels to match style of first figure
axes[0].text(-0.08, 1.05, "a", transform=axes[0].transAxes, fontsize=22, fontweight="bold")
axes[1].text(-0.08, 1.05, "b", transform=axes[1].transAxes, fontsize=22, fontweight="bold")

# Legends
sample_legend = [
    Patch(facecolor=sample_colors[s], edgecolor="black", label=s)
    for s in sample_order
]

type_legend = [
    Patch(facecolor="gray", edgecolor="black", label="Mapped bases"),
    Patch(facecolor="lightgray", edgecolor="black", label="Unmapped bases"),
]

fig.legend(
    handles=sample_legend,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.05),
    ncol=3,
    title="GIAB sample",
    fontsize=13,
    title_fontsize=16,
    frameon=False
)

fig.legend(
    handles=type_legend,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.98),
    ncol=2,
    fontsize=12,
    frameon=False
)

fig.tight_layout(rect=[0, 0, 1, 0.90])
plt.savefig(output_png, dpi=300, bbox_inches="tight")
plt.show()