#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ============================================================
# Paths
# ============================================================

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

input_tsv = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv"
)

output_png = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "input_normalized_cigar_yield_percent.png"
)

output_pdf = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "input_normalized_cigar_yield_percent.pdf"
)

output_tsv = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/input_normalized_cigar_yield_percent.tsv"
)

output_png.parent.mkdir(parents=True, exist_ok=True)
output_tsv.parent.mkdir(parents=True, exist_ok=True)


# ============================================================
# Orders and colors
# ============================================================

aligner_order = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
sample_order = ["HG002", "HG003", "HG004"]
technology_order = ["ONT", "PacBio"]

sample_colors = {
    "HG002": "#0072B2",  # blue
    "HG003": "#D55E00",  # orange
    "HG004": "#009E73",  # green
}


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(input_tsv, sep="\t")

# Standardize aligner naming if needed
df["aligner"] = df["aligner"].replace({"VACMap": "VACmap"})

required_columns = {
    "sample",
    "read_technology",
    "aligner",
    "total_length",
    "bases_mapped_cigar",
}

missing = required_columns - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

# Keep only the desired rows
df = df[
    df["aligner"].isin(aligner_order)
    & df["sample"].isin(sample_order)
    & df["read_technology"].isin(technology_order)
].copy()

# Numeric conversion
df["total_length"] = pd.to_numeric(df["total_length"], errors="raise")
df["bases_mapped_cigar"] = pd.to_numeric(df["bases_mapped_cigar"], errors="raise")


# ============================================================
# Compute common input denominator
# ============================================================

# For each sample and technology, use the maximum total_length
# across aligners as a common denominator.
df["common_input_bases"] = df.groupby(
    ["sample", "read_technology"]
)["total_length"].transform("max")

df["input_normalized_cigar_yield_percent"] = (
    df["bases_mapped_cigar"] / df["common_input_bases"] * 100.0
)

# Save values table
summary_df = df[
    [
        "sample",
        "read_technology",
        "aligner",
        "total_length",
        "common_input_bases",
        "bases_mapped_cigar",
        "input_normalized_cigar_yield_percent",
    ]
].copy()

summary_df.to_csv(output_tsv, sep="\t", index=False)

# Enforce plotting order
df["aligner"] = pd.Categorical(df["aligner"], categories=aligner_order, ordered=True)
df["sample"] = pd.Categorical(df["sample"], categories=sample_order, ordered=True)
df["read_technology"] = pd.Categorical(
    df["read_technology"], categories=technology_order, ordered=True
)

df = df.sort_values(["read_technology", "aligner", "sample"])


# ============================================================
# Plot
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(16, 7.8), sharey=True)

x = np.arange(len(aligner_order))

# small spaces within bars
bar_width = 0.18
sample_spacing = 0.21

sample_offsets = {
    "HG002": -sample_spacing,
    "HG003": 0.0,
    "HG004": sample_spacing,
}

# Dynamic upper limit
ymax = df["input_normalized_cigar_yield_percent"].max()
y_upper = min(100.5, np.ceil((ymax + 0.15) * 10) / 10)

for ax, tech in zip(axes, technology_order):
    sub = df[df["read_technology"] == tech].copy()

    for sample in sample_order:
        sample_df = (
            sub[sub["sample"] == sample]
            .set_index("aligner")
            .reindex(aligner_order)
        )

        xpos = x + sample_offsets[sample]
        vals = sample_df["input_normalized_cigar_yield_percent"].values

        bars = ax.bar(
            xpos,
            vals,
            width=bar_width,
            color=sample_colors[sample],
            edgecolor="white",
            linewidth=1.0,
            zorder=3,
        )

        # Percentage labels above bars
        #for bar, value in zip(bars, vals):
        #    ax.text(
        #        bar.get_x() + bar.get_width() / 2,
        #        value + 0.03,
        #        f"{value:.2f}%",
        #        ha="center",
        #        va="bottom",
        #        fontsize=10,
        #        rotation=45,
        #        rotation_mode="anchor",
        #        clip_on=False,
        #    )

    ax.set_title("PacBio HiFi" if tech == "PacBio" else "ONT", fontsize=20, pad=24)
    ax.set_xticks(x)
    ax.set_xticklabels(aligner_order, fontsize=13)
    ax.set_xlabel("Aligner", fontsize=15)
    ax.tick_params(axis="y", labelsize=12)
    ax.grid(axis="y", alpha=0.25, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(90, y_upper)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[0].set_ylabel("Input-normalized CIGAR-aligned yield (%)", fontsize=16)

# Panel labels
axes[0].text(-0.08, 1.05, "a", transform=axes[0].transAxes,
             fontsize=22, fontweight="bold")
axes[1].text(-0.08, 1.05, "b", transform=axes[1].transAxes,
             fontsize=22, fontweight="bold")

# Legend
sample_legend = [
    Patch(facecolor=sample_colors[s], edgecolor="white", label=s)
    for s in sample_order
]

fig.legend(
    handles=sample_legend,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.02),
    ncol=3,
    title="GIAB sample",
    fontsize=13,
    title_fontsize=16,
    frameon=False
)

fig.tight_layout(rect=[0, 0, 1, 0.93])

plt.savefig(output_png, dpi=300, bbox_inches="tight")
plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
plt.show()

print("\nDirect CIGAR yield percentage figure created:")
print(output_png)
print(output_pdf)

print("\nSummary table:")
print(output_tsv)