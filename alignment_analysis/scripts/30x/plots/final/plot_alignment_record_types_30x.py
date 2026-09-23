#!/usr/bin/env python3
"""Plot secondary and supplementary alignment-record rates for the 30x benchmark."""

from __future__ import annotations

import os
from pathlib import Path
import sys

MPL_CACHE = Path("/tmp/matplotlib-lrs-benchmarking")
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import numpy as np
    import pandas as pd
except ModuleNotFoundError as error:
    raise SystemExit(
        f"Missing plotting package: {error.name}\nPython: {sys.executable}"
    ) from error


def find_project_root() -> Path:
    """Find the repository without depending on a username or script depth."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "alignment_analysis").is_dir():
            return parent
    raise RuntimeError("Could not locate the lrs_benchmarking project root")


PROJECT = find_project_root()
TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "alignment_record_types_30x.tsv"
OUTPUT_PNG = FIGURE_DIR / "07_alignment_record_types_30x.png"
OUTPUT_PDF = FIGURE_DIR / "07_alignment_record_types_30x.pdf"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
TECHNOLOGIES = ["ONT", "PacBio"]
SAMPLES = ["HG002", "HG003", "HG004"]
SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}
METRICS = [
    (
        "supplementary_alignments_per_100_mapped_reads",
        "Supplementary alignments\nper 100 mapped reads",
    ),
    (
        "secondary_alignments_per_100_mapped_reads",
        "Secondary alignments\nper 100 mapped reads",
    ),
]


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.is_file():
        raise FileNotFoundError(f"Benchmark table not found: {INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    required = {
        "sample",
        "read_technology",
        "aligner",
        *(metric for metric, _ in METRICS),
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})
    data["read_technology"] = data["read_technology"].replace(
        {"PacBio HiFi": "PacBio", "PB": "PacBio"}
    )
    data = data[
        data["sample"].isin(SAMPLES)
        & data["read_technology"].isin(TECHNOLOGIES)
        & data["aligner"].isin(ALIGNERS)
    ].copy()

    for metric, _ in METRICS:
        data[metric] = pd.to_numeric(data[metric], errors="raise")

    keys = ["sample", "read_technology", "aligner"]
    duplicates = data.duplicated(keys, keep=False)
    if duplicates.any():
        raise ValueError(
            "Duplicate sample/technology/aligner rows:\n"
            + data.loc[duplicates, keys].to_string(index=False)
        )

    expected = pd.MultiIndex.from_product(
        [SAMPLES, TECHNOLOGIES, ALIGNERS], names=keys
    )
    observed = pd.MultiIndex.from_frame(data[keys])
    missing_combinations = expected.difference(observed)
    if len(missing_combinations):
        raise ValueError(f"Missing benchmark combinations: {list(missing_combinations)}")

    return data


def add_bars(axis, subset: pd.DataFrame, metric: str, y_max: float) -> None:
    x = np.arange(len(ALIGNERS), dtype=float)
    width = 0.23
    offsets = [-width, 0.0, width]

    for sample, offset in zip(SAMPLES, offsets):
        sample_data = (
            subset[subset["sample"] == sample]
            .set_index("aligner")
            .reindex(ALIGNERS)
        )
        values = sample_data[metric].to_numpy(dtype=float)
        bars = axis.bar(
            x + offset,
            values,
            width=width,
            color=SAMPLE_COLORS[sample],
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
        )

        for bar, value in zip(bars, values):
            label = "0" if value == 0 else f"{value:.1f}"
            height = max(value, y_max * 0.008)
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                height + y_max * 0.018,
                label,
                ha="center",
                va="bottom",
                rotation=65 if value > 0 else 0,
                fontsize=7.5,
                color="#222222",
            )

    axis.set_xticks(x, ALIGNERS)
    axis.set_ylim(0, y_max)
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.75, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(axis="x", labelrotation=0)


def main() -> int:
    data = load_data()
    output_columns = [
        "sample",
        "read_technology",
        "aligner",
        *(metric for metric, _ in METRICS),
    ]
    plot_data = data[output_columns].sort_values(
        ["read_technology", "aligner", "sample"]
    )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.2, 8.4),
        sharex="col",
        sharey="row",
        constrained_layout=False,
    )

    row_limits = []
    for metric, _ in METRICS:
        maximum = float(plot_data[metric].max())
        row_limits.append(maximum * 1.18 if maximum > 0 else 1.0)

    panel_letters = iter("abcd")
    for row_index, (metric, y_label) in enumerate(METRICS):
        for column_index, technology in enumerate(TECHNOLOGIES):
            axis = axes[row_index, column_index]
            subset = plot_data[plot_data["read_technology"] == technology]
            add_bars(axis, subset, metric, row_limits[row_index])
            axis.text(
                -0.10,
                1.04,
                next(panel_letters),
                transform=axis.transAxes,
                fontsize=15,
                fontweight="bold",
                va="bottom",
            )
            if row_index == 0:
                title = "PacBio HiFi" if technology == "PacBio" else technology
                axis.set_title(title, pad=10)
            if column_index == 0:
                axis.set_ylabel(y_label)

    legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="white", label=sample)
        for sample in SAMPLES
    ]
    figure.legend(
        handles=legend_handles,
        title="GIAB sample",
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        frameon=False,
    )
    figure.suptitle(
        "Supplementary and secondary alignment records",
        fontsize=15,
        fontweight="bold",
        y=0.925,
    )
    figure.text(
        0.5,
        0.012,
        "Counts reflect SAM flags and aligner output policies; zero is not a universal measure of mapping quality.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#444444",
        style="italic",
    )
    figure.subplots_adjust(left=0.09, right=0.985, top=0.85, bottom=0.12, hspace=0.34, wspace=0.13)

    figure.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(OUTPUT_PDF, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_DATA}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
