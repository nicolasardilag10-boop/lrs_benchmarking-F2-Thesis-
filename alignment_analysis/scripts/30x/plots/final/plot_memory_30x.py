#!/usr/bin/env python3
"""Create a three-panel 30x memory-resource figure."""

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
    for parent in Path(__file__).resolve().parents:
        if (parent / "alignment_analysis").is_dir():
            return parent
    raise RuntimeError("Could not locate the lrs_benchmarking project root")


PROJECT = find_project_root()
TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"
INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "memory_30x_plotting_values.tsv"
OUTPUT_PNG = FIGURE_DIR / "08_memory_30x.png"
OUTPUT_PDF = FIGURE_DIR / "08_memory_30x.pdf"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
MEASURED_ALIGNERS = ["minimap2", "pbmm2"]
TECHNOLOGIES = ["ONT", "PacBio"]
SAMPLES = ["HG002", "HG003", "HG004"]
SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.is_file():
        raise FileNotFoundError(f"Benchmark table not found: {INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    required = {
        "sample",
        "read_technology",
        "aligner",
        "peak_ram_gb",
        "command_ram_limit_gb",
    }
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data["aligner"] = data["aligner"].replace(
        {"VACMap": "VACmap", "VG": "VG Giraffe", "vg": "VG Giraffe"}
    )
    data["read_technology"] = data["read_technology"].replace(
        {"PacBio HiFi": "PacBio", "PB": "PacBio", "pb": "PacBio", "ont": "ONT"}
    )
    data = data[
        data["sample"].isin(SAMPLES)
        & data["read_technology"].isin(TECHNOLOGIES)
        & data["aligner"].isin(ALIGNERS)
    ].copy()

    keys = ["sample", "read_technology", "aligner"]
    if data.duplicated(keys).any():
        raise ValueError("Duplicate sample/technology/aligner rows")
    expected = pd.MultiIndex.from_product(
        [SAMPLES, TECHNOLOGIES, ALIGNERS], names=keys
    )
    observed = pd.MultiIndex.from_frame(data[keys])
    if missing_rows := expected.difference(observed).tolist():
        raise ValueError(f"Missing benchmark combinations: {missing_rows}")

    for column in ["peak_ram_gb", "command_ram_limit_gb"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    if data["command_ram_limit_gb"].isna().any():
        raise ValueError("Configured RAM limits must be present for all rows")

    measured = data[data["peak_ram_gb"].notna()]
    measured_combinations = set(
        measured[["sample", "read_technology", "aligner"]]
        .itertuples(index=False, name=None)
    )
    expected_measured = {
        (sample, technology, aligner)
        for sample in SAMPLES
        for technology in TECHNOLOGIES
        for aligner in MEASURED_ALIGNERS
    }
    if measured_combinations != expected_measured:
        raise ValueError("Peak RSS availability differs from the expected 12 rows")
    return data


def prepare_plot_data(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    limit_counts = data.groupby("aligner")["command_ram_limit_gb"].nunique()
    if not limit_counts.reindex(ALIGNERS).eq(1).all():
        raise ValueError("Configured RAM limits vary between runs for an aligner")

    limits = (
        data.groupby(["sample", "aligner"], as_index=False)["command_ram_limit_gb"]
        .first()
        .set_index(["sample", "aligner"])
        .reindex(pd.MultiIndex.from_product([SAMPLES, ALIGNERS], names=["sample", "aligner"]))
        .reset_index()
    )
    limits["panel"] = "A"
    limits["panel_description"] = "Configured RAM limit"
    limits["read_technology"] = "ONT and PacBio"
    limits["memory_gb"] = limits["command_ram_limit_gb"]
    limits["memory_value_type"] = "Configured RAM limit"
    limits["source_column"] = "command_ram_limit_gb"

    measured = data[data["peak_ram_gb"].notna()].copy()
    measured["panel"] = measured["read_technology"].map(
        {"ONT": "B", "PacBio": "C"}
    )
    measured["panel_description"] = measured["read_technology"].map(
        {
            "ONT": "Measured peak RSS — ONT",
            "PacBio": "Measured peak RSS — PacBio HiFi",
        }
    )
    measured["memory_gb"] = measured["peak_ram_gb"]
    measured["memory_value_type"] = "Measured peak RSS"
    measured["source_column"] = "peak_ram_gb"

    columns = [
        "panel",
        "panel_description",
        "sample",
        "read_technology",
        "aligner",
        "memory_gb",
        "memory_value_type",
        "source_column",
    ]
    export = pd.concat([limits[columns], measured[columns]], ignore_index=True)
    return limits, export


def style_axis(axis: plt.Axes, x_max: float) -> None:
    axis.set_xlim(0, x_max)
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.75, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=9)


def add_limit_panel(axis: plt.Axes, limits: pd.DataFrame, x_max: float) -> None:
    y = np.arange(len(ALIGNERS), dtype=float)
    bar_height = 0.18
    sample_spacing = 0.21
    offsets = {
        "HG002": -sample_spacing,
        "HG003": 0.0,
        "HG004": sample_spacing,
    }
    for sample in SAMPLES:
        sample_data = (
            limits[limits["sample"] == sample]
            .set_index("aligner")
            .reindex(ALIGNERS)
        )
        values = sample_data["memory_gb"].to_numpy(dtype=float)
        bars = axis.barh(
            y + offsets[sample],
            values,
            height=bar_height,
            color=SAMPLE_COLORS[sample],
            edgecolor="#222222",
            linewidth=0.65,
            hatch="////",
            zorder=3,
        )
        for bar, value in zip(bars, values):
            axis.text(
                value + x_max * 0.012,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.0f}",
                ha="left",
                va="center",
                fontsize=8,
            )
    axis.set_yticks(y, ALIGNERS)
    axis.invert_yaxis()
    axis.set_title("Configured RAM limit", pad=10)
    style_axis(axis, x_max)


def add_measured_panel(
    axis: plt.Axes,
    data: pd.DataFrame,
    technology: str,
    x_max: float,
) -> None:
    subset = data[data["read_technology"] == technology]
    y = np.arange(len(MEASURED_ALIGNERS), dtype=float)
    bar_height = 0.18
    sample_spacing = 0.21
    offsets = {
        "HG002": -sample_spacing,
        "HG003": 0.0,
        "HG004": sample_spacing,
    }
    for sample in SAMPLES:
        sample_data = (
            subset[subset["sample"] == sample]
            .set_index("aligner")
            .reindex(MEASURED_ALIGNERS)
        )
        values = sample_data["peak_ram_gb"].to_numpy(dtype=float)
        bars = axis.barh(
            y + offsets[sample],
            values,
            height=bar_height,
            color=SAMPLE_COLORS[sample],
            edgecolor="#222222",
            linewidth=0.45,
            zorder=3,
        )
        for bar, value in zip(bars, values):
            axis.text(
                value + x_max * 0.012,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.1f}",
                ha="left",
                va="center",
                fontsize=8,
            )
    axis.set_yticks(y, MEASURED_ALIGNERS)
    axis.invert_yaxis()
    label = "PacBio HiFi" if technology == "PacBio" else "ONT"
    axis.set_title(f"Measured peak RSS — {label}", pad=10)
    style_axis(axis, x_max)


def main() -> int:
    data = load_data()
    limits, plot_data = prepare_plot_data(data)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(
        3,
        1,
        figsize=(10.8, 12.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.55, 1.0, 1.0]},
    )
    maximum = float(plot_data["memory_gb"].max())
    x_max = np.ceil((maximum * 1.10) / 20.0) * 20.0
    add_limit_panel(axes[0], limits, x_max)
    add_measured_panel(axes[1], data, "ONT", x_max)
    add_measured_panel(axes[2], data, "PacBio", x_max)

    for index, axis in enumerate(axes):
        axis.text(
            -0.08,
            1.02,
            "ABC"[index] + ".",
            transform=axis.transAxes,
            fontsize=15,
            fontweight="bold",
            va="bottom",
        )

    # The hatch pattern only appears in panel A, so its key belongs on
    # that panel rather than in the figure-wide "GIAB sample" legend,
    # where it previously read as if it applied to all three panels.
    limit_handle = Patch(
        facecolor="white",
        edgecolor="#222222",
        hatch="////",
        label="Configured limit",
    )
    axes[0].legend(
        handles=[limit_handle],
        loc="upper right",
        frameon=False,
        fontsize=9,
    )
    figure.supxlabel("Memory (GB)", fontsize=11, y=0.055)
    figure.text(
        0.5,
        0.012,
        "Panel A reports configured limits; panels B–C report observed peak RSS. Limits are resource allocations, not measured consumption.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#444444",
        style="italic",
    )
    figure.subplots_adjust(
        left=0.16,
        right=0.96,
        top=0.90,
        bottom=0.10,
        hspace=0.42,
    )

    # Center the sample legend on the axes' actual horizontal span
    # (not the full figure canvas), since the left margin reserved for
    # the y-axis tick labels shifts the plot area's true center right
    # of figure-x=0.5.
    axes_center_x = (
        axes[0].get_position().x0 + axes[0].get_position().x1
    ) / 2

    sample_legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="#222222", label=sample)
        for sample in SAMPLES
    ]
    figure.legend(
        handles=sample_legend_handles,
        title="GIAB sample",
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(axes_center_x, 0.99),
        frameon=False,
    )

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
