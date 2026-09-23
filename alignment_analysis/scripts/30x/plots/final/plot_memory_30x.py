#!/usr/bin/env python3
"""Create a three-panel 30x memory-resource figure.

Panel a: configured RAM limit (GB) per aligner, one bar per sample.
Panels b/c: measured peak RSS (GB) for minimap2/pbmm2, ONT and PacBio HiFi.
"""

from __future__ import annotations

from pathlib import Path
import sys

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
sys.path.insert(0, str(PROJECT / "alignment_analysis" / "scripts"))

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from utils.plot_style import (
    ALIGNER_ORDER,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    apply_style,
    panel_letter,
    save_figure,
)

MEASURED_ALIGNERS = ["minimap2", "pbmm2"]

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "memory_30x_plotting_values.tsv"
OUTPUT_PNG = FIGURE_DIR / "08_memory_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

BAR_HEIGHT = 0.18
SAMPLE_SPACING = 0.21
SAMPLE_OFFSETS = {"HG002": -SAMPLE_SPACING, "HG003": 0.0, "HG004": SAMPLE_SPACING}


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.is_file():
        raise FileNotFoundError(f"Benchmark table not found: {INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    required = {"sample", "read_technology", "aligner", "peak_ram_gb", "command_ram_limit_gb"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap", "VG": "VG Giraffe", "vg": "VG Giraffe"})
    data["read_technology"] = data["read_technology"].replace(
        {"PacBio HiFi": "PacBio", "PB": "PacBio", "pb": "PacBio", "ont": "ONT"}
    )
    data = data[
        data["sample"].isin(SAMPLE_ORDER)
        & data["read_technology"].isin(TECHNOLOGY_ORDER)
        & data["aligner"].isin(ALIGNER_ORDER)
    ].copy()

    keys = ["sample", "read_technology", "aligner"]
    if data.duplicated(keys).any():
        duplicates = data[data.duplicated(keys, keep=False)][keys]
        raise ValueError("Duplicate sample/technology/aligner rows:\n" + duplicates.to_string(index=False))

    expected = pd.MultiIndex.from_product([SAMPLE_ORDER, TECHNOLOGY_ORDER, ALIGNER_ORDER], names=keys)
    observed = pd.MultiIndex.from_frame(data[keys])
    missing_rows = expected.difference(observed).tolist()
    if missing_rows:
        raise ValueError(f"Missing benchmark combinations: {missing_rows}")

    for column in ["peak_ram_gb", "command_ram_limit_gb"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    if data["command_ram_limit_gb"].isna().any():
        bad = data[data["command_ram_limit_gb"].isna()][["sample", "read_technology", "aligner"]]
        raise ValueError("Configured RAM limits must be present for all rows:\n" + bad.to_string(index=False))

    measured = data[data["peak_ram_gb"].notna()]
    measured_combinations = set(
        measured[["sample", "read_technology", "aligner"]].itertuples(index=False, name=None)
    )
    expected_measured = {
        (sample, technology, aligner)
        for sample in SAMPLE_ORDER
        for technology in TECHNOLOGY_ORDER
        for aligner in MEASURED_ALIGNERS
    }
    if measured_combinations != expected_measured:
        unexpected = measured_combinations - expected_measured
        missing_measured = expected_measured - measured_combinations
        raise ValueError(
            "Peak RSS availability differs from expected.\n"
            f"Missing: {sorted(missing_measured)}\nUnexpected: {sorted(unexpected)}"
        )

    return data


def prepare_plot_data(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    limit_counts = data.groupby("aligner")["command_ram_limit_gb"].nunique()
    if not limit_counts.reindex(ALIGNER_ORDER).eq(1).all():
        raise ValueError("Configured RAM limits vary between runs for an aligner")

    limits = (
        data.groupby(["sample", "aligner"], as_index=False)["command_ram_limit_gb"]
        .first()
        .set_index(["sample", "aligner"])
        .reindex(pd.MultiIndex.from_product([SAMPLE_ORDER, ALIGNER_ORDER], names=["sample", "aligner"]))
        .reset_index()
    )
    limits["panel"] = "a"
    limits["panel_description"] = "Configured RAM limit"
    limits["read_technology"] = "ONT and PacBio"
    limits["memory_gb"] = limits["command_ram_limit_gb"]
    limits["memory_value_type"] = "Configured RAM limit"
    limits["source_column"] = "command_ram_limit_gb"

    measured = data[data["peak_ram_gb"].notna()].copy()
    measured["panel"] = measured["read_technology"].map({"ONT": "b", "PacBio": "c"})
    measured["panel_description"] = measured["read_technology"].map(
        {"ONT": "Measured peak RSS — ONT", "PacBio": "Measured peak RSS — PacBio HiFi"}
    )
    measured["memory_gb"] = measured["peak_ram_gb"]
    measured["memory_value_type"] = "Measured peak RSS"
    measured["source_column"] = "peak_ram_gb"

    columns = [
        "panel", "panel_description", "sample", "read_technology", "aligner",
        "memory_gb", "memory_value_type", "source_column",
    ]
    export = pd.concat([limits[columns], measured[columns]], ignore_index=True)
    return limits, export


def style_axis(axis: plt.Axes, x_max: float) -> None:
    axis.set_xlim(0, x_max)
    axis.grid(axis="x", color="#E5E5E5", linewidth=0.45, zorder=0)
    axis.set_axisbelow(True)
    for spine_name, spine in axis.spines.items():
        spine.set_visible(spine_name in ("left", "bottom"))
    axis.spines["left"].set_linewidth(0.4)
    axis.spines["bottom"].set_linewidth(0.4)
    axis.tick_params(axis="both", labelsize=12)


def add_limit_panel(axis: plt.Axes, limits: pd.DataFrame, x_max: float) -> None:
    y = np.arange(len(ALIGNER_ORDER), dtype=float)

    for sample in SAMPLE_ORDER:
        sample_data = limits[limits["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
        values = sample_data["memory_gb"].to_numpy(dtype=float)
        bars = axis.barh(
            y + SAMPLE_OFFSETS[sample], values, height=BAR_HEIGHT,
            color=SAMPLE_COLORS[sample], edgecolor="#222222", linewidth=0.5, hatch="//", zorder=3,
        )
        for bar, value in zip(bars, values):
            axis.text(
                value + x_max * 0.012, bar.get_y() + bar.get_height() / 2,
                f"{value:.0f}", ha="left", va="center", fontsize=6.5,
            )

    axis.set_yticks(y, ALIGNER_ORDER)
    axis.invert_yaxis()
    axis.set_title("Configured RAM limit", pad=6)
    style_axis(axis, x_max)


def add_measured_panel(axis: plt.Axes, data: pd.DataFrame, technology: str, x_max: float) -> None:
    subset = data[data["read_technology"] == technology]
    y = np.arange(len(MEASURED_ALIGNERS), dtype=float)

    for sample in SAMPLE_ORDER:
        sample_data = (
            subset[subset["sample"] == sample].set_index("aligner").reindex(MEASURED_ALIGNERS)
        )
        values = sample_data["peak_ram_gb"].to_numpy(dtype=float)
        bars = axis.barh(
            y + SAMPLE_OFFSETS[sample], values, height=BAR_HEIGHT,
            color=SAMPLE_COLORS[sample], edgecolor="#222222", linewidth=0.4, zorder=3,
        )
        for bar, value in zip(bars, values):
            axis.text(
                value + x_max * 0.012, bar.get_y() + bar.get_height() / 2,
                f"{value:.1f}", ha="left", va="center", fontsize=6.5,
            )

    axis.set_yticks(y, MEASURED_ALIGNERS)
    axis.invert_yaxis()
    label = "PacBio HiFi" if technology == "PacBio" else "ONT"
    axis.set_title(f"Measured peak RSS — {label}", pad=6)
    style_axis(axis, x_max)


def main() -> int:
    data = load_data()
    limits, plot_data = prepare_plot_data(data)

    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    apply_style()
    matplotlib.rcParams["hatch.linewidth"] = 0.6
    matplotlib.rcParams["hatch.color"] = "#868181"

    figure, axes = plt.subplots(
        3, 1, figsize=(7.16, 6.6), sharex=True,
        gridspec_kw={"height_ratios": [1.55, 1.0, 1.0]},
    )

    maximum = float(plot_data["memory_gb"].max())
    x_max = np.ceil((maximum * 1.10) / 20.0) * 20.0

    add_limit_panel(axes[0], limits, x_max)
    add_measured_panel(axes[1], data, "ONT", x_max)
    add_measured_panel(axes[2], data, "PacBio", x_max)

    for index, axis in enumerate(axes):
        panel_letter(axis, "abc"[index], x=-0.10, y=1.03, fontsize=12)

    limit_handle = Patch(facecolor="white", edgecolor="#868181", hatch="//", label="Configured limit")
    axes[0].legend(handles=[limit_handle], loc="upper right", frameon=False)

    figure.supxlabel("Memory (GB)", y=0.02, fontsize=12)
    figure.subplots_adjust(left=0.16, right=0.96, top=0.84, bottom=0.10, hspace=0.55)

    axes_center_x = (axes[0].get_position().x0 + axes[0].get_position().x1) / 2
    sample_legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="#222222", linewidth=0.4, label=sample)
        for sample in SAMPLE_ORDER
    ]
    figure.legend(
        handles=sample_legend_handles, ncol=3, loc="upper center",
        bbox_to_anchor=(axes_center_x, 0.99), frameon=False,
    )

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_DATA}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
