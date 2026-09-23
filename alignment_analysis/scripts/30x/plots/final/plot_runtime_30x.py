#!/usr/bin/env python3
"""Plot observed wall-clock runtime (hours) for the 30x benchmark aligners.

Grouped horizontal bars per sample per aligner, from a zero baseline.
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

from utils.plot_style import (
    ALIGNER_ORDER,
    FULL_WIDTH_IN,
    SAMPLE_COLORS,
    SAMPLE_ORDER,
    TECHNOLOGY_ORDER,
    TECHNOLOGY_TITLES,
    apply_style,
    panel_letter,
    sample_legend_handles,
    save_figure,
)

TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
FIGURE_DIR = PROJECT / "alignment_analysis" / "figures" / "30x" / "final"

INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "wallclock_runtime_30x_plotting_values.tsv"
OUTPUT_PNG = FIGURE_DIR / "05_runtime_30x.png"
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

BAR_HEIGHT = 0.18
SAMPLE_OFFSETS = {"HG002": 0.21, "HG003": 0.0, "HG004": -0.21}


def load_data() -> pd.DataFrame:
    if not INPUT_TSV.is_file():
        raise FileNotFoundError(f"Benchmark table not found: {INPUT_TSV}")

    data = pd.read_csv(INPUT_TSV, sep="\t")
    required = {
        "sample", "read_technology", "aligner", "configuration",
        "wallclock_runtime_hours", "threads", "wallclock_runtime_source_log",
    }
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

    data["wallclock_runtime_hours"] = pd.to_numeric(data["wallclock_runtime_hours"], errors="raise")
    data["threads"] = pd.to_numeric(data["threads"], errors="raise").astype(int)
    if data["wallclock_runtime_hours"].isna().any():
        rows = data.loc[data["wallclock_runtime_hours"].isna(), ["sample", "read_technology", "aligner"]]
        raise ValueError("Missing runtime values:\n" + rows.to_string(index=False))
    if (data["wallclock_runtime_hours"] < 0).any():
        raise ValueError("Wall-clock runtime cannot be negative")

    keys = ["sample", "read_technology", "aligner"]
    duplicates = data.duplicated(keys, keep=False)
    if duplicates.any():
        raise ValueError("Duplicate sample/technology/aligner rows:\n" + data.loc[duplicates, keys].to_string(index=False))

    expected = pd.MultiIndex.from_product([SAMPLE_ORDER, TECHNOLOGY_ORDER, ALIGNER_ORDER], names=keys)
    observed = pd.MultiIndex.from_frame(data[keys])
    missing_combinations = expected.difference(observed)
    if len(missing_combinations):
        raise ValueError(f"Missing benchmark combinations: {list(missing_combinations)}")

    if "runtime_configuration" in data.columns:
        runtime_configuration = data["runtime_configuration"].astype("string")
        data["runtime_matches_benchmark_configuration"] = runtime_configuration.eq(data["configuration"].astype("string"))
    else:
        data["runtime_configuration"] = pd.NA
        data["runtime_matches_benchmark_configuration"] = pd.NA

    return data


def add_panel(axis: plt.Axes, subset: pd.DataFrame, x_max: float, show_labels: bool) -> None:
    y_positions = np.arange(len(ALIGNER_ORDER), dtype=float)

    for sample in SAMPLE_ORDER:
        sample_data = (
            subset[subset["sample"] == sample].set_index("aligner").reindex(ALIGNER_ORDER)
        )
        values = sample_data["wallclock_runtime_hours"].to_numpy(dtype=float)
        valid = np.isfinite(values)
        y = y_positions[valid] + SAMPLE_OFFSETS[sample]

        axis.barh(
            y, values[valid],
            height=BAR_HEIGHT,
            color=SAMPLE_COLORS[sample], edgecolor="white", linewidth=0.6,
            zorder=3,
        )

    axis.set_yticks(y_positions, ALIGNER_ORDER)
    axis.invert_yaxis()
    axis.set_ylim(len(ALIGNER_ORDER) - 0.5, -0.5)
    axis.set_xlim(0, x_max)
    axis.grid(axis="x", color="#E5E5E5", linewidth=0.45, zorder=0)
    axis.set_axisbelow(True)
    for spine_name, spine in axis.spines.items():
        spine.set_visible(spine_name in ("left", "bottom"))
    axis.spines["left"].set_linewidth(0.7)
    axis.spines["bottom"].set_linewidth(0.7)
    axis.tick_params(axis="both", labelsize=12)
    if not show_labels:
        axis.tick_params(axis="y", left=False, labelleft=False)


def main() -> int:
    data = load_data()
    output_columns = [
        "sample", "read_technology", "aligner", "configuration",
        "runtime_configuration", "runtime_matches_benchmark_configuration",
        "wallclock_runtime_hours", "wallclock_runtime_seconds",
        "wallclock_runtime_hms", "wallclock_runtime_source_log", "threads",
    ]
    output_columns = [column for column in output_columns if column in data.columns]
    plot_data = data[output_columns].copy()
    plot_data["read_technology"] = pd.Categorical(plot_data["read_technology"], TECHNOLOGY_ORDER, ordered=True)
    plot_data["aligner"] = pd.Categorical(plot_data["aligner"], ALIGNER_ORDER, ordered=True)
    plot_data["sample"] = pd.Categorical(plot_data["sample"], SAMPLE_ORDER, ordered=True)
    plot_data = plot_data.sort_values(["read_technology", "aligner", "sample"])

    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    apply_style()

    figure, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(FULL_WIDTH_IN, 3.0))
    x_max = float(plot_data["wallclock_runtime_hours"].max()) * 1.15

    for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
        axis = axes[panel_index]
        subset = plot_data[plot_data["read_technology"] == technology]
        add_panel(axis, subset, x_max, show_labels=technology == "ONT")
        axis.set_title(TECHNOLOGY_TITLES[technology], pad=6)
        panel_letter(axis, "a" if technology == "ONT" else "b", x=-0.05, fontsize=12)

    figure.supxlabel("Wall-clock runtime (hours)", x=0.55,y=0.02, fontsize=12)
    figure.legend(
        handles=sample_legend_handles(),
        frameon=False,
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.08),
        columnspacing=1.3,
        handletextpad=0.4,
    )
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(left=0.17, right=0.97, top=0.86, bottom=0.15, wspace=0.15)

    save_figure(figure, OUTPUT_PNG, OUTPUT_PDF)

    print(f"Input: {INPUT_TSV}")
    print(f"Plot data: {OUTPUT_DATA}")
    print(f"PNG: {OUTPUT_PNG}")
    print(f"PDF: {OUTPUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
