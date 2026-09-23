#!/usr/bin/env python3
"""Plot observed wall-clock runtime for all four 30x benchmark aligners."""

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
OUTPUT_DATA = TABLE_DIR / "derived" / "plot_data" / "wallclock_runtime_30x_plotting_values.tsv"
OUTPUT_PNG = FIGURE_DIR / "05_runtime_30x.png"
OUTPUT_PDF = FIGURE_DIR / "05_runtime_30x.pdf"

ALIGNERS = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
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
        "configuration",
        "wallclock_runtime_hours",
        "threads",
        "wallclock_runtime_source_log",
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

    data["wallclock_runtime_hours"] = pd.to_numeric(
        data["wallclock_runtime_hours"], errors="raise"
    )
    data["threads"] = pd.to_numeric(data["threads"], errors="raise").astype(int)
    if data["wallclock_runtime_hours"].isna().any():
        rows = data.loc[
            data["wallclock_runtime_hours"].isna(),
            ["sample", "read_technology", "aligner"],
        ]
        raise ValueError("Missing runtime values:\n" + rows.to_string(index=False))
    if (data["wallclock_runtime_hours"] < 0).any():
        raise ValueError("Wall-clock runtime cannot be negative")

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

    # Keep configuration provenance in the exported plotting table while
    # plotting every runtime value recorded in the canonical benchmark table.
    if "runtime_configuration" in data.columns:
        runtime_configuration = data["runtime_configuration"].astype("string")
        data["runtime_matches_benchmark_configuration"] = (
            runtime_configuration.eq(data["configuration"].astype("string"))
        )
    else:
        data["runtime_configuration"] = pd.NA
        data["runtime_matches_benchmark_configuration"] = pd.NA

    return data


def thread_labels(data: pd.DataFrame) -> list[str]:
    return list(ALIGNERS)


def add_panel(
    axis: plt.Axes,
    subset: pd.DataFrame,
    x_max: float,
    labels: list[str],
) -> None:
    y = np.arange(len(ALIGNERS), dtype=float)
    bar_height = 0.18
    sample_spacing = 0.21
    offsets = [-sample_spacing, 0.0, sample_spacing]

    for sample, offset in zip(SAMPLES, offsets):
        sample_data = (
            subset[subset["sample"] == sample]
            .set_index("aligner")
            .reindex(ALIGNERS)
        )
        values = sample_data["wallclock_runtime_hours"].to_numpy(dtype=float)
        valid = np.isfinite(values)
        bars = axis.barh(
            y[valid] + offset,
            values[valid],
            height=bar_height,
            color=SAMPLE_COLORS[sample],
            edgecolor="#222222",
            linewidth=0.45,
            zorder=3,
        )

        for bar, value in zip(bars, values[valid]):
            axis.text(
                value + x_max * 0.012,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.1f}",
                ha="left",
                va="center",
                fontsize=8,
                color="#222222",
            )

    axis.set_yticks(y, labels)
    axis.invert_yaxis()
    axis.set_xlim(0, x_max)
    axis.grid(axis="x", color="#D9D9D9", linewidth=0.7, alpha=0.75, zorder=0)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=9)


def main() -> int:
    data = load_data()
    output_columns = [
        "sample",
        "read_technology",
        "aligner",
        "configuration",
        "runtime_configuration",
        "runtime_matches_benchmark_configuration",
        "wallclock_runtime_hours",
        "wallclock_runtime_seconds",
        "wallclock_runtime_hms",
        "wallclock_runtime_source_log",
        "threads",
    ]
    output_columns = [column for column in output_columns if column in data.columns]
    plot_data = data[output_columns].copy()
    plot_data["read_technology"] = pd.Categorical(
        plot_data["read_technology"], TECHNOLOGIES, ordered=True
    )
    plot_data["aligner"] = pd.Categorical(
        plot_data["aligner"], ALIGNERS, ordered=True
    )
    plot_data["sample"] = pd.Categorical(
        plot_data["sample"], SAMPLES, ordered=True
    )
    plot_data = plot_data.sort_values(["read_technology", "aligner", "sample"])

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT_DATA, sep="\t", index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 12,
            "axes.titlesize": 15,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(2, 1, figsize=(7.9, 7.6), sharex=True)
    x_max = float(plot_data["wallclock_runtime_hours"].max()) * 1.18
    labels = thread_labels(data)

    for panel_index, (axis, technology) in enumerate(zip(axes, TECHNOLOGIES)):
        subset = plot_data[plot_data["read_technology"] == technology]
        add_panel(axis, subset, x_max, labels)
        axis.set_title("PacBio HiFi" if technology == "PacBio" else "ONT", pad=10)
        axis.text(
            -0.10,
            1.03,
            "ab"[panel_index],
            transform=axis.transAxes,
            fontsize=15,
            fontweight="bold",
            va="bottom",
        )

    figure.supxlabel("Wall-clock runtime (hours)", fontsize=11, y=0.035)
    legend_handles = [
        Patch(facecolor=SAMPLE_COLORS[sample], edgecolor="#222222", label=sample)
        for sample in SAMPLES
    ]
    figure.legend(
        handles=legend_handles,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        frameon=False,
    )
    figure.subplots_adjust(left=0.16, right=0.96, top=0.87, bottom=0.10, hspace=0.4)

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
