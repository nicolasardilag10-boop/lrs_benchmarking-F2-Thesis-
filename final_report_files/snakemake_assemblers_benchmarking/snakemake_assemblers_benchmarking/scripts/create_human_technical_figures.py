#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLE_DIR = RESULTS / "tables"
FIGURE_DIR = RESULTS / "figures" / "human_technical"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

HUMAN_SAMPLES = {"HG002", "HG003", "HG004"}


def calculate_n50(lengths: list[int]) -> int:
    if not lengths:
        return 0

    total = sum(lengths)
    halfway = total / 2
    cumulative = 0

    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= halfway:
            return length

    return 0


def fasta_metrics(path: Path) -> dict[str, int]:
    lengths: list[int] = []
    current = 0
    found_header = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if found_header and current > 0:
                    lengths.append(current)

                found_header = True
                current = 0
            else:
                if not found_header:
                    raise ValueError(f"Invalid FASTA: {path}")

                current += len(line)

    if found_header and current > 0:
        lengths.append(current)

    if not lengths:
        raise ValueError(f"No non-empty sequences: {path}")

    return {
        "contig_count": len(lengths),
        "total_length_bp": sum(lengths),
        "longest_contig_bp": max(lengths),
        "n50_bp": calculate_n50(lengths),
    }


def infer_metadata(path: Path) -> tuple[str, str, str] | None:
    text = path.as_posix()

    match = re.search(
        r"(HG00[234])[._](ont|pb)",
        text,
        flags=re.IGNORECASE,
    )

    if not match:
        return None

    sample = match.group(1).upper()
    technology = match.group(2).lower()

    lower = text.lower()

    if "/ntlink/" in lower:
        tool = "ntLink"
    elif "/goldrush/" in lower:
        tool = "GoldRush"
    elif "/flye/" in lower or ".flye/" in lower:
        tool = "Flye"
    elif "/verkko/" in lower:
        tool = "Verkko"
    else:
        return None

    return tool, sample, technology


def discover_human_outputs() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    for path in RESULTS.rglob("assembly.fasta"):
        if not path.is_file() or path.stat().st_size == 0:
            continue

        relative_parts = {
            part.lower() for part in path.relative_to(ROOT).parts
        }

        if "work" in relative_parts:
            continue

        if "manual_test" in relative_parts:
            continue

        metadata = infer_metadata(path)

        if metadata is None:
            continue

        tool, sample, technology = metadata

        if sample not in HUMAN_SAMPLES:
            continue

        try:
            metrics = fasta_metrics(path)
        except ValueError as error:
            print(f"Skipping {path}: {error}")
            continue

        rows.append(
            {
                "tool": tool,
                "sample": sample,
                "technology": technology,
                "dataset_type": "human_1k_or_smoke_test",
                "assembly_path": str(path.relative_to(ROOT)),
                "file_size_bytes": path.stat().st_size,
                **metrics,
            }
        )

    rows.sort(
        key=lambda row: (
            str(row["technology"]),
            str(row["sample"]),
            str(row["tool"]),
        )
    )

    return rows


def format_bp(value: float, _: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"

    if value >= 1_000:
        return f"{value / 1_000:.0f}K"

    return f"{value:.0f}"


def save(fig: plt.Figure, name: str) -> None:
    png = FIGURE_DIR / f"{name}.png"
    pdf = FIGURE_DIR / f"{name}.pdf"

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Created {png.relative_to(ROOT)}")
    print(f"Created {pdf.relative_to(ROOT)}")


def label(row: dict[str, object]) -> str:
    technology = str(row["technology"]).upper()

    return (
        f"{row['sample']} | "
        f"{technology} | "
        f"{row['tool']}"
    )


def plot_metric(
    rows: list[dict[str, object]],
    metric: str,
    title: str,
    axis_label: str,
    filename: str,
    log_scale: bool = False,
) -> None:
    if not rows:
        return

    labels = [label(row) for row in rows]
    values = [int(row[metric]) for row in rows]

    height = max(5.5, len(rows) * 0.48 + 2.5)
    fig, ax = plt.subplots(figsize=(11.5, height))

    bars = ax.barh(labels, values)

    if log_scale and all(value > 0 for value in values):
        ax.set_xscale("log")
        axis_label += " — logarithmic scale"

    ax.set_xlabel(axis_label)
    ax.set_ylabel("Human technical output")
    ax.set_title(
        f"{title}\n"
        "HG002/HG003/HG004 smoke-test outputs only"
    )

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_major_formatter(FuncFormatter(format_bp))

    for bar, value in zip(bars, values):
        ax.annotate(
            f"{value:,}",
            xy=(bar.get_width(), bar.get_y() + bar.get_height() / 2),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=8,
        )

    fig.text(
        0.01,
        0.01,
        (
            "Technical smoke-test results. These outputs are not valid "
            "chromosome-level biological performance estimates."
        ),
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, filename)


def plot_human_status(rows: list[dict[str, object]]) -> None:
    combinations = [
        ("HG002", "ont"),
        ("HG002", "pb"),
        ("HG003", "ont"),
        ("HG003", "pb"),
        ("HG004", "ont"),
        ("HG004", "pb"),
    ]

    tools = ["Flye", "GoldRush", "ntLink", "Verkko"]

    status_rows = []

    for sample, technology in combinations:
        for tool in tools:
            available = any(
                row["sample"] == sample
                and row["technology"] == technology
                and row["tool"] == tool
                for row in rows
            )

            status_rows.append(
                {
                    "label": f"{sample} {technology.upper()} | {tool}",
                    "available": 1 if available else 0,
                }
            )

    fig_height = max(8, len(status_rows) * 0.32)
    fig, ax = plt.subplots(figsize=(11, fig_height))

    bars = ax.barh(
        [row["label"] for row in status_rows],
        [row["available"] for row in status_rows],
    )

    ax.set_xlim(0, 1.25)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["No human output", "Human output found"])
    ax.set_xlabel("Technical-output availability")
    ax.set_ylabel("Sample, technology and tool")
    ax.set_title(
        "Human technical-output availability\n"
        "No E. coli results included"
    )

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, row in zip(bars, status_rows):
        text = "AVAILABLE" if row["available"] else "NOT AVAILABLE"

        ax.text(
            1.03 if row["available"] else 0.03,
            bar.get_y() + bar.get_height() / 2,
            text,
            va="center",
            fontsize=8,
        )

    fig.tight_layout()
    save(fig, "human_output_availability")


def main() -> None:
    rows = discover_human_outputs()

    if not rows:
        raise SystemExit(
            "ERROR: No valid HG002/HG003/HG004 assembly.fasta outputs found."
        )

    table = TABLE_DIR / "human_technical_assembly_metrics.tsv"

    fields = [
        "tool",
        "sample",
        "technology",
        "dataset_type",
        "assembly_path",
        "file_size_bytes",
        "contig_count",
        "total_length_bp",
        "longest_contig_bp",
        "n50_bp",
    ]

    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Created {table.relative_to(ROOT)}")
    print(f"Human outputs included: {len(rows)}")

    plot_human_status(rows)

    plot_metric(
        rows,
        "total_length_bp",
        "Total assembled sequence length",
        "Total length (bp)",
        "human_total_length",
        log_scale=True,
    )

    plot_metric(
        rows,
        "contig_count",
        "Number of assembled sequences",
        "Contig or scaffold count",
        "human_contig_count",
    )

    plot_metric(
        rows,
        "n50_bp",
        "Assembly N50",
        "N50 (bp)",
        "human_n50",
        log_scale=True,
    )

    plot_metric(
        rows,
        "longest_contig_bp",
        "Longest assembled sequence",
        "Longest sequence (bp)",
        "human_longest_sequence",
        log_scale=True,
    )


if __name__ == "__main__":
    main()
