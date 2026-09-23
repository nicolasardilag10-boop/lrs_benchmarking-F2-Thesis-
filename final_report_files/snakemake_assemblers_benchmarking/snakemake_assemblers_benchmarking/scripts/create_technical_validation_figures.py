#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
TABLE_DIR = RESULTS_DIR / "tables"
FIGURE_DIR = RESULTS_DIR / "figures"
SAMPLE_SHEET = ROOT / "samples.tsv"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def read_sample_sheet() -> dict[tuple[str, str], str]:
    """Return {(sample, technology): FASTQ path} from samples.tsv."""
    sample_inputs: dict[tuple[str, str], str] = {}

    if not SAMPLE_SHEET.exists():
        return sample_inputs

    with SAMPLE_SHEET.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"sample", "technology", "fastq"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                f"{SAMPLE_SHEET} must contain: sample, technology, fastq"
            )

        for row in reader:
            sample = row["sample"].strip()
            technology = row["technology"].strip().lower()
            fastq = row["fastq"].strip()

            sample_inputs[(sample, technology)] = fastq

    return sample_inputs


def calculate_n50(lengths: list[int]) -> int:
    """Calculate N50 from a list of sequence lengths."""
    if not lengths:
        return 0

    total_length = sum(lengths)
    threshold = total_length / 2
    cumulative = 0

    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length

    return 0


def read_fasta_metrics(path: Path) -> dict[str, int]:
    """Calculate basic assembly metrics without external dependencies."""
    sequence_lengths: list[int] = []
    current_length = 0
    found_header = False

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if found_header:
                    sequence_lengths.append(current_length)

                found_header = True
                current_length = 0
            else:
                if not found_header:
                    raise ValueError(
                        f"{path} does not appear to be a valid FASTA file"
                    )

                current_length += len(line)

    if found_header:
        sequence_lengths.append(current_length)

    sequence_lengths = [
        length for length in sequence_lengths if length > 0
    ]

    if not sequence_lengths:
        raise ValueError(f"No non-empty FASTA sequences found in {path}")

    return {
        "contig_count": len(sequence_lengths),
        "total_length_bp": sum(sequence_lengths),
        "longest_contig_bp": max(sequence_lengths),
        "n50_bp": calculate_n50(sequence_lengths),
    }


def infer_assembler(path: Path) -> str | None:
    """Infer assembler from the final assembly output path."""
    path_string = path.as_posix().lower()

    if "/ntlink/" in path_string:
        return "ntLink"

    if "/goldrush/" in path_string:
        return "GoldRush"

    if "/verkko/" in path_string:
        return "Verkko"

    if (
        "/flye/" in path_string
        or ".flye/" in path_string
        or path.parent.name.lower().endswith(".flye")
    ):
        return "Flye"

    return None


def infer_sample_and_technology(path: Path) -> tuple[str, str]:
    """Extract sample and technology from common project path formats."""
    path_string = path.as_posix()

    human_match = re.search(
        r"(HG00[234])[._](ont|pb)",
        path_string,
        flags=re.IGNORECASE,
    )

    if human_match:
        sample = human_match.group(1).upper()
        technology = human_match.group(2).lower()
        return sample, technology

    if "ecoli" in path_string.lower():
        return "E_coli", "HiFi+ONT"

    return "unknown", "unknown"


def infer_dataset_type(
    sample: str,
    technology: str,
    input_fastq: str,
    assembly_path: Path,
) -> str:
    """Assign an explicit technical-validation dataset category."""
    combined = f"{input_fastq} {assembly_path.as_posix()}".lower()

    if sample == "E_coli":
        return "official_ecoli_test"

    if "1k" in combined:
        return "human_1k_smoke_test"

    if "smoke" in combined:
        return "human_smoke_test"

    if sample.startswith("HG00"):
        return "human_technical_test_unknown_subset"

    return "technical_test_unknown_dataset"


def discover_final_assemblies() -> list[Path]:
    """Find only standardized final assembly.fasta files."""
    assemblies: list[Path] = []

    for path in RESULTS_DIR.rglob("assembly.fasta"):
        relative_parts = {
            part.lower() for part in path.relative_to(ROOT).parts
        }

        # Exclude internal temporary and manual-test assemblies.
        if "work" in relative_parts:
            continue

        if "manual_test" in relative_parts:
            continue

        if not path.is_file() or path.stat().st_size == 0:
            continue

        if infer_assembler(path) is None:
            continue

        assemblies.append(path)

    return sorted(set(assemblies))


def assembly_label(row: dict[str, object]) -> str:
    """Create an explicit label that includes the dataset category."""
    return (
        f"{row['assembler']} | "
        f"{row['sample']} | "
        f"{row['technology']} | "
        f"{row['dataset_type']}"
    )


def save_figure(fig: plt.Figure, basename: str) -> None:
    """Save every figure as PNG and PDF."""
    png_path = FIGURE_DIR / f"{basename}.png"
    pdf_path = FIGURE_DIR / f"{basename}.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Created: {png_path.relative_to(ROOT)}")
    print(f"Created: {pdf_path.relative_to(ROOT)}")


def human_number(value: float, _: int) -> str:
    """Format base-pair counts for plot axes."""
    if value >= 1_000_000_000:
        return f"{value / 1_000_000_000:.1f}G"

    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"

    if value >= 1_000:
        return f"{value / 1_000:.1f}K"

    return f"{value:.0f}"


def create_status_figure(rows: list[dict[str, object]]) -> None:
    """Show whether each assembler has at least one valid final output."""
    assemblers = ["Flye", "GoldRush", "Verkko", "ntLink"]

    status_values = []
    status_text = []

    for assembler in assemblers:
        passed = any(
            row["assembler"] == assembler and row["status"] == "PASS"
            for row in rows
        )

        status_values.append(1 if passed else 0)
        status_text.append("PASS" if passed else "NO VALID OUTPUT")

    fig, ax = plt.subplots(figsize=(9, 5.5))

    bars = ax.barh(assemblers, status_values)

    ax.set_xlim(0, 1.20)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Not validated", "Validated"])
    ax.set_xlabel("Technical execution status")
    ax.set_ylabel("Assembler")
    ax.set_title(
        "Assembler technical-validation status\n"
        "Smoke tests only — not biological performance ranking"
    )

    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, text in zip(bars, status_text):
        x_position = 1.03 if text == "PASS" else 0.03
        ax.text(
            x_position,
            bar.get_y() + bar.get_height() / 2,
            text,
            va="center",
            fontsize=10,
        )

    fig.tight_layout()
    save_figure(fig, "technical_validation_status")


def create_metric_figure(
    rows: list[dict[str, object]],
    metric: str,
    title: str,
    x_label: str,
    basename: str,
    use_log_scale: bool,
) -> None:
    """Create one horizontal bar chart for an assembly metric."""
    valid_rows = [
        row
        for row in rows
        if row["status"] == "PASS" and int(row[metric]) > 0
    ]

    if not valid_rows:
        print(f"Skipping {basename}: no valid data")
        return

    valid_rows.sort(
        key=lambda row: (
            str(row["dataset_type"]),
            str(row["assembler"]),
            str(row["sample"]),
            str(row["technology"]),
        )
    )

    labels = [assembly_label(row) for row in valid_rows]
    values = [int(row[metric]) for row in valid_rows]

    figure_height = max(5.5, 0.52 * len(valid_rows) + 2.5)
    fig, ax = plt.subplots(figsize=(13, figure_height))

    bars = ax.barh(labels, values)

    ax.set_xlabel(x_label)
    ax.set_ylabel("Technical-validation output")
    ax.set_title(
        f"{title}\n"
        "Mixed smoke-test datasets — do not interpret as assembler ranking"
    )

    if use_log_scale and min(values) > 0:
        ax.set_xscale("log")
        ax.set_xlabel(f"{x_label} — logarithmic scale")

    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

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
            "Dataset categories are shown in every label. "
            "These figures demonstrate technical execution only."
        ),
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.035, 1, 1))
    save_figure(fig, basename)


def main() -> int:
    sample_inputs = read_sample_sheet()
    assembly_paths = discover_final_assemblies()

    if not assembly_paths:
        print(
            "ERROR: No standardized final assembly.fasta files were found.",
            file=sys.stderr,
        )
        return 1

    rows: list[dict[str, object]] = []

    for path in assembly_paths:
        assembler = infer_assembler(path)

        if assembler is None:
            continue

        sample, technology = infer_sample_and_technology(path)

        input_fastq = sample_inputs.get(
            (sample, technology.lower()),
            "",
        )

        dataset_type = infer_dataset_type(
            sample=sample,
            technology=technology,
            input_fastq=input_fastq,
            assembly_path=path,
        )

        row: dict[str, object] = {
            "assembler": assembler,
            "sample": sample,
            "technology": technology,
            "dataset_type": dataset_type,
            "input_fastq": input_fastq,
            "assembly_path": str(path.relative_to(ROOT)),
            "file_size_bytes": path.stat().st_size,
            "contig_count": 0,
            "total_length_bp": 0,
            "longest_contig_bp": 0,
            "n50_bp": 0,
            "status": "FAIL",
        }

        try:
            metrics = read_fasta_metrics(path)
            row.update(metrics)
            row["status"] = "PASS"
        except (OSError, ValueError) as error:
            print(f"WARNING: {path}: {error}", file=sys.stderr)

        rows.append(row)

    output_table = TABLE_DIR / "assembler_technical_validation.tsv"

    fieldnames = [
        "assembler",
        "sample",
        "technology",
        "dataset_type",
        "input_fastq",
        "assembly_path",
        "file_size_bytes",
        "contig_count",
        "total_length_bp",
        "longest_contig_bp",
        "n50_bp",
        "status",
    ]

    with output_table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Created: {output_table.relative_to(ROOT)}")
    print(f"Valid final assemblies found: {len(rows)}")

    create_status_figure(rows)

    create_metric_figure(
        rows=rows,
        metric="total_length_bp",
        title="Total assembled length",
        x_label="Total assembly length (bp)",
        basename="assembly_total_length",
        use_log_scale=True,
    )

    create_metric_figure(
        rows=rows,
        metric="contig_count",
        title="Assembly sequence count",
        x_label="Number of contigs or scaffolds",
        basename="assembly_contig_count",
        use_log_scale=False,
    )

    create_metric_figure(
        rows=rows,
        metric="n50_bp",
        title="Assembly N50",
        x_label="N50 (bp)",
        basename="assembly_n50",
        use_log_scale=True,
    )

    create_metric_figure(
        rows=rows,
        metric="longest_contig_bp",
        title="Longest assembled sequence",
        x_label="Longest contig or scaffold (bp)",
        basename="assembly_longest_contig",
        use_log_scale=True,
    )

    print("\nTechnical-validation figure generation completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
