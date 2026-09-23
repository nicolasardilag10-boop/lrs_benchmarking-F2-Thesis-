#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import time
from pathlib import Path


ROOT = Path("benchmark_chr21_real")
MANIFEST = ROOT / "config" / "sources.resolved.tsv"
START_FILE = (
    ROOT
    / "results"
    / "provenance"
    / "current_input_run_started.epoch"
)

WINDOW_SIZE = 1_000_000


def format_seconds(seconds: float) -> str:
    seconds = max(int(seconds), 0)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def memory_information() -> tuple[float, float]:
    values = {}

    with Path("/proc/meminfo").open() as handle:
        for line in handle:
            key, value = line.split(":", 1)
            values[key] = int(value.strip().split()[0])

    total = values["MemTotal"] / 1024 / 1024
    available = values["MemAvailable"] / 1024 / 1024
    used = total - available

    return used, total


def read_single_tsv(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        return None

    with path.open(encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--interval",
        type=int,
        default=15,
    )

    parser.add_argument(
        "--once",
        action="store_true",
    )

    args = parser.parse_args()

    with MANIFEST.open(encoding="utf-8") as handle:
        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    datasets = [
        (
            row["sample"],
            row["technology"],
            int(row["chromosome_length"]),
        )
        for row in rows
    ]

    total_windows = sum(
        math.ceil(length / WINDOW_SIZE)
        for _, _, length in datasets
    )

    while True:
        now = time.time()

        if START_FILE.is_file():
            start_time = float(
                START_FILE.read_text().strip()
            )
        else:
            start_time = now

        elapsed = now - start_time
        completed_windows = 0
        extracted_count = 0
        normalized_count = 0
        dataset_lines = []

        for sample, technology, chromosome_length in datasets:
            expected_windows = math.ceil(
                chromosome_length / WINDOW_SIZE
            )

            chunk_directory = (
                ROOT
                / "work"
                / "final_inputs"
                / sample
                / technology
                / "chunks"
            )

            completed = sum(
                1
                for window_number
                in range(1, expected_windows + 1)
                if (
                    chunk_directory
                    / f"window_{window_number:03d}.bam"
                ).is_file()
            )

            completed_windows += completed

            extraction_checkpoint = (
                ROOT
                / "results"
                / "checkpoints"
                / "extraction"
                / (
                    f"{sample}.{technology}."
                    "chr21.primary.mapq20.ok"
                )
            )

            normalization_checkpoint = (
                ROOT
                / "results"
                / "checkpoints"
                / "normalization"
                / (
                    f"{sample}.{technology}."
                    "chr21.primary.mapq20.30x.ok"
                )
            )

            extraction_summary = (
                ROOT
                / "results"
                / "validation"
                / "final_extracted"
                / (
                    f"{sample}.{technology}."
                    "chr21.primary.mapq20.tsv"
                )
            )

            normalization_summary = (
                ROOT
                / "results"
                / "validation"
                / "final_normalized"
                / (
                    f"{sample}.{technology}."
                    "chr21.primary.mapq20.30x.tsv"
                )
            )

            extracted = extraction_checkpoint.is_file()
            normalized = normalization_checkpoint.is_file()

            extracted_count += int(extracted)
            normalized_count += int(normalized)

            stage = "EXTRACTING"

            if completed == expected_windows:
                stage = "MERGING/VALIDATING"

            if extracted:
                stage = "EXTRACTED PASS"

            if normalized:
                stage = "30X PASS"

            details = ""

            extraction_row = read_single_tsv(
                extraction_summary
            )

            if extraction_row:
                details += (
                    f" | source={float(extraction_row['sequence_coverage']):.3f}x"
                    f" | N50={int(extraction_row['read_n50']):,}"
                )

            normalization_row = read_single_tsv(
                normalization_summary
            )

            if normalization_row:
                details += (
                    f" | final="
                    f"{float(normalization_row['observed_sequence_coverage']):.3f}x"
                )

            dataset_lines.append(
                (
                    f"{sample} {technology.upper():<4} "
                    f"| windows {completed:02d}/{expected_windows:02d} "
                    f"| {stage}{details}"
                )
            )

        # Weighted progress:
        # 85% extraction windows, 5% merged validation, 10% normalization.
        window_progress = (
            completed_windows / total_windows
            if total_windows
            else 0
        )

        extraction_progress = extracted_count / len(datasets)
        normalization_progress = normalized_count / len(datasets)

        overall_fraction = (
            0.85 * window_progress
            + 0.05 * extraction_progress
            + 0.10 * normalization_progress
        )

        overall_percent = 100 * overall_fraction

        if overall_fraction > 0:
            estimated_total = elapsed / overall_fraction
            remaining = estimated_total - elapsed
        else:
            remaining = 0

        used_memory, total_memory = memory_information()
        disk = shutil.disk_usage(ROOT)
        disk_free_gib = disk.free / 1024**3
        load_1, load_5, load_15 = os.getloadavg()

        active_jobs = subprocess.run(
            [
                "bash",
                "-c",
                (
                    "pgrep -fc "
                    "'[s]nakemake.*Snakefile.inputs' || true"
                ),
            ],
            capture_output=True,
            text=True,
        ).stdout.strip()

        print("\033[2J\033[H", end="")
        print("=" * 78)
        print("MATCHED CHROMOSOME-21 INPUT WORKFLOW")
        print("=" * 78)
        print(f"Overall progress:       {overall_percent:6.2f}%")
        print(
            f"Extraction windows:    "
            f"{completed_windows}/{total_windows}"
        )
        print(
            f"Extracted datasets:    "
            f"{extracted_count}/{len(datasets)}"
        )
        print(
            f"Normalized datasets:   "
            f"{normalized_count}/{len(datasets)}"
        )
        print(f"Elapsed:                {format_seconds(elapsed)}")

        if overall_fraction > 0 and normalized_count < len(datasets):
            print(
                f"Estimated remaining:    "
                f"{format_seconds(remaining)}"
            )
        else:
            print("Estimated remaining:    calculating")

        print(
            f"RAM used:               "
            f"{used_memory:.2f}/{total_memory:.2f} GiB"
        )
        print(f"Disk free:              {disk_free_gib:.1f} GiB")
        print(
            f"Load average:           "
            f"{load_1:.2f}, {load_5:.2f}, {load_15:.2f}"
        )
        print(f"Snakemake processes:    {active_jobs}")
        print()
        print("-" * 78)

        for line in dataset_lines:
            print(line)

        print("-" * 78)
        print(
            f"Updated: "
            f"{time.strftime('%Y-%m-%d %H:%M:%S')}"
        )

        final_checkpoint = (
            ROOT
            / "results"
            / "checkpoints"
            / "final_inputs_ready.ok"
        )

        if final_checkpoint.is_file():
            print()
            print("FINAL STATUS: PASS")
            break

        if args.once:
            break

        time.sleep(args.interval)


if __name__ == "__main__":
    main()
