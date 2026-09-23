#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SAMPLES_TSV = ROOT / "samples.tsv"
PROGRESS = ROOT / "results" / "progress"


def status(
    done: Path,
    running: Path | None = None,
    failed: Path | None = None,
) -> str:
    if done.exists():
        return "DONE"

    if failed is not None and failed.exists():
        return "FAILED"

    if running is not None and running.exists():
        return "RUNNING"

    return "WAIT"


def chr21_status(sample: str, tech: str, stage: str) -> str:
    base = PROGRESS / "chr21"

    done = base / f"{sample}.{tech}.{stage}.done"

    running = None
    failed = None

    if stage == "mapped":
        running = base / f"{sample}.{tech}.mapping.running"

    return status(
        done=done,
        running=running,
        failed=failed,
    )


def assembler_status(
    tool: str,
    sample: str,
    tech: str | None = None,
) -> str:
    base = PROGRESS / tool

    if tech is None:
        stem = sample
    else:
        stem = f"{sample}.{tech}"

    return status(
        done=base / f"{stem}.done",
        running=base / f"{stem}.running",
        failed=base / f"{stem}.failed",
    )


def read_samples() -> list[tuple[str, str]]:
    rows = []

    with SAMPLES_TSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            rows.append(
                (
                    row["sample"].strip(),
                    row["technology"].strip(),
                )
            )

    return rows


def detect_mode() -> str:
    logs = ROOT / "results" / "logs"

    if (
        logs
        / "validate_inputs.whole_genome.ok"
    ).exists():
        return "whole_genome"

    if (
        logs
        / "validate_inputs.preextracted.ok"
    ).exists():
        return "preextracted"

    return "unknown"


def main() -> None:
    samples = read_samples()
    mode = detect_mode()

    headers = [
        "SAMPLE",
        "TECH",
        "MAP",
        "FILTER",
        "30X",
        "FLYE",
        "GOLDRUSH",
        "NTLINK",
        "VERKKO",
    ]

    table = []

    for sample, tech in samples:
        if mode == "preextracted":
            mapping = "N/A"
            filtering = "N/A"
            normalization = "PRESET"
        else:
            mapping = chr21_status(
                sample,
                tech,
                "mapped",
            )

            filtering = chr21_status(
                sample,
                tech,
                "filtered",
            )

            normalization = chr21_status(
                sample,
                tech,
                "normalized",
            )

        table.append(
            [
                sample,
                tech,
                mapping,
                filtering,
                normalization,
                assembler_status(
                    "flye",
                    sample,
                    tech,
                ),
                assembler_status(
                    "goldrush",
                    sample,
                    tech,
                ),
                assembler_status(
                    "ntlink",
                    sample,
                    tech,
                ),
                assembler_status(
                    "verkko",
                    sample,
                ),
            ]
        )

    widths = []

    for index, header in enumerate(headers):
        widths.append(
            max(
                len(header),
                *[
                    len(str(row[index]))
                    for row in table
                ],
            )
        )

    def format_row(row):
        return "  ".join(
            str(value).ljust(widths[index])
            for index, value in enumerate(row)
        )

    print()
    print("ASSEMBLER WORKFLOW STATUS")
    print("=" * 90)
    print(f"Mode: {mode}")
    print(f"Project: {ROOT}")
    print()

    print(format_row(headers))

    print(
        format_row(
            [
                "-" * width
                for width in widths
            ]
        )
    )

    for row in table:
        print(format_row(row))

    print()
    print("Status meanings:")
    print("  DONE    validated/completed output checkpoint exists")
    print("  RUNNING active-stage marker exists")
    print("  FAILED  failure marker exists")
    print("  WAIT    stage has not completed")
    print("  PRESET  pre-extracted normalized Chr21 input is being used")
    print("  N/A     stage is unnecessary in pre-extracted mode")
    print()


if __name__ == "__main__":
    main()
