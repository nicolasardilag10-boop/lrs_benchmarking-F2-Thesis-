#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
from pathlib import Path


EXCLUDED_FLAGS = 2308


def fail(message: str) -> None:
    raise RuntimeError(message)


def n50(lengths: list[int]) -> int:
    if not lengths:
        return 0

    target = sum(lengths) / 2
    cumulative = 0

    for read_length in sorted(lengths, reverse=True):
        cumulative += read_length

        if cumulative >= target:
            return read_length

    return 0


def reference_targets(bam: Path) -> dict[str, int]:
    result = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )

    targets: dict[str, int] = {}

    for line in result.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        name: str | None = None
        target_length: int | None = None

        for field in line.split("\t"):
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                target_length = int(field[3:])

        if name is not None and target_length is not None:
            targets[name] = target_length

    return targets


def mapped_cigar_bases(bam: Path) -> int:
    result = subprocess.run(
        ["samtools", "stats", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )

    for line in result.stdout.splitlines():
        fields = line.split("\t")

        if (
            len(fields) >= 3
            and fields[0] == "SN"
            and fields[1] == "bases mapped (cigar):"
        ):
            return int(fields[2].split()[0])

    return 0


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--technology", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--chromosome-length", required=True, type=int)
    parser.add_argument("--minimum-mapq", required=True, type=int)
    parser.add_argument("--target-coverage", required=True, type=float)

    args = parser.parse_args()

    bam = args.bam.resolve()
    summary = args.summary.resolve()
    checkpoint = args.checkpoint.resolve()

    if not bam.is_file() or bam.stat().st_size == 0:
        fail(f"BAM is missing or empty: {bam}")

    subprocess.run(
        ["samtools", "quickcheck", str(bam)],
        check=True,
    )

    subprocess.run(
        ["samtools", "idxstats", str(bam)],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    targets = reference_targets(bam)

    if args.chromosome not in targets:
        fail(
            f"{args.chromosome} was not found in the BAM header."
        )

    observed_length = targets[args.chromosome]

    if observed_length != args.chromosome_length:
        fail(
            "Chromosome-length mismatch: "
            f"observed={observed_length}, "
            f"expected={args.chromosome_length}"
        )

    process = subprocess.Popen(
        ["samtools", "view", str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        fail("Could not read the BAM.")

    read_names: set[str] = set()
    read_lengths: list[int] = []

    for line_number, line in enumerate(process.stdout, start=1):
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            fail(f"Malformed SAM record at line {line_number}")

        read_name = fields[0]
        flag = int(fields[1])
        reference_name = fields[2]
        mapq = int(fields[4])
        sequence = fields[9]

        if flag & EXCLUDED_FLAGS:
            fail(
                f"Excluded alignment flag found for read {read_name}"
            )

        if reference_name != args.chromosome:
            fail(
                f"Unexpected target {reference_name} for {read_name}"
            )

        if mapq < args.minimum_mapq:
            fail(
                f"MAPQ below threshold for read {read_name}: {mapq}"
            )

        if sequence == "*":
            fail(f"Read sequence is absent for {read_name}")

        if read_name in read_names:
            fail(f"Duplicate primary read name: {read_name}")

        read_names.add(read_name)
        read_lengths.append(len(sequence))

    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()

    if return_code != 0:
        fail(f"samtools view failed:\n{stderr}")

    if not read_lengths:
        fail("The merged extracted BAM contains zero reads.")

    total_read_bases = sum(read_lengths)
    sequence_coverage = total_read_bases / args.chromosome_length
    cigar_bases = mapped_cigar_bases(bam)
    cigar_coverage = cigar_bases / args.chromosome_length

    if sequence_coverage < args.target_coverage:
        fail(
            "Source coverage is insufficient for normalization: "
            f"observed={sequence_coverage:.6f}x, "
            f"required={args.target_coverage:.6f}x"
        )

    summary.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "sample",
        "technology",
        "chromosome",
        "chromosome_length",
        "minimum_mapq",
        "reads",
        "total_read_bases",
        "mean_read_length",
        "median_read_length",
        "read_n50",
        "maximum_read_length",
        "sequence_coverage",
        "mapped_cigar_bases",
        "mapped_cigar_coverage",
        "target_coverage",
        "bam",
        "validation_status",
    ]

    values = {
        "sample": args.sample,
        "technology": args.technology,
        "chromosome": args.chromosome,
        "chromosome_length": args.chromosome_length,
        "minimum_mapq": args.minimum_mapq,
        "reads": len(read_lengths),
        "total_read_bases": total_read_bases,
        "mean_read_length": f"{statistics.mean(read_lengths):.3f}",
        "median_read_length": f"{statistics.median(read_lengths):.3f}",
        "read_n50": n50(read_lengths),
        "maximum_read_length": max(read_lengths),
        "sequence_coverage": f"{sequence_coverage:.6f}",
        "mapped_cigar_bases": cigar_bases,
        "mapped_cigar_coverage": f"{cigar_coverage:.6f}",
        "target_coverage": f"{args.target_coverage:.6f}",
        "bam": str(bam),
        "validation_status": "PASS",
    }

    temporary_summary = Path(str(summary) + ".partial")

    with temporary_summary.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(values)

    temporary_summary.replace(summary)

    checkpoint.write_text(
        "\n".join(
            [
                "EXTRACTED CHROMOSOME-21 BAM READY",
                f"Sample: {args.sample}",
                f"Technology: {args.technology}",
                f"Reads: {len(read_lengths)}",
                f"Sequence coverage: {sequence_coverage:.6f}x",
                f"Read N50: {n50(read_lengths)}",
                f"BAM: {bam}",
                "STATUS: PASS",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(
        (
            f"PASS: {args.sample} {args.technology}; "
            f"reads={len(read_lengths):,}; "
            f"coverage={sequence_coverage:.6f}x; "
            f"N50={n50(read_lengths):,}"
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
