#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


EXCLUDED_FLAGS = 2308  # unmapped + secondary + supplementary


def progress(percent: int, message: str) -> None:
    print(f"[{percent:3d}%] {message}", flush=True)


def fail(message: str) -> None:
    raise RuntimeError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def require_program(program: str) -> None:
    if shutil.which(program) is None:
        fail(f"Required program was not found: {program}")


def check_output_absent(paths: list[Path]) -> None:
    existing = [str(path) for path in paths if path.exists()]

    if existing:
        fail(
            "Refusing to overwrite existing output files:\n"
            + "\n".join(existing)
        )


def reference_targets(bam: Path) -> dict[str, int]:
    completed = subprocess.run(
        ["samtools", "view", "-H", str(bam)],
        check=True,
        capture_output=True,
        text=True,
    )

    targets: dict[str, int] = {}

    for line in completed.stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        name: str | None = None
        reference_length: int | None = None

        for field in line.split("\t"):
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                reference_length = int(field[3:])

        if name is not None and reference_length is not None:
            targets[name] = reference_length

    return targets


def collect_reads(
    bam: Path,
    chromosome: str,
    seed: int,
) -> list[tuple[bytes, str, int]]:
    command = [
        "samtools",
        "view",
        "-F",
        str(EXCLUDED_FLAGS),
        str(bam),
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        fail("Could not read samtools output.")

    seen_names: set[str] = set()
    ranked_reads: list[tuple[bytes, str, int]] = []

    for line_number, line in enumerate(process.stdout, start=1):
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            fail(
                f"Malformed SAM record at output line {line_number}."
            )

        read_name = fields[0]
        reference_name = fields[2]
        sequence = fields[9]

        if reference_name != chromosome:
            fail(
                f"Unexpected reference {reference_name!r} for "
                f"read {read_name!r}; expected {chromosome!r}."
            )

        if sequence == "*":
            fail(
                f"Read {read_name!r} has no sequence in the BAM."
            )

        if read_name in seen_names:
            fail(
                f"Duplicate primary read name found: {read_name}"
            )

        seen_names.add(read_name)

        rank = hashlib.sha256(
            f"{seed}\0{read_name}".encode("utf-8")
        ).digest()

        ranked_reads.append(
            (
                rank,
                read_name,
                len(sequence),
            )
        )

    stderr = process.stderr.read() if process.stderr else ""
    return_code = process.wait()

    if return_code != 0:
        fail(
            "samtools view failed while reading the source BAM:\n"
            + stderr
        )

    if not ranked_reads:
        fail("No eligible reads were found in the source BAM.")

    return ranked_reads


def select_nearest_target(
    ranked_reads: list[tuple[bytes, str, int]],
    target_bases: int,
) -> tuple[list[str], int, int]:
    ranked_reads.sort(key=lambda record: record[0])

    selected_names: list[str] = []
    selected_bases = 0
    maximum_read_length = 0

    for _, read_name, read_length in ranked_reads:
        maximum_read_length = max(
            maximum_read_length,
            read_length,
        )

        if selected_bases >= target_bases:
            break

        candidate_bases = selected_bases + read_length

        if candidate_bases < target_bases:
            selected_names.append(read_name)
            selected_bases = candidate_bases
            continue

        distance_without = abs(target_bases - selected_bases)
        distance_with = abs(candidate_bases - target_bases)

        if distance_with <= distance_without:
            selected_names.append(read_name)
            selected_bases = candidate_bases

        break

    if not selected_names:
        fail("The normalization algorithm selected no reads.")

    return (
        selected_names,
        selected_bases,
        maximum_read_length,
    )


def validate_fastq(
    fastq: Path,
) -> tuple[int, int, set[str]]:
    read_count = 0
    total_bases = 0
    names: set[str] = set()

    with gzip.open(
        fastq,
        "rt",
        encoding="utf-8",
    ) as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")

            if not sequence or not plus or not quality:
                fail(f"Truncated FASTQ file: {fastq}")

            if not header.startswith("@"):
                fail(f"Invalid FASTQ header: {header.rstrip()}")

            if not plus.startswith("+"):
                fail(
                    f"Invalid FASTQ separator for "
                    f"{header.rstrip()}"
                )

            if len(sequence) != len(quality):
                fail(
                    f"Sequence/quality length mismatch for "
                    f"{header.rstrip()}"
                )

            read_name = header[1:].strip().split()[0]

            if read_name in names:
                fail(
                    f"Duplicate read name in normalized FASTQ: "
                    f"{read_name}"
                )

            names.add(read_name)
            read_count += 1
            total_bases += len(sequence)

    return read_count, total_bases, names


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Deterministically normalize a chromosome-specific "
            "long-read BAM to a target sequence coverage."
        )
    )

    parser.add_argument("--input-bam", required=True, type=Path)
    parser.add_argument("--output-fastq", required=True, type=Path)
    parser.add_argument("--selected-names", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--log", required=True, type=Path)

    parser.add_argument("--sample", required=True)
    parser.add_argument("--technology", required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument(
        "--chromosome-length",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--target-coverage",
        required=True,
        type=float,
    )
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument(
        "--selection-policy",
        default="primary_MAPQ20_source",
    )

    args = parser.parse_args()

    require_program("samtools")
    require_program("gzip")

    input_bam = args.input_bam.resolve()
    output_fastq = args.output_fastq.resolve()
    selected_names = args.selected_names.resolve()
    summary = args.summary.resolve()
    checkpoint = args.checkpoint.resolve()
    log_path = args.log.resolve()

    temporary_fastq = Path(str(output_fastq) + ".partial")
    temporary_names = Path(str(selected_names) + ".partial")

    check_output_absent(
        [
            output_fastq,
            selected_names,
            summary,
            checkpoint,
            temporary_fastq,
            temporary_names,
        ]
    )

    for path in [
        output_fastq,
        selected_names,
        summary,
        checkpoint,
        log_path,
    ]:
        path.parent.mkdir(parents=True, exist_ok=True)

    try:
        progress(
            0,
            (
                f"Starting normalization for "
                f"{args.sample} {args.technology}."
            ),
        )

        if not input_bam.is_file() or input_bam.stat().st_size == 0:
            fail(f"Input BAM is missing or empty: {input_bam}")

        progress(5, "Checking BAM integrity.")

        subprocess.run(
            ["samtools", "quickcheck", str(input_bam)],
            check=True,
        )

        targets = reference_targets(input_bam)

        if args.chromosome not in targets:
            fail(
                f"{args.chromosome} was not found in the BAM header."
            )

        observed_reference_length = targets[args.chromosome]

        if observed_reference_length != args.chromosome_length:
            fail(
                "Chromosome length mismatch: "
                f"BAM={observed_reference_length}, "
                f"expected={args.chromosome_length}"
            )

        target_bases = round(
            args.chromosome_length
            * args.target_coverage
        )

        progress(
            15,
            (
                f"Scanning eligible reads. "
                f"Target: {target_bases:,} bases."
            ),
        )

        ranked_reads = collect_reads(
            input_bam,
            args.chromosome,
            args.seed,
        )

        source_bases = sum(
            read_length
            for _, _, read_length in ranked_reads
        )

        source_coverage = (
            source_bases
            / args.chromosome_length
        )

        if source_bases < target_bases:
            fail(
                "Source coverage is below the requested target: "
                f"source={source_coverage:.6f}x, "
                f"target={args.target_coverage:.6f}x"
            )

        progress(
            40,
            (
                f"Source contains {len(ranked_reads):,} unique reads "
                f"and {source_coverage:.6f}x sequence coverage."
            ),
        )

        (
            chosen_names,
            predicted_bases,
            maximum_read_length,
        ) = select_nearest_target(
            ranked_reads,
            target_bases,
        )

        with temporary_names.open(
            "w",
            encoding="utf-8",
        ) as handle:
            for read_name in chosen_names:
                handle.write(read_name + "\n")

        progress(
            55,
            (
                f"Selected {len(chosen_names):,} whole reads "
                f"using seed {args.seed}."
            ),
        )

        quoted_input = shlex.quote(str(input_bam))
        quoted_names = shlex.quote(str(temporary_names))
        quoted_fastq = shlex.quote(str(temporary_fastq))
        quoted_log = shlex.quote(str(log_path))

        shell_command = (
            "set -o pipefail; "
            f"samtools view -F {EXCLUDED_FLAGS} "
            f"-N {quoted_names} -u {quoted_input} "
            f"2> {quoted_log} "
            "| samtools fastq -n - "
            f"2>> {quoted_log} "
            "| gzip -1 "
            f"> {quoted_fastq}"
        )

        progress(
            65,
            "Writing deterministic normalized FASTQ.",
        )

        subprocess.run(
            ["bash", "-c", shell_command],
            check=True,
        )

        subprocess.run(
            ["gzip", "-t", str(temporary_fastq)],
            check=True,
        )

        progress(85, "Validating normalized FASTQ.")

        (
            observed_reads,
            observed_bases,
            observed_names,
        ) = validate_fastq(temporary_fastq)

        expected_names = set(chosen_names)

        if observed_names != expected_names:
            missing = expected_names - observed_names
            foreign = observed_names - expected_names

            fail(
                "Normalized FASTQ read-name mismatch. "
                f"Missing={len(missing)}, foreign={len(foreign)}"
            )

        if observed_reads != len(chosen_names):
            fail(
                "Normalized FASTQ read count differs from "
                "the selected-name count."
            )

        if observed_bases != predicted_bases:
            fail(
                "Normalized FASTQ base count differs from "
                f"predicted bases: FASTQ={observed_bases}, "
                f"predicted={predicted_bases}"
            )

        observed_coverage = (
            observed_bases
            / args.chromosome_length
        )

        maximum_expected_difference = max(
            maximum_read_length
            / args.chromosome_length,
            0.001,
        )

        coverage_difference = abs(
            observed_coverage
            - args.target_coverage
        )

        if coverage_difference > maximum_expected_difference:
            fail(
                "Observed coverage is outside the expected "
                "whole-read tolerance: "
                f"observed={observed_coverage:.6f}x, "
                f"target={args.target_coverage:.6f}x"
            )

        os.replace(temporary_fastq, output_fastq)
        os.replace(temporary_names, selected_names)

        fastq_sha256 = sha256_file(output_fastq)
        names_sha256 = sha256_file(selected_names)
        source_bam_sha256 = sha256_file(input_bam)

        fields = [
            "sample",
            "technology",
            "chromosome",
            "chromosome_length",
            "selection_policy",
            "seed",
            "source_reads",
            "source_read_bases",
            "source_sequence_coverage",
            "target_coverage",
            "target_bases",
            "selected_reads",
            "selected_read_bases",
            "observed_sequence_coverage",
            "coverage_difference",
            "input_bam",
            "output_fastq",
            "selected_names",
            "input_bam_sha256",
            "output_fastq_sha256",
            "selected_names_sha256",
            "validation_status",
        ]

        values = {
            "sample": args.sample,
            "technology": args.technology,
            "chromosome": args.chromosome,
            "chromosome_length": args.chromosome_length,
            "selection_policy": args.selection_policy,
            "seed": args.seed,
            "source_reads": len(ranked_reads),
            "source_read_bases": source_bases,
            "source_sequence_coverage": f"{source_coverage:.6f}",
            "target_coverage": f"{args.target_coverage:.6f}",
            "target_bases": target_bases,
            "selected_reads": observed_reads,
            "selected_read_bases": observed_bases,
            "observed_sequence_coverage": f"{observed_coverage:.6f}",
            "coverage_difference": f"{coverage_difference:.6f}",
            "input_bam": str(input_bam),
            "output_fastq": str(output_fastq),
            "selected_names": str(selected_names),
            "input_bam_sha256": source_bam_sha256,
            "output_fastq_sha256": fastq_sha256,
            "selected_names_sha256": names_sha256,
            "validation_status": "PASS",
        }

        with summary.open(
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

        checkpoint.write_text(
            "\n".join(
                [
                    "NORMALIZED INPUT READY",
                    f"Sample: {args.sample}",
                    f"Technology: {args.technology}",
                    f"Chromosome: {args.chromosome}",
                    f"Seed: {args.seed}",
                    (
                        "Observed coverage: "
                        f"{observed_coverage:.6f}x"
                    ),
                    f"FASTQ: {output_fastq}",
                    f"SHA256: {fastq_sha256}",
                    "STATUS: PASS",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        progress(
            100,
            (
                f"Normalization completed: "
                f"{observed_reads:,} reads, "
                f"{observed_coverage:.6f}x coverage, PASS."
            ),
        )

        print()
        print(f"FASTQ:      {output_fastq}")
        print(f"Names:      {selected_names}")
        print(f"Summary:    {summary}")
        print(f"Checkpoint: {checkpoint}")

    except Exception:
        temporary_fastq.unlink(missing_ok=True)
        temporary_names.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        sys.exit(1)
