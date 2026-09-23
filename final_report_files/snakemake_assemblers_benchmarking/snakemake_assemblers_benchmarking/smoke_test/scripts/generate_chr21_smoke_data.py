#!/usr/bin/env python3

from pathlib import Path
import gzip
import sys


PROJECT = Path.cwd()

OUTPUT_REFERENCE = (
    PROJECT
    / "smoke_test"
    / "data"
    / "reference"
    / "chr21_1mb.fasta"
)

OUTPUT_READS = (
    PROJECT
    / "smoke_test"
    / "data"
    / "reads_chr21"
)

SAMPLES = ["HG002", "HG003", "HG004"]

REGION_LENGTH = 1_000_000
WINDOW_STEP = 100_000
MIN_ACGT_FRACTION = 0.98


def find_chr21_reference() -> Path:
    candidates = []

    for root in [PROJECT / "data", PROJECT / "smoke_test"]:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            name = path.name.lower()

            if "chr21" not in name:
                continue

            if name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
                if path.resolve() != OUTPUT_REFERENCE.resolve():
                    candidates.append(path)

    if not candidates:
        raise FileNotFoundError(
            "No chromosome 21 FASTA was found under data/ or smoke_test/."
        )

    candidates.sort(key=lambda path: len(str(path)))
    return candidates[0]


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")

    return path.open("r")


def read_fasta(path: Path) -> str:
    sequence_parts = []

    with open_text(path) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith(">"):
                continue

            sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)

    if not sequence:
        raise ValueError(f"No sequence found in {path}")

    return sequence


def acgt_fraction(sequence: str) -> float:
    acgt_count = sum(base in "ACGT" for base in sequence)
    return acgt_count / len(sequence)


def choose_high_quality_region(sequence: str) -> tuple[int, str, float]:
    if len(sequence) < REGION_LENGTH:
        raise ValueError(
            f"Reference has only {len(sequence):,} bp; "
            f"{REGION_LENGTH:,} bp are required."
        )

    best_start = None
    best_region = None
    best_fraction = -1.0

    for start in range(
        0,
        len(sequence) - REGION_LENGTH + 1,
        WINDOW_STEP,
    ):
        region = sequence[start : start + REGION_LENGTH]
        fraction = acgt_fraction(region)

        if fraction > best_fraction:
            best_start = start
            best_region = region
            best_fraction = fraction

        if fraction >= MIN_ACGT_FRACTION:
            break

    if best_region is None:
        raise RuntimeError("Could not select a chromosome 21 region.")

    if best_fraction < 0.90:
        raise RuntimeError(
            f"Best 1 Mb region contains only "
            f"{best_fraction:.2%} A/C/G/T bases."
        )

    return best_start, best_region, best_fraction


def write_reference(
    region: str,
    start: int,
    fraction: float,
) -> None:
    OUTPUT_REFERENCE.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_REFERENCE.open("w") as handle:
        handle.write(
            f">chr21_smoke_1mb_"
            f"{start + 1}_{start + len(region)}_"
            f"acgt_{fraction:.4f}\n"
        )

        for position in range(0, len(region), 80):
            handle.write(region[position : position + 80] + "\n")


def write_fastq(
    region: str,
    sample: str,
    technology: str,
    read_length: int,
    step: int,
) -> None:
    output = (
        OUTPUT_READS
        / f"{sample}.{technology}.chr21.smoke.fastq"
    )

    output.parent.mkdir(parents=True, exist_ok=True)

    read_count = 0
    total_bases = 0

    with output.open("w") as handle:
        for start in range(
            0,
            len(region) - read_length + 1,
            step,
        ):
            read = region[start : start + read_length]

            # Skip reads containing too many ambiguous bases.
            if acgt_fraction(read) < 0.98:
                continue

            read_count += 1
            total_bases += len(read)

            read_name = (
                f"{sample}_{technology}_"
                f"{read_count:06d}_"
                f"{start + 1}_{start + read_length}"
            )

            handle.write(f"@{read_name}\n")
            handle.write(f"{read}\n")
            handle.write("+\n")
            handle.write("I" * len(read) + "\n")

    if read_count == 0:
        raise RuntimeError(
            f"No valid reads generated for {sample} {technology}."
        )

    coverage = total_bases / len(region)

    print(
        f"{sample} {technology}: "
        f"{read_count:,} reads, "
        f"{total_bases:,} bases, "
        f"{coverage:.2f}x coverage"
    )


def main() -> None:
    reference = find_chr21_reference()
    sequence = read_fasta(reference)

    start, region, fraction = choose_high_quality_region(sequence)

    print(f"Source reference: {reference}")
    print(
        f"Selected region: "
        f"{start + 1:,}-{start + len(region):,}"
    )
    print(f"A/C/G/T fraction: {fraction:.4%}")

    write_reference(region, start, fraction)

    for sample in SAMPLES:
        write_fastq(
            region,
            sample,
            "ont",
            read_length=20_000,
            step=2_000,
        )

        write_fastq(
            region,
            sample,
            "pb",
            read_length=15_000,
            step=1_500,
        )

    print()
    print("Smoke data regenerated successfully.")
    print(f"Reference: {OUTPUT_REFERENCE}")
    print(f"Reads: {OUTPUT_READS}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
