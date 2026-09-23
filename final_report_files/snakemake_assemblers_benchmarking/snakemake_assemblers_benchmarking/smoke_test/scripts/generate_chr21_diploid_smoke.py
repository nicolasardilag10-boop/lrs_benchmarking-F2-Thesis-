#!/usr/bin/env python3

from pathlib import Path
import random

PROJECT = Path.cwd()

REFERENCE = (
    PROJECT
    / "smoke_test"
    / "data"
    / "reference"
    / "chr21_1mb.fasta"
)

OUTPUT_DIR = (
    PROJECT
    / "smoke_test"
    / "data"
    / "reads_chr21"
)

SAMPLES = ["HG002", "HG003", "HG004"]

BASES = "ACGT"
SEED = 20260804

HIFI_LENGTH = 15000
HIFI_STEP = 1500
HIFI_ERROR_RATE = 0.002

ONT_LENGTH = 20000
ONT_STEP = 2000
ONT_ERROR_RATE = 0.03

HETEROZYGOUS_RATE = 0.01


def read_fasta(path):
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith(">"):
                continue

            sequence.append(line.upper())

    result = "".join(sequence)

    if not result:
        raise RuntimeError(f"Empty reference: {path}")

    if any(base not in BASES for base in result):
        raise RuntimeError(
            "Reference contains characters other than A, C, G and T."
        )

    return result


def mutate_sequence(sequence, rate, rng):
    result = list(sequence)

    for index, base in enumerate(result):
        if rng.random() < rate:
            alternatives = [
                candidate
                for candidate in BASES
                if candidate != base
            ]
            result[index] = rng.choice(alternatives)

    return "".join(result)


def add_errors(sequence, rate, rng):
    result = list(sequence)

    for index, base in enumerate(result):
        if rng.random() < rate:
            alternatives = [
                candidate
                for candidate in BASES
                if candidate != base
            ]
            result[index] = rng.choice(alternatives)

    return "".join(result)


def reverse_complement(sequence):
    table = str.maketrans("ACGT", "TGCA")
    return sequence.translate(table)[::-1]


def write_fastq(
    haplotypes,
    sample,
    technology,
    read_length,
    step,
    error_rate,
    rng,
):
    output = (
        OUTPUT_DIR
        / f"{sample}.{technology}.chr21.smoke.fastq"
    )

    count = 0
    total_bases = 0

    with output.open("w") as handle:
        for hap_number, haplotype in enumerate(
            haplotypes,
            start=1,
        ):
            offset = 0 if hap_number == 1 else step // 2

            for start in range(
                offset,
                len(haplotype) - read_length + 1,
                step,
            ):
                read = haplotype[
                    start : start + read_length
                ]

                if rng.random() < 0.5:
                    read = reverse_complement(read)
                    orientation = "rev"
                else:
                    orientation = "fwd"

                read = add_errors(
                    read,
                    error_rate,
                    rng,
                )

                count += 1
                total_bases += len(read)

                name = (
                    f"{sample}_{technology}_"
                    f"hap{hap_number}_"
                    f"{count:06d}_"
                    f"{start + 1}_"
                    f"{orientation}"
                )

                handle.write(f"@{name}\n")
                handle.write(f"{read}\n")
                handle.write("+\n")
                handle.write("I" * len(read) + "\n")

    coverage = total_bases / len(haplotypes[0])

    print(
        f"{sample} {technology}: "
        f"{count} reads, "
        f"{coverage:.2f}x combined coverage"
    )


def main():
    if not REFERENCE.exists():
        raise FileNotFoundError(
            f"Reference not found: {REFERENCE}"
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    reference = read_fasta(REFERENCE)

    print(
        f"Reference length: "
        f"{len(reference):,} bp"
    )

    for sample_index, sample in enumerate(SAMPLES):
        rng = random.Random(
            SEED + sample_index
        )

        haplotype_1 = reference

        haplotype_2 = mutate_sequence(
            reference,
            HETEROZYGOUS_RATE,
            rng,
        )

        differences = sum(
            first != second
            for first, second in zip(
                haplotype_1,
                haplotype_2,
            )
        )

        print()
        print(
            f"{sample}: "
            f"{differences:,} haplotype differences "
            f"({differences / len(reference):.2%})"
        )

        write_fastq(
            [haplotype_1, haplotype_2],
            sample,
            "pb",
            HIFI_LENGTH,
            HIFI_STEP,
            HIFI_ERROR_RATE,
            rng,
        )

        write_fastq(
            [haplotype_1, haplotype_2],
            sample,
            "ont",
            ONT_LENGTH,
            ONT_STEP,
            ONT_ERROR_RATE,
            rng,
        )

    print()
    print("Diploid smoke reads created successfully.")


if __name__ == "__main__":
    main()
