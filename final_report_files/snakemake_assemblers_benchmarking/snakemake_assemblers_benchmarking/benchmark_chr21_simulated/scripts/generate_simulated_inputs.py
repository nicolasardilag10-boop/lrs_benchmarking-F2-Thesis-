#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REFERENCE = (
    ROOT
    / "inputs"
    / "reference"
    / "chr21_region_1mb.fasta"
)

TRUTH_DIR = ROOT / "inputs" / "truth"
READ_DIR = ROOT / "inputs" / "reads"
CONFIG_DIR = ROOT / "configs"
VALIDATION_DIR = ROOT / "results" / "validation"
LOG_DIR = ROOT / "results" / "logs"

SAMPLES = ["SIM01", "SIM02", "SIM03"]

BASE_SEED = 20260805
HETEROZYGOUS_RATE = 0.001
TARGET_COVERAGE_PER_HAPLOTYPE = 30.0

TECHNOLOGIES = {
    "pb": {
        "platform": "PacBio HiFi",
        "mean_length": 15_000,
        "length_cv": 0.20,
        "minimum_length": 5_000,
        "maximum_length": 25_000,
        "error_rate": 0.002,
        "substitution_fraction": 0.80,
        "deletion_fraction": 0.10,
        "insertion_fraction": 0.10,
        "seed_offset": 10,
    },
    "ont": {
        "platform": "Oxford Nanopore",
        "mean_length": 20_000,
        "length_cv": 0.60,
        "minimum_length": 3_000,
        "maximum_length": 80_000,
        "error_rate": 0.030,
        "substitution_fraction": 0.55,
        "deletion_fraction": 0.25,
        "insertion_fraction": 0.20,
        "seed_offset": 20,
    },
}

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_fasta(path: Path) -> str:
    sequence_parts: list[str] = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="strict",
    ) as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith(">"):
                continue

            sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)

    if not sequence:
        raise ValueError(f"No sequence found in {path}")

    unexpected = sorted(set(sequence) - set(BASES))

    if unexpected:
        raise ValueError(
            f"Reference contains unsupported bases: {unexpected}"
        )

    return sequence


def write_fasta(
    path: Path,
    identifier: str,
    sequence: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{identifier}\n")

        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def alternate_base(
    base: str,
    rng: random.Random,
) -> str:
    choices = [candidate for candidate in BASES if candidate != base]
    return rng.choice(choices)


def create_haplotypes(
    reference: str,
    rng: random.Random,
) -> tuple[str, str, list[dict[str, object]]]:
    """
    Create exactly the requested pairwise heterozygosity.

    At every selected heterozygous position, either haplotype 1
    or haplotype 2 receives one alternate allele.
    """
    heterozygous_sites = max(
        1,
        round(len(reference) * HETEROZYGOUS_RATE),
    )

    selected_positions = sorted(
        rng.sample(
            range(len(reference)),
            heterozygous_sites,
        )
    )

    haplotype_1 = list(reference)
    haplotype_2 = list(reference)
    variants: list[dict[str, object]] = []

    for position in selected_positions:
        reference_base = reference[position]
        alternate = alternate_base(reference_base, rng)

        if rng.random() < 0.5:
            haplotype_1[position] = alternate
            genotype = "1|0"
            haplotype_changed = "hap1"
        else:
            haplotype_2[position] = alternate
            genotype = "0|1"
            haplotype_changed = "hap2"

        variants.append(
            {
                "position_1_based": position + 1,
                "reference": reference_base,
                "alternate": alternate,
                "genotype": genotype,
                "haplotype_changed": haplotype_changed,
            }
        )

    hap1 = "".join(haplotype_1)
    hap2 = "".join(haplotype_2)

    observed_differences = sum(
        first != second
        for first, second in zip(hap1, hap2)
    )

    if observed_differences != heterozygous_sites:
        raise RuntimeError(
            "Generated haplotypes do not contain the expected "
            "number of differences."
        )

    return hap1, hap2, variants


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def sample_read_length(
    parameters: dict[str, object],
    rng: random.Random,
) -> int:
    mean_length = float(parameters["mean_length"])
    coefficient_of_variation = float(parameters["length_cv"])

    sigma = math.sqrt(
        math.log(1 + coefficient_of_variation**2)
    )

    mu = math.log(mean_length) - (sigma**2 / 2)

    sampled = round(
        rng.lognormvariate(mu, sigma)
    )

    minimum = int(parameters["minimum_length"])
    maximum = int(parameters["maximum_length"])

    return max(minimum, min(sampled, maximum))


def apply_sequencing_errors(
    sequence: str,
    parameters: dict[str, object],
    rng: random.Random,
) -> str:
    total_error = float(parameters["error_rate"])

    substitution_rate = (
        total_error
        * float(parameters["substitution_fraction"])
    )

    deletion_rate = (
        total_error
        * float(parameters["deletion_fraction"])
    )

    insertion_rate = (
        total_error
        * float(parameters["insertion_fraction"])
    )

    output: list[str] = []

    for base in sequence:
        event = rng.random()

        if event < deletion_rate:
            pass
        elif event < deletion_rate + substitution_rate:
            output.append(alternate_base(base, rng))
        else:
            output.append(base)

        if rng.random() < insertion_rate:
            output.append(rng.choice(BASES))

    if not output:
        raise RuntimeError("Sequencing simulation produced an empty read.")

    return "".join(output)


def quality_character(error_rate: float) -> str:
    phred = round(-10 * math.log10(error_rate))
    phred = max(2, min(phred, 40))

    return chr(phred + 33)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def generate_fastq(
    sample: str,
    technology: str,
    haplotypes: tuple[str, str],
    seed: int,
) -> dict[str, object]:
    parameters = TECHNOLOGIES[technology]
    rng = random.Random(seed)

    output_path = READ_DIR / f"{sample}.{technology}.fastq.gz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    target_bases = [
        round(
            len(haplotype)
            * TARGET_COVERAGE_PER_HAPLOTYPE
        )
        for haplotype in haplotypes
    ]

    source_bases = [0, 0]
    emitted_bases = [0, 0]
    read_counts = [0, 0]
    emitted_lengths: list[int] = []
    sequence_hashes: set[str] = set()

    qchar = quality_character(
        float(parameters["error_rate"])
    )

    read_number = 0

    with gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        compresslevel=6,
    ) as handle:
        while any(
            source_bases[index] < target_bases[index]
            for index in (0, 1)
        ):
            available_haplotypes = [
                index
                for index in (0, 1)
                if source_bases[index] < target_bases[index]
            ]

            haplotype_index = rng.choice(
                available_haplotypes
            )

            haplotype = haplotypes[haplotype_index]

            source_length = sample_read_length(
                parameters,
                rng,
            )

            source_length = min(
                source_length,
                len(haplotype),
            )

            start = rng.randrange(
                0,
                len(haplotype) - source_length + 1,
            )

            source_sequence = haplotype[
                start : start + source_length
            ]

            strand = "+"

            if rng.random() < 0.5:
                source_sequence = reverse_complement(
                    source_sequence
                )
                strand = "-"

            emitted_sequence = apply_sequencing_errors(
                source_sequence,
                parameters,
                rng,
            )

            read_number += 1
            read_counts[haplotype_index] += 1
            source_bases[haplotype_index] += source_length
            emitted_bases[haplotype_index] += len(
                emitted_sequence
            )

            emitted_lengths.append(
                len(emitted_sequence)
            )

            sequence_hashes.add(
                hashlib.sha256(
                    emitted_sequence.encode()
                ).hexdigest()
            )

            identifier = (
                f"{sample}_{technology}_"
                f"read{read_number:06d}_"
                f"hap{haplotype_index + 1}_"
                f"start{start + 1}_"
                f"source{source_length}_"
                f"strand{strand}"
            )

            quality = qchar * len(emitted_sequence)

            handle.write(f"@{identifier}\n")
            handle.write(emitted_sequence + "\n")
            handle.write("+\n")
            handle.write(quality + "\n")

    total_reads = sum(read_counts)
    total_source_bases = sum(source_bases)
    total_emitted_bases = sum(emitted_bases)

    return {
        "sample": sample,
        "technology": technology,
        "platform": parameters["platform"],
        "seed": seed,
        "read_path": str(output_path.relative_to(ROOT)),
        "target_coverage_per_haplotype": (
            TARGET_COVERAGE_PER_HAPLOTYPE
        ),
        "hap1_read_count": read_counts[0],
        "hap2_read_count": read_counts[1],
        "total_read_count": total_reads,
        "hap1_source_bases": source_bases[0],
        "hap2_source_bases": source_bases[1],
        "hap1_nominal_coverage": (
            source_bases[0] / len(haplotypes[0])
        ),
        "hap2_nominal_coverage": (
            source_bases[1] / len(haplotypes[1])
        ),
        "total_source_bases": total_source_bases,
        "total_emitted_bases": total_emitted_bases,
        "mean_emitted_read_length": (
            total_emitted_bases / total_reads
        ),
        "minimum_emitted_read_length": min(
            emitted_lengths
        ),
        "maximum_emitted_read_length": max(
            emitted_lengths
        ),
        "unique_emitted_sequences": len(
            sequence_hashes
        ),
        "fastq_sha256": sha256_file(output_path),
    }


def remove_previous_outputs() -> None:
    patterns = [
        TRUTH_DIR / "SIM*.fasta",
        TRUTH_DIR / "SIM*.tsv",
        READ_DIR / "SIM*.fastq.gz",
    ]

    for pattern in patterns:
        for path in pattern.parent.glob(pattern.name):
            path.unlink()

    for path in [
        CONFIG_DIR / "samples.tsv",
        VALIDATION_DIR / "simulation_manifest.tsv",
        VALIDATION_DIR / "simulation_provenance.json",
        VALIDATION_DIR / "checksums.sha256",
    ]:
        if path.exists():
            path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate reproducible synthetic diploid long-read "
            "benchmark inputs."
        )
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace previously generated simulation outputs.",
    )

    arguments = parser.parse_args()

    TRUTH_DIR.mkdir(parents=True, exist_ok=True)
    READ_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    existing_outputs = (
        list(TRUTH_DIR.glob("SIM*.fasta"))
        + list(READ_DIR.glob("SIM*.fastq.gz"))
    )

    if existing_outputs and not arguments.force:
        raise SystemExit(
            "Generated outputs already exist. Use --force "
            "to regenerate them."
        )

    if arguments.force:
        remove_previous_outputs()

    reference = read_fasta(REFERENCE)

    if len(reference) != 1_000_000:
        raise ValueError(
            f"Expected a 1,000,000-bp reference, found "
            f"{len(reference):,} bp."
        )

    manifest_rows: list[dict[str, object]] = []
    truth_summary: list[dict[str, object]] = []

    for sample_index, sample in enumerate(SAMPLES):
        truth_seed = BASE_SEED + sample_index * 100

        truth_rng = random.Random(truth_seed)

        hap1, hap2, variants = create_haplotypes(
            reference,
            truth_rng,
        )

        hap1_path = TRUTH_DIR / f"{sample}.hap1.fasta"
        hap2_path = TRUTH_DIR / f"{sample}.hap2.fasta"
        variant_path = (
            TRUTH_DIR
            / f"{sample}.heterozygous_snps.tsv"
        )

        write_fasta(
            hap1_path,
            f"{sample}_haplotype_1",
            hap1,
        )

        write_fasta(
            hap2_path,
            f"{sample}_haplotype_2",
            hap2,
        )

        with variant_path.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            fieldnames = [
                "position_1_based",
                "reference",
                "alternate",
                "genotype",
                "haplotype_changed",
            ]

            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=fieldnames,
            )

            writer.writeheader()
            writer.writerows(variants)

        observed_difference = sum(
            first != second
            for first, second in zip(hap1, hap2)
        )

        truth_summary.append(
            {
                "sample": sample,
                "truth_seed": truth_seed,
                "haplotype_length": len(hap1),
                "heterozygous_sites": observed_difference,
                "observed_heterozygosity": (
                    observed_difference / len(hap1)
                ),
                "hap1_sha256": sha256_file(hap1_path),
                "hap2_sha256": sha256_file(hap2_path),
                "variant_table_sha256": sha256_file(
                    variant_path
                ),
            }
        )

        for technology in ("pb", "ont"):
            simulation_seed = (
                truth_seed
                + int(
                    TECHNOLOGIES[technology][
                        "seed_offset"
                    ]
                )
            )

            manifest_rows.append(
                generate_fastq(
                    sample=sample,
                    technology=technology,
                    haplotypes=(hap1, hap2),
                    seed=simulation_seed,
                )
            )

    samples_path = CONFIG_DIR / "samples.tsv"

    with samples_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            ["sample", "technology", "fastq"]
        )

        for sample in SAMPLES:
            for technology in ("ont", "pb"):
                writer.writerow(
                    [
                        sample,
                        technology,
                        (
                            "benchmark_chr21_simulated/"
                            "inputs/reads/"
                            f"{sample}.{technology}.fastq.gz"
                        ),
                    ]
                )

    manifest_path = (
        VALIDATION_DIR
        / "simulation_manifest.tsv"
    )

    manifest_fields = list(manifest_rows[0].keys())

    with manifest_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=manifest_fields,
        )

        writer.writeheader()
        writer.writerows(manifest_rows)

    provenance_path = (
        VALIDATION_DIR
        / "simulation_provenance.json"
    )

    provenance = {
        "experiment": (
            "chr21_1mb_synthetic_diploid_benchmark"
        ),
        "reference": str(
            REFERENCE.relative_to(ROOT)
        ),
        "reference_length": len(reference),
        "base_seed": BASE_SEED,
        "heterozygous_rate": HETEROZYGOUS_RATE,
        "target_coverage_per_haplotype": (
            TARGET_COVERAGE_PER_HAPLOTYPE
        ),
        "samples": SAMPLES,
        "technologies": TECHNOLOGIES,
        "truth_summary": truth_summary,
        "real_biological_samples": False,
        "allowed_description": (
            "Controlled synthetic diploid benchmark "
            "using a 1-Mb region derived from "
            "human chromosome 21."
        ),
    }

    with provenance_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            provenance,
            handle,
            indent=2,
            sort_keys=True,
        )

        handle.write("\n")

    checksum_paths = [
        REFERENCE,
        samples_path,
        manifest_path,
        provenance_path,
        *sorted(TRUTH_DIR.glob("SIM*")),
        *sorted(READ_DIR.glob("SIM*.fastq.gz")),
    ]

    checksum_path = (
        VALIDATION_DIR
        / "checksums.sha256"
    )

    with checksum_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for path in checksum_paths:
            handle.write(
                f"{sha256_file(path)}  "
                f"{path.relative_to(ROOT)}\n"
            )

    print("=== SIMULATION COMPLETED ===")
    print(f"Reference length: {len(reference):,} bp")
    print(
        "Target coverage: "
        f"{TARGET_COVERAGE_PER_HAPLOTYPE:.0f}x "
        "per haplotype"
    )
    print()

    for truth in truth_summary:
        print(
            f"{truth['sample']}: "
            f"{truth['heterozygous_sites']:,} "
            "heterozygous SNPs "
            f"({float(truth['observed_heterozygosity']):.3%})"
        )

    print()

    for row in manifest_rows:
        print(
            f"{row['sample']} "
            f"{str(row['technology']).upper()}: "
            f"{row['total_read_count']} reads, "
            f"{int(row['total_emitted_bases']):,} emitted bases, "
            f"{float(row['hap1_nominal_coverage']):.2f}x hap1, "
            f"{float(row['hap2_nominal_coverage']):.2f}x hap2, "
            f"{row['unique_emitted_sequences']} unique reads"
        )

    print()
    print("Created:")
    print(f"  {samples_path.relative_to(ROOT)}")
    print(f"  {manifest_path.relative_to(ROOT)}")
    print(f"  {provenance_path.relative_to(ROOT)}")
    print(f"  {checksum_path.relative_to(ROOT)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
