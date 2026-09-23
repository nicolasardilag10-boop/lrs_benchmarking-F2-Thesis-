#!/usr/bin/env python3

from __future__ import annotations

import csv
import html
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path("benchmark_chr21_real")

OUTPUT = ROOT / "config" / "sources.resolved.tsv"

LOG = (
    ROOT
    / "results"
    / "logs"
    / "source_validation"
    / "resolve_sources.log"
)

CHECKPOINT = (
    ROOT
    / "results"
    / "checkpoints"
    / "source"
    / "all_sources_validated.ok"
)

EXPECTED_CHR21_LENGTH = 46_709_983
MINIMUM_MAPQ = 20
TARGET_COVERAGE = 30


HIFI_DIRECTORIES = {
    "HG002": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG002_NA24385_son/"
        "PacBio_HiFi-Revio_20231031/"
    ),
    "HG003": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG003_NA24149_father/"
        "PacBio_HiFi-Revio_20231031/"
    ),
    "HG004": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG004_NA24143_mother/"
        "PacBio_HiFi-Revio_20231031/"
    ),
}


ONT_BAMS = {
    "HG002": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG002_NA24385_son/"
        "Ultralong_OxfordNanopore/"
        "guppy-V3.2.4_2020-01-22/"
        "HG002_GRCh38_ONT-UL_GIAB_20200122.phased.bam"
    ),
    "HG003": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG003_NA24149_father/"
        "UCSC_Ultralong_OxfordNanopore_Promethion/"
        "HG003_GRCh38_ONT-UL_UCSC_20200508.bam"
    ),
    "HG004": (
        "https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG004_NA24143_mother/"
        "UCSC_Ultralong_OxfordNanopore_Promethion/"
        "HG004_GRCh38_ONT-UL_UCSC_20200508.bam"
    ),
}


SEEDS = {
    ("HG002", "hifi"): 21002,
    ("HG003", "hifi"): 21003,
    ("HG004", "hifi"): 21004,
    ("HG002", "ont"): 22002,
    ("HG003", "ont"): 22003,
    ("HG004", "ont"): 22004,
}


def report(message: str) -> None:
    print(message, flush=True)

    LOG.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with LOG.open(
        "a",
        encoding="utf-8",
    ) as handle:
        handle.write(message + "\n")


def fetch_text(url: str) -> str:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "chr21-source-resolver/1.0",
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=120,
    ) as response:
        return response.read().decode(
            "utf-8",
            errors="replace",
        )


def discover_hifi_bam(
    sample: str,
    directory_url: str,
) -> str:
    report(
        f"Discovering GRCh38-GIABv3 Revio BAM for {sample}."
    )

    page = fetch_text(directory_url)

    hrefs = re.findall(
        r'href=["\']([^"\']+)["\']',
        page,
        flags=re.IGNORECASE,
    )

    candidates: list[str] = []

    for raw_href in hrefs:
        href = html.unescape(raw_href)

        filename = urllib.parse.unquote(
            Path(
                urllib.parse.urlparse(href).path
            ).name
        )

        if not filename.endswith(".bam"):
            continue

        if "PacBio-HiFi-Revio" not in filename:
            continue

        if "GRCh38-GIABv3" not in filename:
            continue

        candidates.append(
            urllib.parse.urljoin(
                directory_url,
                href,
            )
        )

    candidates = sorted(set(candidates))

    if len(candidates) != 1:
        raise RuntimeError(
            f"{sample}: expected exactly one matching HiFi BAM, "
            f"but found {len(candidates)}:\n"
            + "\n".join(candidates)
        )

    report(
        f"Resolved {sample} HiFi BAM: "
        f"{Path(candidates[0]).name}"
    )

    return candidates[0]


def run_samtools(
    arguments: list[str],
    timeout: int = 900,
) -> str:
    completed = subprocess.run(
        arguments,
        check=True,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    return completed.stdout


def parse_chr21_header(
    header: str,
) -> tuple[str, int]:
    targets: dict[str, int] = {}

    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue

        name: str | None = None
        reference_length: int | None = None

        for field in line.split("\t"):
            if field.startswith("SN:"):
                name = field[3:]

            elif field.startswith("LN:"):
                reference_length = int(field[3:])

        if (
            name is not None
            and reference_length is not None
        ):
            targets[name] = reference_length

    if "chr21" in targets:
        return "chr21", targets["chr21"]

    if "21" in targets:
        return "21", targets["21"]

    raise RuntimeError(
        "Chromosome 21 was not found in the remote BAM header."
    )


def validate_remote_bam(
    sample: str,
    technology: str,
    bam_url: str,
) -> tuple[str, int, int]:
    report(
        f"Reading remote header for {sample} {technology}."
    )

    header = run_samtools(
        [
            "samtools",
            "view",
            "-H",
            bam_url,
        ]
    )

    chromosome_name, chromosome_length = (
        parse_chr21_header(header)
    )

    if chromosome_length != EXPECTED_CHR21_LENGTH:
        raise RuntimeError(
            f"{sample} {technology}: chromosome-length mismatch. "
            f"Observed={chromosome_length}; "
            f"expected={EXPECTED_CHR21_LENGTH}."
        )

    report(
        f"Reading remote index for {sample} {technology}."
    )

    idxstats = run_samtools(
        [
            "samtools",
            "idxstats",
            bam_url,
        ]
    )

    mapped_records: int | None = None

    for line in idxstats.splitlines():
        fields = line.split("\t")

        if len(fields) < 4:
            continue

        if fields[0] == chromosome_name:
            mapped_records = int(fields[2])
            break

    if mapped_records is None:
        raise RuntimeError(
            f"{sample} {technology}: chromosome "
            f"{chromosome_name} was absent from idxstats."
        )

    if mapped_records <= 0:
        raise RuntimeError(
            f"{sample} {technology}: the remote index reports "
            "zero chromosome-21 records."
        )

    report(
        f"PASS: {sample} {technology}; "
        f"{chromosome_name}; "
        f"{chromosome_length:,} bp; "
        f"{mapped_records:,} indexed records."
    )

    return (
        chromosome_name,
        chromosome_length,
        mapped_records,
    )


def main() -> None:
    if OUTPUT.exists():
        raise RuntimeError(
            f"Refusing to overwrite existing manifest: {OUTPUT}"
        )

    LOG.unlink(missing_ok=True)

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    CHECKPOINT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report("=" * 72)
    report("CHROMOSOME-21 SOURCE RESOLUTION")
    report("=" * 72)

    started = time.monotonic()

    hifi_bams = {
        sample: discover_hifi_bam(
            sample,
            directory_url,
        )
        for sample, directory_url
        in HIFI_DIRECTORIES.items()
    }

    candidates: list[tuple[str, str, str, str]] = []

    for sample in ["HG002", "HG003", "HG004"]:
        candidates.append(
            (
                sample,
                "hifi",
                hifi_bams[sample],
                "PacBio_HiFi_Revio_20231031",
            )
        )

        candidates.append(
            (
                sample,
                "ont",
                ONT_BAMS[sample],
                "ONT_ultralong_GRCh38",
            )
        )

    rows: list[dict[str, object]] = []
    total = len(candidates)

    for index, (
        sample,
        technology,
        bam_url,
        source_family,
    ) in enumerate(
        candidates,
        start=1,
    ):
        elapsed = time.monotonic() - started
        completed = index - 1

        if completed > 0:
            estimated_total = (
                elapsed
                * total
                / completed
            )

            remaining = max(
                estimated_total - elapsed,
                0,
            )
        else:
            remaining = 0

        percentage = round(
            100
            * completed
            / total
        )

        report("")
        report(
            f"[{percentage:3d}%] Dataset {index}/{total}: "
            f"{sample} {technology}; "
            f"elapsed={elapsed / 60:.1f} min; "
            f"estimated remaining={remaining / 60:.1f} min."
        )

        (
            chromosome_name,
            chromosome_length,
            indexed_records,
        ) = validate_remote_bam(
            sample,
            technology,
            bam_url,
        )

        rows.append(
            {
                "sample": sample,
                "technology": technology,
                "source_family": source_family,
                "source_bam_url": bam_url,
                "source_bai_url": bam_url + ".bai",
                "chromosome_name": chromosome_name,
                "chromosome_length": chromosome_length,
                "minimum_mapq": MINIMUM_MAPQ,
                "target_coverage": TARGET_COVERAGE,
                "seed": SEEDS[(sample, technology)],
                "remote_chr21_indexed_records": indexed_records,
                "source_validation": "PASS",
            }
        )

    rows.sort(
        key=lambda row: (
            str(row["sample"]),
            str(row["technology"]),
        )
    )

    fields = [
        "sample",
        "technology",
        "source_family",
        "source_bam_url",
        "source_bai_url",
        "chromosome_name",
        "chromosome_length",
        "minimum_mapq",
        "target_coverage",
        "seed",
        "remote_chr21_indexed_records",
        "source_validation",
    ]

    temporary_output = Path(
        str(OUTPUT) + ".partial"
    )

    with temporary_output.open(
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
        writer.writerows(rows)

    temporary_output.replace(OUTPUT)

    CHECKPOINT.write_text(
        "\n".join(
            [
                "ALL REMOTE SOURCES VALIDATED",
                f"Datasets: {len(rows)}",
                (
                    "Reference chromosome length: "
                    f"{EXPECTED_CHR21_LENGTH}"
                ),
                f"Manifest: {OUTPUT.resolve()}",
                "STATUS: PASS",
                "",
            ]
        ),
        encoding="utf-8",
    )

    elapsed = time.monotonic() - started

    report("")
    report("=" * 72)
    report("ALL SIX SOURCES VALIDATED")
    report("=" * 72)
    report(f"Elapsed: {elapsed / 60:.2f} minutes")
    report(f"Manifest: {OUTPUT}")
    report(f"Checkpoint: {CHECKPOINT}")


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        LOG.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with LOG.open(
            "a",
            encoding="utf-8",
        ) as handle:
            handle.write(f"\nERROR: {error}\n")

        print(
            f"\nERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)
