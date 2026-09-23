from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
ROOT = PROJECT / "benchmark_chr21_real"

BAM_URL = (
    "https://downloads.pacbcloud.com/public/dataset/"
    "HG002-CpG-methylation-202202/"
    "HG002.GRCh38.haplotagged.bam"
)

OUTPUT_DIR = ROOT / "inputs" / "pacbio" / "HG002"
CHUNK_DIR = OUTPUT_DIR / "chr21_chunks"

LOG_DIR = (
    ROOT
    / "results"
    / "logs"
    / "extraction"
    / "HG002_hifi_chr21_chunks"
)

VALIDATION_DIR = ROOT / "results" / "validation"

FILTER_SCRIPT = (
    ROOT
    / "scripts"
    / "filter_sam_window.py"
)

FINAL_BAM = (
    OUTPUT_DIR
    / "HG002.hifi.chr21.primary.q20.bam"
)

SUMMARY = (
    VALIDATION_DIR
    / "HG002.hifi.chr21.extraction.tsv"
)

CHROMOSOME = "chr21"
CHROMOSOME_LENGTH = 46_709_983
WINDOW_SIZE = 1_000_000

# Unmapped + secondary + supplementary.
EXCLUDED_FLAGS = 2308
MINIMUM_MAPQ = 20

MAX_ATTEMPTS = 5


def format_seconds(value: float | None) -> str:
    if value is None or value < 0:
        return "--:--:--"

    value = int(value)

    hours, remainder = divmod(value, 3600)
    minutes, seconds = divmod(remainder, 60)

    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def quickcheck(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False

    completed = subprocess.run(
        ["samtools", "quickcheck", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    return completed.returncode == 0


def read_count(path: Path) -> int:
    result = subprocess.run(
        ["samtools", "view", "-c", str(path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )

    return int(result.stdout.strip())


def chunk_path(
    window_number: int,
    window_start: int,
    window_end: int,
) -> Path:
    return CHUNK_DIR / (
        f"chunk_{window_number:03d}_"
        f"{window_start:09d}_"
        f"{window_end:09d}.bam"
    )


def download_window(
    *,
    window_number: int,
    window_start: int,
    window_end: int,
    attempt: int,
    partial: Path,
) -> tuple[bool, Path]:
    region = (
        f"{CHROMOSOME}:"
        f"{window_start}-{window_end}"
    )

    log_path = LOG_DIR / (
        f"window_{window_number:03d}_"
        f"attempt_{attempt}.log"
    )

    partial.unlink(missing_ok=True)

    with log_path.open("wb") as log_handle:
        remote_process = subprocess.Popen(
            [
                "samtools",
                "view",
                "-h",
                BAM_URL,
                region,
            ],
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )

        assert remote_process.stdout is not None

        filter_process = subprocess.Popen(
            [
                sys.executable,
                str(FILTER_SCRIPT),
                CHROMOSOME,
                str(window_start),
                str(window_end),
                str(EXCLUDED_FLAGS),
                str(MINIMUM_MAPQ),
            ],
            stdin=remote_process.stdout,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )

        remote_process.stdout.close()

        assert filter_process.stdout is not None

        bam_process = subprocess.Popen(
            [
                "samtools",
                "view",
                "-@",
                "2",
                "-b",
                "-o",
                str(partial),
                "-",
            ],
            stdin=filter_process.stdout,
            stderr=log_handle,
        )

        filter_process.stdout.close()

        bam_status = bam_process.wait()
        filter_status = filter_process.wait()
        remote_status = remote_process.wait()

    successful = (
        remote_status == 0
        and filter_status == 0
        and bam_status == 0
        and quickcheck(partial)
    )

    return successful, log_path


def calculate_total_sequence_bases(
    bam_path: Path,
) -> int:
    process = subprocess.Popen(
        ["samtools", "view", str(bam_path)],
        text=True,
        stdout=subprocess.PIPE,
    )

    assert process.stdout is not None

    total_bases = 0

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            continue

        sequence = fields[9]

        if sequence != "*":
            total_bases += len(sequence)

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            "samtools view failed while calculating bases"
        )

    return total_bases


def mapped_cigar_bases(bam_path: Path) -> int:
    result = subprocess.run(
        ["samtools", "stats", str(bam_path)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )

    for line in result.stdout.splitlines():
        if line.startswith("SN\tbases mapped (cigar):"):
            return int(line.split("\t")[2])

    return 0


def main() -> None:
    for program in ["samtools", sys.executable]:
        if shutil.which(program) is None:
            raise SystemExit(
                f"ERROR: required program not found: {program}"
            )

    if not FILTER_SCRIPT.is_file():
        raise SystemExit(
            f"ERROR: missing filter: {FILTER_SCRIPT}"
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

    total_windows = (
        CHROMOSOME_LENGTH
        + WINDOW_SIZE
        - 1
    ) // WINDOW_SIZE

    windows: list[tuple[int, int, int, Path]] = []

    for window_number in range(1, total_windows + 1):
        window_start = (
            (window_number - 1) * WINDOW_SIZE
            + 1
        )

        window_end = min(
            window_number * WINDOW_SIZE,
            CHROMOSOME_LENGTH,
        )

        windows.append(
            (
                window_number,
                window_start,
                window_end,
                chunk_path(
                    window_number,
                    window_start,
                    window_end,
                ),
            )
        )

    valid_existing = {
        window_number
        for (
            window_number,
            window_start,
            window_end,
            chunk,
        ) in windows
        if quickcheck(chunk)
    }

    remaining_at_start = (
        total_windows - len(valid_existing)
    )

    print(
        "[  0%] Starting resumable HG002 HiFi "
        "chromosome-21 extraction.",
        flush=True,
    )

    print(
        f"       Windows: {total_windows}",
        flush=True,
    )

    print(
        f"       Valid existing windows: "
        f"{len(valid_existing)}",
        flush=True,
    )

    print(
        f"       Windows remaining: "
        f"{remaining_at_start}",
        flush=True,
    )

    print(
        "       ETA becomes more reliable after "
        "3–5 new windows.",
        flush=True,
    )

    started = time.monotonic()
    completed_new = 0
    total_selected = 0

    for (
        window_number,
        window_start,
        window_end,
        chunk,
    ) in windows:
        region = (
            f"{CHROMOSOME}:"
            f"{window_start}-{window_end}"
        )

        if quickcheck(chunk):
            selected = read_count(chunk)
            total_selected += selected

            done = sum(
                quickcheck(existing_chunk)
                for (
                    existing_number,
                    existing_start,
                    existing_end,
                    existing_chunk,
                ) in windows
            )

            percentage = 10 + int(
                done / total_windows * 70
            )

            print(
                f"[{percentage:3d}%] "
                f"Reused window "
                f"{window_number}/{total_windows} | "
                f"{region} | "
                f"selected {selected:,}",
                flush=True,
            )

            continue

        success = False

        for attempt in range(1, MAX_ATTEMPTS + 1):
            elapsed = time.monotonic() - started

            if completed_new > 0:
                rate = completed_new / elapsed

                remaining_new = (
                    remaining_at_start
                    - completed_new
                )

                eta = (
                    remaining_new / rate
                    if rate > 0
                    else None
                )
            else:
                eta = None

            done_before = (
                len(valid_existing)
                + completed_new
            )

            percentage = 10 + int(
                done_before / total_windows * 70
            )

            print(
                f"[{percentage:3d}%] "
                f"Window {window_number}/{total_windows} | "
                f"{region} | "
                f"attempt {attempt}/{MAX_ATTEMPTS} | "
                f"elapsed {format_seconds(elapsed)} | "
                f"ETA {format_seconds(eta)}",
                flush=True,
            )

            partial = Path(f"{chunk}.partial")

            successful, log_path = download_window(
                window_number=window_number,
                window_start=window_start,
                window_end=window_end,
                attempt=attempt,
                partial=partial,
            )

            if successful:
                os.replace(partial, chunk)
                success = True
                break

            partial.unlink(missing_ok=True)

            delay = attempt * 10

            print(
                f"       Window failed; retrying in "
                f"{delay} seconds.",
                flush=True,
            )

            print(
                f"       Log: {log_path}",
                flush=True,
            )

            time.sleep(delay)

        if not success:
            raise SystemExit(
                "\nERROR: window "
                f"{window_number} failed after "
                f"{MAX_ATTEMPTS} attempts.\n"
                "Rerun the same command to resume."
            )

        selected = read_count(chunk)

        total_selected += selected
        completed_new += 1

        elapsed = time.monotonic() - started
        rate = completed_new / elapsed

        remaining_new = (
            remaining_at_start - completed_new
        )

        eta = (
            remaining_new / rate
            if rate > 0
            else None
        )

        completed_total = (
            len(valid_existing)
            + completed_new
        )

        percentage = 10 + int(
            completed_total / total_windows * 70
        )

        print(
            f"[{percentage:3d}%] "
            f"Completed window "
            f"{window_number}/{total_windows} | "
            f"selected {selected:,} | "
            f"total selected {total_selected:,} | "
            f"elapsed {format_seconds(elapsed)} | "
            f"ETA {format_seconds(eta)}",
            flush=True,
        )

    print(
        "[ 82%] All chromosome windows completed.",
        flush=True,
    )

    chunks = [
        chunk
        for (
            window_number,
            window_start,
            window_end,
            chunk,
        ) in windows
    ]

    invalid_chunks = [
        str(chunk)
        for chunk in chunks
        if not quickcheck(chunk)
    ]

    if invalid_chunks:
        raise SystemExit(
            "ERROR: invalid chunks remain:\n"
            + "\n".join(invalid_chunks)
        )

    print(
        "[ 85%] Concatenating BAM chunks.",
        flush=True,
    )

    partial_final = Path(f"{FINAL_BAM}.partial")

    partial_final.unlink(missing_ok=True)
    FINAL_BAM.unlink(missing_ok=True)
    Path(f"{FINAL_BAM}.bai").unlink(
        missing_ok=True
    )

    subprocess.run(
        [
            "samtools",
            "cat",
            "-o",
            str(partial_final),
            *[str(chunk) for chunk in chunks],
        ],
        check=True,
    )

    if not quickcheck(partial_final):
        raise SystemExit(
            "ERROR: concatenated BAM failed validation"
        )

    os.replace(partial_final, FINAL_BAM)

    print(
        "[ 90%] Indexing final chromosome-21 BAM.",
        flush=True,
    )

    subprocess.run(
        [
            "samtools",
            "index",
            "-@",
            "2",
            str(FINAL_BAM),
        ],
        check=True,
    )

    print(
        "[ 94%] Calculating read and coverage metrics.",
        flush=True,
    )

    reads = read_count(FINAL_BAM)

    total_read_bases = (
        calculate_total_sequence_bases(FINAL_BAM)
    )

    cigar_bases = mapped_cigar_bases(FINAL_BAM)

    input_coverage = (
        total_read_bases / CHROMOSOME_LENGTH
    )

    cigar_coverage = (
        cigar_bases / CHROMOSOME_LENGTH
    )

    status = (
        "READY_FOR_30X_DOWNSAMPLING"
        if input_coverage >= 30
        else "BELOW_30X"
    )

    with SUMMARY.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "sample",
                "technology",
                "chromosome",
                "filter",
                "reads",
                "total_read_bases",
                "chromosome_length",
                "estimated_input_coverage",
                "mapped_cigar_bases",
                "mapped_cigar_coverage",
                "target_coverage",
                "coverage_status",
                "bam",
            ]
        )

        writer.writerow(
            [
                "HG002",
                "PacBio_HiFi",
                "chr21",
                "primary_MAPQ20",
                reads,
                total_read_bases,
                CHROMOSOME_LENGTH,
                f"{input_coverage:.4f}",
                cigar_bases,
                f"{cigar_coverage:.4f}",
                "30",
                status,
                FINAL_BAM,
            ]
        )

    print(
        "[ 98%] Running final integrity check.",
        flush=True,
    )

    if not quickcheck(FINAL_BAM):
        raise SystemExit(
            "ERROR: final BAM integrity check failed"
        )

    elapsed = time.monotonic() - started

    print(
        "[100%] HG002 HiFi chromosome-21 "
        "extraction completed.",
        flush=True,
    )

    print(
        f"Total elapsed time: "
        f"{format_seconds(elapsed)}",
        flush=True,
    )

    print(
        f"Reads: {reads:,}",
        flush=True,
    )

    print(
        f"Total read bases: "
        f"{total_read_bases:,}",
        flush=True,
    )

    print(
        f"Estimated input coverage: "
        f"{input_coverage:.4f}x",
        flush=True,
    )

    print(
        f"Mapped-CIGAR coverage: "
        f"{cigar_coverage:.4f}x",
        flush=True,
    )

    print(
        f"Status: {status}",
        flush=True,
    )

    print(
        f"Summary: {SUMMARY}",
        flush=True,
    )


if __name__ == "__main__":
    main()
