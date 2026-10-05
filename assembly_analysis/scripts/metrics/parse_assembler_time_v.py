#!/usr/bin/env python3
"""Parse ``/usr/bin/time -v`` output file(s) into one row of a performance
table (wall-clock time, peak RAM, CPU time, CPU efficiency) for the
assembler RAM/time re-run.

Used by the five assembler *.ram_time.smk rule variants
(assemblers/whole_genome_asm/ram_time/{ont,pb}.assembly.{flye2,goldrush}.ram_time.smk,
hybrid.assembly.verkko.ram_time.smk) and the four aligner ones
(alignment_analysis/run_metrics/ram_time/{ont,pb}.read_mapping.{vg,vacmap}.ram_time.smk).
Pass ``--assembler`` or ``--aligner``; the first output column is named after
whichever was given. The aligners have three timed stages (map, sort, index)
that are summed like GoldRush stages. Snakemake's own ``benchmark:``
directive is deliberately not used for this: it samples the host-side
process tree, but every assembler here runs inside ``docker run`` --
Docker's real memory usage lives in a separate cgroup that a host-side
sampler cannot see, so ``benchmark:`` would report a near-zero, meaningless
number. Wrapping the actual tool binary with ``/usr/bin/time -v`` *inside*
the container avoids that; this is the same mechanism GoldRush's own
Makefile already uses internally (``track_time=1``), reused here uniformly
for Flye and Verkko too.

One time file per call for Flye/Verkko (the whole tool invocation is
wrapped once). Multiple time files per call for GoldRush, one per internal
Makefile stage (silver_path, golden_path, goldpolish, tigmint, 5x ntLink
rounds, ...) -- GoldRush has no single end-to-end number, so this script
sums each stage's elapsed wall-clock time and takes the max of each
stage's peak RSS, which is the standard way to report a multi-stage
pipeline's total wall-clock time and peak (not summed) memory footprint.
CPU time (user + system) is summed across stages like wall-clock time.

cpu_efficiency = CPU seconds / (wall-clock seconds * threads): 1.0 means all
requested threads were busy for the whole run; single-threaded stages pull
it down. CPU-hours, unlike wall-clock time, barely depend on the thread
count, so they stay comparable with the 32-thread production runs.

Every time file must report ``Exit status: 0``; a non-zero status means a
stage failed, and the script refuses to write a row for that run.

``--resource-samples`` takes the TSV written by
assemblers/whole_genome_asm/ram_time/sample_container_resources.sh during the
run (epoch_seconds, container_mem_bytes, scratch_disk_bytes; NA = not
measured). It adds:
  container_peak_mem_gb  -- highest memory of the whole container (all
                            processes together, page cache excluded). The
                            right peak for Verkko, which runs many jobs in
                            parallel; ``peak_rss_gb`` (largest single process)
                            undercounts it. Sampled every 30 s, so very short
                            spikes can be missed.
  peak_scratch_disk_gb   -- largest size of the scratch directory (sampled
                            every 5 min plus once at the end).

The elapsed-time line's format, confirmed from a real GoldRush time file
in this repository, is NOT the standard GNU coreutils colon format
(``h:mm:ss`` / ``m:ss``); it is unit-letter-suffixed and space-separated,
e.g. ``0m 30.22s`` or ``7m 57.79s``. The regex below accepts that format
and, defensively, the standard colon format too, in case a different
container's ``time`` build is used for Flye/Verkko.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def find_repo_root(start: Path) -> Path:
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "CONSTITUTION.md").is_file():
            return candidate
    raise RuntimeError("Could not locate lrs_benchmarking repository root")


PROJECT_ROOT = find_repo_root(Path(__file__))
DEFAULT_OUTPUT = (
    PROJECT_ROOT / "assembly_analysis" / "tables" / "assembler_performance_30x.tsv"
)

ELAPSED_LABEL = re.compile(
    r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)$", re.MULTILINE
)
MAX_RSS_LABEL = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)", re.MULTILINE)
USER_TIME_LABEL = re.compile(r"User time \(seconds\):\s*([\d.]+)", re.MULTILINE)
SYSTEM_TIME_LABEL = re.compile(r"System time \(seconds\):\s*([\d.]+)", re.MULTILINE)
EXIT_STATUS_LABEL = re.compile(r"Exit status:\s*(-?\d+)", re.MULTILINE)

# Unit-letter-suffixed format actually observed in this repo's time files,
# e.g. "0m 30.22s", "7m 57.79s", or (untested but consistent) "1h 5m 12.34s".
_LETTER_FORMAT = re.compile(
    r"^\s*(?:(?P<hours>\d+)h\s*)?(?:(?P<minutes>\d+)m\s*)?(?P<seconds>[\d.]+)s\s*$"
)
# Standard GNU coreutils colon format, e.g. "1:15:33" or "0:30.22", kept as
# a fallback in case a different `time` build is used for Flye/Verkko.
_COLON_FORMAT = re.compile(
    r"^\s*(?:(?P<hours>\d+):)?(?P<minutes>\d+):(?P<seconds>[\d.]+)\s*$"
)


def parse_elapsed_seconds(raw: str) -> float:
    raw = raw.strip()
    match = _LETTER_FORMAT.match(raw) or _COLON_FORMAT.match(raw)
    if not match:
        raise ValueError(f"Could not parse elapsed-time value: {raw!r}")
    parts = match.groupdict()
    hours = float(parts["hours"]) if parts["hours"] else 0.0
    minutes = float(parts["minutes"]) if parts["minutes"] else 0.0
    seconds = float(parts["seconds"])
    return hours * 3600.0 + minutes * 60.0 + seconds


def _require(pattern: re.Pattern, text: str, path: Path, label: str) -> str:
    match = pattern.search(text)
    if not match:
        raise ValueError(f"{path}: no '{label}' line found")
    return match.group(1)


def parse_time_file(path: Path) -> tuple[float, float, float, int]:
    """Return (elapsed_seconds, peak_rss_kbytes, cpu_seconds, exit_status)
    for one /usr/bin/time -v file."""
    text = path.read_text(encoding="utf-8")

    elapsed_seconds = parse_elapsed_seconds(
        _require(ELAPSED_LABEL, text, path, "Elapsed (wall clock) time")
    )
    peak_rss_kbytes = float(_require(MAX_RSS_LABEL, text, path, "Maximum resident set size"))
    cpu_seconds = float(_require(USER_TIME_LABEL, text, path, "User time")) + float(
        _require(SYSTEM_TIME_LABEL, text, path, "System time")
    )
    exit_status = int(_require(EXIT_STATUS_LABEL, text, path, "Exit status"))
    return elapsed_seconds, peak_rss_kbytes, cpu_seconds, exit_status


def parse_resource_samples(path: Path) -> tuple[float | None, float | None, int]:
    """Return (peak_mem_bytes, peak_disk_bytes, n_mem_samples); None = no value."""
    peak_mem = peak_disk = None
    n_mem = 0
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mem, disk = row["container_mem_bytes"], row["scratch_disk_bytes"]
            if mem not in ("", "NA"):
                n_mem += 1
                peak_mem = max(peak_mem or 0.0, float(mem))
            if disk not in ("", "NA"):
                peak_disk = max(peak_disk or 0.0, float(disk))
    return peak_mem, peak_disk, n_mem


def _gb(value: float | None) -> float | str:
    return "NA" if value is None else round(value / 1024**3, 4)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    tool = parser.add_mutually_exclusive_group(required=True)
    tool.add_argument("--assembler", choices=["flye", "goldrush", "verkko"])
    tool.add_argument("--aligner", choices=["vg_giraffe", "vacmap"])
    parser.add_argument("--sample", required=True, help="e.g. HG002")
    parser.add_argument("--technology", required=True, choices=["ont", "pb", "hybrid"])
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument(
        "--time-file", required=True, nargs="+", type=Path,
        help="One /usr/bin/time -v output file (Flye/Verkko), or several -- "
        "one per internal stage (GoldRush).",
    )
    parser.add_argument(
        "--resource-samples", type=Path,
        help="TSV from sample_container_resources.sh (container memory and "
        "scratch disk samples). Without it the two columns are NA.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    missing = [str(f) for f in args.time_file if not f.is_file()]
    if missing:
        raise FileNotFoundError(f"time file(s) not found: {missing}")

    elapsed_total = 0.0
    peak_rss_max = 0.0
    cpu_total = 0.0
    failed = []
    for time_file in args.time_file:
        elapsed_seconds, peak_rss_kbytes, cpu_seconds, exit_status = parse_time_file(time_file)
        elapsed_total += elapsed_seconds
        peak_rss_max = max(peak_rss_max, peak_rss_kbytes)
        cpu_total += cpu_seconds
        if exit_status != 0:
            failed.append(f"{time_file.name} (exit {exit_status})")
        print(
            f"  {time_file.name}: elapsed={elapsed_seconds:.2f}s "
            f"cpu={cpu_seconds:.2f}s "
            f"peak_rss={peak_rss_kbytes / 1_048_576:.3f} GB "
            f"exit={exit_status}"
        )

    if failed:
        raise SystemExit(
            "ERROR: non-zero exit status in time file(s), not recording a failed run: "
            + ", ".join(failed)
        )
    if elapsed_total <= 0:
        raise SystemExit("ERROR: total elapsed time is zero, cannot compute CPU efficiency")

    peak_mem = peak_disk = None
    n_mem = 0
    if args.resource_samples is not None:
        if not args.resource_samples.is_file():
            raise FileNotFoundError(f"resource samples not found: {args.resource_samples}")
        peak_mem, peak_disk, n_mem = parse_resource_samples(args.resource_samples)
        if peak_mem is None:
            print("WARNING: no container memory samples; container_peak_mem_gb = NA")
        if peak_disk is None:
            print("WARNING: no scratch disk samples; peak_scratch_disk_gb = NA")

    tool_column = "assembler" if args.assembler else "aligner"
    tool_name = args.assembler or args.aligner

    row = {
        tool_column: tool_name,
        "sample": args.sample,
        "technology": args.technology,
        "threads": args.threads,
        "wall_clock_seconds": round(elapsed_total, 2),
        "wall_clock_hours": round(elapsed_total / 3600.0, 4),
        "peak_rss_gb": round(peak_rss_max / 1_048_576, 4),
        "cpu_hours": round(cpu_total / 3600.0, 4),
        "cpu_efficiency": round(cpu_total / (elapsed_total * args.threads), 4),
        "exit_status": 0,
        "container_peak_mem_gb": _gb(peak_mem),
        "peak_scratch_disk_gb": _gb(peak_disk),
        "n_mem_samples": n_mem,
        "n_stages_summed": len(args.time_file),
        "source_time_files": ";".join(f.name for f in args.time_file),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    file_exists = args.output.is_file()
    with args.output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"\n{tool_name} {args.sample} {args.technology}:")
    print(f"  wall_clock: {row['wall_clock_seconds']} s ({row['wall_clock_hours']} h)")
    print(f"  peak_rss:   {row['peak_rss_gb']} GB")
    print(f"  cpu:        {row['cpu_hours']} CPU-h (efficiency {row['cpu_efficiency']})")
    print(f"  container:  {row['container_peak_mem_gb']} GB peak memory ({n_mem} samples)")
    print(f"  scratch:    {row['peak_scratch_disk_gb']} GB peak disk")
    print(f"  appended to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
