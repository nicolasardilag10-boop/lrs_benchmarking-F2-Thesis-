from __future__ import annotations

import sys
import time


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: filter_sam_progress.py "
        "TOTAL START_PERCENT END_PERCENT "
        "EXCLUDED_FLAGS MINIMUM_MAPQ"
    )


total = int(sys.argv[1])
start_percent = int(sys.argv[2])
end_percent = int(sys.argv[3])
excluded_flags = int(sys.argv[4])
minimum_mapq = int(sys.argv[5])

processed = 0
selected = 0

start_time = time.monotonic()
last_report_time = 0.0
last_displayed_percent = -1


def format_seconds(seconds: float | None) -> str:
    if seconds is None or seconds < 0:
        return "--:--:--"

    seconds = int(seconds)

    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def report(force: bool = False) -> None:
    global last_report_time
    global last_displayed_percent

    now = time.monotonic()
    elapsed = now - start_time

    fraction = (
        processed / total
        if total > 0
        else 0.0
    )

    fraction = min(max(fraction, 0.0), 1.0)

    displayed_percent = start_percent + int(
        fraction * (end_percent - start_percent)
    )

    rate = (
        processed / elapsed
        if elapsed > 0
        else 0.0
    )

    remaining_records = max(total - processed, 0)

    eta = (
        remaining_records / rate
        if rate > 0
        else None
    )

    should_report = (
        force
        or displayed_percent != last_displayed_percent
        or now - last_report_time >= 5
    )

    if not should_report:
        return

    print(
        (
            f"\r[{displayed_percent:3d}%] "
            f"processed {processed:,}/{total:,} | "
            f"selected {selected:,} | "
            f"elapsed {format_seconds(elapsed)} | "
            f"ETA {format_seconds(eta)}"
        ),
        file=sys.stderr,
        end="",
        flush=True,
    )

    last_report_time = now
    last_displayed_percent = displayed_percent


for line in sys.stdin:
    if line.startswith("@"):
        sys.stdout.write(line)
        continue

    processed += 1

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 11:
        report()
        continue

    flag = int(fields[1])
    mapq = int(fields[4])

    if (
        flag & excluded_flags == 0
        and mapq >= minimum_mapq
    ):
        sys.stdout.write(line)
        selected += 1

    report()


sys.stdout.flush()
report(force=True)
print(file=sys.stderr)


if processed < total:
    completion = (
        processed / total * 100
        if total
        else 0
    )

    print(
        (
            "ERROR: remote stream ended early: "
            f"{processed:,}/{total:,} records "
            f"({completion:.3f}%)."
        ),
        file=sys.stderr,
        flush=True,
    )

    raise SystemExit(3)


print(
    (
        f"[{end_percent:3d}%] "
        f"chr21 stream completed successfully."
    ),
    file=sys.stderr,
    flush=True,
)
