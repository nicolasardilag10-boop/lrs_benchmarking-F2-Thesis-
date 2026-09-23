from __future__ import annotations

import sys
from pathlib import Path


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: stream_file_progress.py "
        "FILE START_PERCENT END_PERCENT"
    )

path = Path(sys.argv[1])
start_percent = int(sys.argv[2])
end_percent = int(sys.argv[3])

file_size = path.stat().st_size
processed = 0
last_percent = -1
block_size = 1024 * 1024

with path.open("rb") as handle:
    while True:
        block = handle.read(block_size)

        if not block:
            break

        sys.stdout.buffer.write(block)
        processed += len(block)

        fraction = (
            processed / file_size
            if file_size
            else 1.0
        )

        percent = start_percent + int(
            fraction * (end_percent - start_percent)
        )

        percent = min(percent, end_percent)

        if percent != last_percent:
            print(
                (
                    f"\r[{percent:3d}%] "
                    f"BAM processed: "
                    f"{processed:,}/{file_size:,} bytes"
                ),
                file=sys.stderr,
                end="",
                flush=True,
            )

            last_percent = percent

sys.stdout.buffer.flush()

print(
    (
        f"\r[{end_percent:3d}%] "
        f"BAM streaming completed: "
        f"{processed:,}/{file_size:,} bytes"
    ),
    file=sys.stderr,
    flush=True,
)
