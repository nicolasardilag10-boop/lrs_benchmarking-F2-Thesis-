#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from difflib import SequenceMatcher
from hashlib import sha256
from pathlib import Path
import datetime
import os
import re
import subprocess
import sys


ROOT = Path(
    subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"],
        text=True,
    ).strip()
).resolve()

os.chdir(ROOT)

STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

OUT = ROOT / "docs" / "cleanup_audit" / STAMP
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# Repository safety rules
# ============================================================

# Never content-hash large/biological/generated areas.
HASH_SKIP_PREFIXES = (
    ".git",
    ".snakemake",
    "fastq",
    "reference",
    "results",
    "assemblies",
    "happy_results",
    "truvari",
    "sawfish",
    "vcf_called",
    "alignment_analysis/alignments",
    "alignment_analysis/tmp",
    "SV aligners call/results",
    "assemblers/results",
    "assemblers/debug",
    "final_report_files/snakemake_aligners_benchmarking/results",
    "final_report_files/snakemake_assemblers_benchmarking/results",
)

# Even if zero-byte, don't automatically delete things from active areas.
DELETE_PROTECTED_PREFIXES = (
    ".git",
    ".snakemake",
    "fastq",
    "reference",
    "assemblers",
    "final_report_files",
    "alignment_analysis",
    "SV aligners call",
    "samples_try",
    "results",
    "assemblies",
    "happy_results",
    "truvari",
    "sawfish",
    "vcf_called",
    "workflow",
)

# Workflow marker files can legitimately be empty.
MARKER_SUFFIXES = (
    ".done",
    ".ok",
    ".running",
    ".failed",
    ".lock",
    ".touch",
    ".sentinel",
)

SOURCE_SUFFIXES = {
    ".smk",
    ".py",
    ".sh",
    ".yaml",
    ".yml",
    ".json",
    ".toml",
    ".ini",
    ".cfg",
    ".md",
    ".txt",
    ".tsv",
    ".csv",
    ".tex",
    ".r",
    ".rmd",
}

MAX_HASH_SIZE = 20 * 1024 * 1024
MAX_TEXT_COMPARE_SIZE = 2 * 1024 * 1024


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def starts_with_any(path_string: str, prefixes) -> bool:
    for prefix in prefixes:
        if (
            path_string == prefix
            or path_string.startswith(prefix + "/")
        ):
            return True
    return False


def is_source_candidate(path: Path) -> bool:
    name = path.name.lower()

    if name.startswith("snakefile"):
        return True

    return path.suffix.lower() in SOURCE_SUFFIXES


def sha(path: Path) -> str:
    h = sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def git_tracked() -> set[str]:
    result = subprocess.check_output(
        ["git", "ls-files", "-z"]
    )

    return {
        item.decode(errors="replace")
        for item in result.split(b"\0")
        if item
    }


TRACKED = git_tracked()


# ============================================================
# Collect files and directories
# ============================================================

all_files: list[Path] = []
all_dirs: list[Path] = []

for path in ROOT.rglob("*"):

    try:
        relative = rel(path)
    except ValueError:
        continue

    if relative.startswith(".git/"):
        continue

    if path.is_dir():
        all_dirs.append(path)

    elif path.is_file():
        all_files.append(path)


# ============================================================
# 1. Zero-byte files
# ============================================================

zero_files = []

safe_zero_files = []

for path in all_files:

    try:
        size = path.stat().st_size
    except OSError:
        continue

    if size != 0:
        continue

    relative = rel(path)

    tracked = relative in TRACKED

    marker = path.name.lower().endswith(MARKER_SUFFIXES)

    protected = starts_with_any(
        relative,
        DELETE_PROTECTED_PREFIXES,
    )

    zero_files.append(
        (
            relative,
            "TRACKED" if tracked else "UNTRACKED",
            "MARKER" if marker else "",
            "PROTECTED" if protected else "",
        )
    )

    if (
        not tracked
        and not marker
        and not protected
        and path.name != ".gitkeep"
    ):
        safe_zero_files.append(relative)


with (OUT / "zero_byte_files.tsv").open("w") as out:
    out.write("path\tgit_status\tmarker\tprotection\n")

    for row in sorted(zero_files):
        out.write("\t".join(row) + "\n")


with (OUT / "safe_zero_byte_candidates.txt").open("w") as out:
    for path in sorted(safe_zero_files):
        out.write(path + "\n")


# ============================================================
# 2. Empty directories
#
# REPORT ONLY.
# Do not delete automatically.
# ============================================================

empty_dirs = []

for path in all_dirs:

    try:
        if not any(path.iterdir()):
            empty_dirs.append(rel(path))
    except OSError:
        pass


with (OUT / "empty_directories.txt").open("w") as out:
    for path in sorted(empty_dirs):
        out.write(path + "\n")


# ============================================================
# 3. Exact duplicate SOURCE files
# ============================================================

hash_groups: dict[str, list[Path]] = defaultdict(list)

hashed_count = 0

for path in all_files:

    relative = rel(path)

    if starts_with_any(relative, HASH_SKIP_PREFIXES):
        continue

    if not is_source_candidate(path):
        continue

    try:
        size = path.stat().st_size
    except OSError:
        continue

    if size == 0 or size > MAX_HASH_SIZE:
        continue

    try:
        digest = sha(path)
    except OSError:
        continue

    hash_groups[digest].append(path)
    hashed_count += 1


duplicate_groups = [
    (digest, paths)
    for digest, paths in hash_groups.items()
    if len(paths) > 1
]


with (OUT / "exact_duplicate_files.txt").open("w") as out:

    for number, (digest, paths) in enumerate(
        sorted(
            duplicate_groups,
            key=lambda item: sorted(rel(p) for p in item[1])[0],
        ),
        start=1,
    ):
        out.write(
            f"==================================================\n"
        )
        out.write(f"GROUP {number}\n")
        out.write(f"SHA256: {digest}\n")
        out.write(
            f"SIZE: {paths[0].stat().st_size} bytes\n"
        )

        for path in sorted(paths, key=lambda p: rel(p)):
            flag = []

            if rel(path) in TRACKED:
                flag.append("tracked")
            else:
                flag.append("untracked")

            if re.search(
                r"(before|backup|\.bak|old|copy)",
                path.name,
                flags=re.IGNORECASE,
            ):
                flag.append("backup-like")

            out.write(
                f"  {rel(path)}"
                f"\t[{', '.join(flag)}]\n"
            )

        out.write("\n")


# ============================================================
# 4. Same filename in multiple locations
# ============================================================

same_names: dict[str, list[Path]] = defaultdict(list)

for path in all_files:

    relative = rel(path)

    if starts_with_any(relative, HASH_SKIP_PREFIXES):
        continue

    if not is_source_candidate(path):
        continue

    same_names[path.name.lower()].append(path)


with (OUT / "same_filename_multiple_locations.txt").open("w") as out:

    for name, paths in sorted(same_names.items()):

        if len(paths) < 2:
            continue

        out.write(
            "==================================================\n"
        )
        out.write(f"{name}\n")

        hashes = {}

        for path in paths:

            try:
                if path.stat().st_size <= MAX_HASH_SIZE:
                    hashes[path] = sha(path)
                else:
                    hashes[path] = "TOO_LARGE"
            except OSError:
                hashes[path] = "ERROR"

        unique_hashes = set(hashes.values())

        if len(unique_hashes) == 1:
            status = "IDENTICAL_CONTENT"
        else:
            status = "DIFFERENT_CONTENT"

        out.write(f"STATUS: {status}\n")

        for path in sorted(paths, key=lambda p: rel(p)):
            out.write(
                f"  {rel(path)}"
                f"\t{hashes[path]}\n"
            )

        out.write("\n")


# ============================================================
# 5. Backup/version families
# ============================================================

def normalize_backup_name(name: str) -> str:

    value = name.lower()

    # alignment_metrics.before_barplots.py
    value = re.sub(
        r"([._-])before([._-]).*",
        "",
        value,
    )

    # config.yaml.backup_2026...
    value = re.sub(
        r"([._-])backup([._-]?).*",
        "",
        value,
    )

    # file.bak
    value = re.sub(
        r"\.bak$",
        "",
        value,
    )

    # file.old
    value = re.sub(
        r"\.old$",
        "",
        value,
    )

    # copy variants
    value = re.sub(
        r"([._-])copy([._-]?\d*)?",
        "",
        value,
    )

    # timestamps
    value = re.sub(
        r"[._-]20\d{6}[_-]\d{6}",
        "",
        value,
    )

    return value


families: dict[str, list[Path]] = defaultdict(list)

for path in all_files:

    relative = rel(path)

    if starts_with_any(relative, HASH_SKIP_PREFIXES):
        continue

    if not is_source_candidate(path):
        continue

    normalized = normalize_backup_name(path.name)

    if normalized != path.name.lower():
        families[normalized].append(path)

        # Look for potential current sibling.
        sibling = path.parent / normalized

        if sibling.exists() and sibling.is_file():
            families[normalized].append(sibling)


with (OUT / "backup_version_families.txt").open("w") as out:

    for family, paths in sorted(families.items()):

        paths = sorted(
            set(paths),
            key=lambda p: rel(p),
        )

        if len(paths) < 2:
            continue

        out.write(
            "==================================================\n"
        )
        out.write(f"FAMILY: {family}\n")

        for path in paths:
            out.write(f"  {rel(path)}\n")

        out.write("\n")


# ============================================================
# 6. Similar TOP-LEVEL directory names
# ============================================================

top_dirs = sorted(
    [
        path
        for path in ROOT.iterdir()
        if path.is_dir()
        and path.name not in {".git", ".snakemake"}
    ],
    key=lambda p: p.name.lower(),
)


similar_top_dirs = []

for i, left in enumerate(top_dirs):

    for right in top_dirs[i + 1:]:

        ratio = SequenceMatcher(
            None,
            left.name.lower(),
            right.name.lower(),
        ).ratio()

        if ratio >= 0.72:
            similar_top_dirs.append(
                (
                    ratio,
                    left.name,
                    right.name,
                )
            )


with (OUT / "similar_top_level_directories.tsv").open("w") as out:
    out.write("similarity\tdirectory_1\tdirectory_2\n")

    for ratio, left, right in sorted(
        similar_top_dirs,
        reverse=True,
    ):
        out.write(
            f"{ratio:.3f}\t{left}\t{right}\n"
        )


# ============================================================
# 7. Near-duplicate text content inside backup families
# ============================================================

with (OUT / "backup_content_similarity.tsv").open("w") as out:

    out.write(
        "family\tfile_1\tfile_2\tsimilarity\n"
    )

    for family, paths in sorted(families.items()):

        paths = sorted(
            set(paths),
            key=lambda p: rel(p),
        )

        if len(paths) < 2:
            continue

        for i, left in enumerate(paths):

            for right in paths[i + 1:]:

                try:
                    if (
                        left.stat().st_size
                        > MAX_TEXT_COMPARE_SIZE
                        or right.stat().st_size
                        > MAX_TEXT_COMPARE_SIZE
                    ):
                        continue

                    left_text = left.read_text(
                        errors="replace"
                    )

                    right_text = right.read_text(
                        errors="replace"
                    )

                except OSError:
                    continue

                ratio = SequenceMatcher(
                    None,
                    left_text,
                    right_text,
                ).ratio()

                out.write(
                    f"{family}\t"
                    f"{rel(left)}\t"
                    f"{rel(right)}\t"
                    f"{ratio:.4f}\n"
                )


# ============================================================
# 8. Summary
# ============================================================

summary = f"""
LRS BENCHMARKING CLEANUP AUDIT

Repository:
  {ROOT}

Generated:
  {STAMP}

Source/config/doc files hashed:
  {hashed_count}

Exact duplicate groups:
  {len(duplicate_groups)}

Zero-byte files:
  {len(zero_files)}

Safe untracked zero-byte candidates:
  {len(safe_zero_files)}

Empty directories:
  {len(empty_dirs)}

Similar top-level directory pairs:
  {len(similar_top_dirs)}

Reports:
  {OUT}
""".strip()


(OUT / "SUMMARY.txt").write_text(
    summary + "\n"
)

print()
print("=" * 70)
print(summary)
print("=" * 70)


# ============================================================
# Optional zero-byte deletion
# ============================================================

if "--delete-safe-empty" in sys.argv:

    print()
    print("SAFE ZERO-BYTE CLEANUP")
    print("-" * 70)

    deleted = 0

    for relative in safe_zero_files:

        path = ROOT / relative

        # Revalidate at deletion time.
        if not path.is_file():
            continue

        if path.stat().st_size != 0:
            continue

        if relative in TRACKED:
            continue

        if starts_with_any(
            relative,
            DELETE_PROTECTED_PREFIXES,
        ):
            continue

        if path.name.lower().endswith(
            MARKER_SUFFIXES
        ):
            continue

        print(f"DELETE: {relative}")
        path.unlink()

        deleted += 1

    print()
    print(
        f"Deleted {deleted} safe untracked "
        "zero-byte file(s)."
    )
