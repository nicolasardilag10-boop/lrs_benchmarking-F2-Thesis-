#!/usr/bin/env python3
"""
Recover only missing insertion/deletion metrics in an existing alignment summary.

Purpose
-------
This script is intentionally lightweight. It DOES NOT:
  - rerun minimap2, pbmm2, VACMap, or VG Giraffe
  - recompute soft clipping
  - recompute coverage/breadth
  - recompute mapping/error metrics
  - modify already-populated indel values

It reads an existing TSV, identifies rows whose indel metrics are missing,
and tries, in this order:

  1) Parse a matching full samtools *.cram.stats file if one exists.
  2) Otherwise, locate the matching CRAM and run `samtools stats` on that
     existing CRAM using a pinned Docker image and the reference FASTA.

Recovered metrics:
  insertion_events
  deletion_events
  inserted_bases
  deleted_bases
  insertion_events_per_100kb
  deletion_events_per_100kb

Normalization:
  events_per_100kb = events / bases_mapped_cigar * 100000

All other columns and existing values in the input TSV are preserved exactly.
A separate recovery report is written for provenance.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

NA_VALUES = {"", "NA", "N/A", "NAN", "NONE"}
SAMTOOLS_IMAGE_DEFAULT = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"

INDEL_COLUMNS = [
    "insertion_events",
    "deletion_events",
    "inserted_bases",
    "deleted_bases",
    "insertion_events_per_100kb",
    "deletion_events_per_100kb",
]

METHOD_SPECS = {
    "mm2-ont": ("mm2-ont", "minimap2"),
    "mm2-pb": ("mm2-pb", "minimap2"),
    "pbmm2-ont": ("pbmm2-ont", "pbmm2"),
    "pbmm2-subread": ("pbmm2-ont", "pbmm2"),
    "pbmm2-pb": ("pbmm2-pb", "pbmm2"),
    "pbmm2-ccs": ("pbmm2-pb", "pbmm2"),
    "vacmap-ont": ("vacmap-ont", "VACMap"),
    "vacmap-pb": ("vacmap-pb", "VACMap"),
    "vg-ont": ("vg-ont", "VG Giraffe"),
    "vg-pb": ("vg-pb", "VG Giraffe"),
}

METHOD_TAGS = sorted(METHOD_SPECS, key=len, reverse=True)

DATASET_RE = re.compile(
    r"(?P<sample>HG00[234])[._](?P<technology>ont|pb)[._](?P<coverage>1k|30x)(?=[._])",
    re.IGNORECASE,
)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def is_missing(value) -> bool:
    if value is None:
        return True
    return str(value).strip().upper() in NA_VALUES


def num(value) -> Optional[float]:
    if is_missing(value):
        return None
    try:
        return float(str(value).strip().replace(",", ""))
    except ValueError:
        return None


def fmt(value, digits: int = 6) -> str:
    if value is None:
        return "NA"
    value = float(value)
    if abs(value - round(value)) < 1e-12:
        return str(int(round(value)))
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def per100k(events, cigar_bases) -> Optional[float]:
    events = num(events)
    cigar_bases = num(cigar_bases)
    if events is None or cigar_bases in (None, 0):
        return None
    return events / cigar_bases * 100000.0


def normalize_technology(value: str) -> Optional[str]:
    v = str(value).strip().lower()
    if v in {"ont", "nanopore", "oxford nanopore"}:
        return "ONT"
    if v in {"pb", "pacbio", "pacbio hifi", "hifi"}:
        return "PacBio"
    return None


def infer_mapper_tag(row: Dict[str, str]) -> Optional[str]:
    """
    Prefer mapper_tag/configuration already present in the table.
    Otherwise infer it from aligner + technology.
    """
    for key in ("mapper_tag", "configuration"):
        raw = str(row.get(key, "")).strip().lower()
        if raw in METHOD_SPECS:
            return METHOD_SPECS[raw][0]

    aligner = str(row.get("aligner", "")).strip().lower()
    tech = normalize_technology(row.get("read_technology", ""))
    if tech is None:
        return None

    suffix = "ont" if tech == "ONT" else "pb"

    if aligner == "minimap2":
        return f"mm2-{suffix}"
    if aligner == "pbmm2":
        return f"pbmm2-{suffix}"
    if aligner in {"vacmap", "vacmap"}:
        return f"vacmap-{suffix}"
    if aligner in {"vg giraffe", "vg", "giraffe"}:
        return f"vg-{suffix}"

    return None


def row_identity(row: Dict[str, str]) -> Optional[Tuple[str, str, str, str]]:
    sample = str(row.get("sample", "")).strip().upper()
    tech = normalize_technology(row.get("read_technology", ""))
    coverage = str(row.get("coverage", "30x")).strip().lower() or "30x"
    mapper_tag = infer_mapper_tag(row)

    if not sample or tech is None or mapper_tag is None:
        return None

    tech_token = "ont" if tech == "ONT" else "pb"
    return sample, tech_token, coverage, mapper_tag


def identity_from_filename(path: Path) -> Optional[Tuple[str, str, str, str]]:
    m = DATASET_RE.search(path.name)
    if not m:
        return None

    sample = m.group("sample").upper()
    tech = m.group("technology").lower()
    coverage = m.group("coverage").lower()

    lower = path.name.lower()
    tag = None
    for candidate in METHOD_TAGS:
        if re.search(rf"(?:^|\.){re.escape(candidate)}(?=\.|$)", lower):
            tag = METHOD_SPECS[candidate][0]
            break

    if tag is None:
        return None

    return sample, tech, coverage, tag


def parse_id_records(lines: Iterable[str]) -> Dict[str, Optional[float]]:
    """
    Parse samtools stats ID lines.

    ID columns:
      ID <indel_length> <number_of_insertions> <number_of_deletions>
    """
    insertion_events = 0
    deletion_events = 0
    inserted_bases = 0
    deleted_bases = 0
    saw_id = False

    for line in lines:
        if not line.startswith("ID\t"):
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue

        length = num(fields[1])
        n_ins = num(fields[2])
        n_del = num(fields[3])

        if length is None or n_ins is None or n_del is None:
            continue

        saw_id = True
        length = int(length)
        n_ins = int(n_ins)
        n_del = int(n_del)

        insertion_events += n_ins
        deletion_events += n_del
        inserted_bases += length * n_ins
        deleted_bases += length * n_del

    if not saw_id:
        return {}

    return {
        "insertion_events": insertion_events,
        "deletion_events": deletion_events,
        "inserted_bases": inserted_bases,
        "deleted_bases": deleted_bases,
    }


def parse_full_stats(path: Path) -> Dict[str, Optional[float]]:
    """
    Parse ID records and, as a fallback denominator, SN 'bases mapped (cigar)'.

    The file is streamed because whole-genome stats files can be hundreds of
    megabytes; memory usage therefore stays essentially constant.
    """
    insertion_events = 0
    deletion_events = 0
    inserted_bases = 0
    deleted_bases = 0
    saw_id = False
    cigar_bases = None

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                length = num(fields[1])
                n_ins = num(fields[2])
                n_del = num(fields[3])
                if length is None or n_ins is None or n_del is None:
                    continue
                saw_id = True
                insertion_events += int(n_ins)
                deletion_events += int(n_del)
                inserted_bases += int(length) * int(n_ins)
                deleted_bases += int(length) * int(n_del)
            elif line.startswith("SN\t"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                label = fields[1].strip().rstrip(":").lower()
                if label == "bases mapped (cigar)":
                    cigar_bases = num(fields[2])

    indels: Dict[str, Optional[float]] = {}
    if saw_id:
        indels.update(
            {
                "insertion_events": insertion_events,
                "deletion_events": deletion_events,
                "inserted_bases": inserted_bases,
                "deleted_bases": deleted_bases,
            }
        )

    if cigar_bases is not None:
        indels["_cigar_bases_from_stats"] = cigar_bases

    return indels


def discover_full_stats(
    project: Path,
    stats_dir: Optional[Path] = None,
) -> Dict[Tuple[str, str, str, str], Path]:
    """
    Index full stats files once. SN-only extracts are explicitly excluded.
    """
    roots = (
        [stats_dir]
        if stats_dir is not None
        else [
            project / "cram",
            project / "alignment_analysis" / "tables",
            project / "samtools_stats_30x_Christian",
        ]
    )

    result: Dict[Tuple[str, str, str, str], Path] = {}

    for root in roots:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if not path.is_file():
                continue

            name = path.name
            if name.endswith(".cram.stats.SN.txt"):
                continue

            if not (
                name.endswith(".cram.stats")
                or name.endswith(".samtools_stats.txt")
                or name.endswith(".stats.txt")
                or name.endswith(".stats")
            ):
                continue

            ident = identity_from_filename(path)
            if ident is None:
                continue

            # Prefer files under project/cram because they are closest to
            # the original alignment workflow output.
            current = result.get(ident)
            if current is None:
                result[ident] = path
            elif str(path).startswith(str(project / "cram")) and not str(current).startswith(
                str(project / "cram")
            ):
                result[ident] = path

    return result


def discover_crams(cram_root: Path) -> Dict[Tuple[str, str, str, str], Path]:
    """
    Recursively index CRAM files below an explicitly selected directory.
    """
    result: Dict[Tuple[str, str, str, str], Path] = {}

    for path in cram_root.rglob("*.cram"):
        if not path.is_file():
            continue

        ident = identity_from_filename(path)
        if ident is None:
            continue

        current = result.get(ident)
        if current is None:
            result[ident] = path
        else:
            # Prefer the conventional *.hg38.<tag>.cram form.
            preferred_new = ".hg38." in path.name.lower()
            preferred_old = ".hg38." in current.name.lower()
            if preferred_new and not preferred_old:
                result[ident] = path

    return result


def docker_base(
    project: Path,
    reference: Path,
    cram: Path,
    image: str,
) -> Tuple[List[str], str, str]:
    """
    Return the Docker command prefix plus container paths for the reference
    and CRAM. Inputs outside the project are mounted separately, read-only.
    """
    project = project.resolve()
    reference = reference.resolve()
    cram = cram.resolve()

    mounts = ["-v", f"{project}:/work:ro"]

    try:
        ref_rel = reference.relative_to(project)
        ref_in = f"/work/{ref_rel}"
    except ValueError:
        mounts.extend(["-v", f"{reference.parent}:/reference:ro"])
        ref_in = f"/reference/{reference.name}"

    try:
        cram_rel = cram.relative_to(project)
        cram_in = f"/work/{cram_rel}"
    except ValueError:
        mounts.extend(["-v", f"{cram.parent}:/cram_input:ro"])
        cram_in = f"/cram_input/{cram.name}"

    cmd = [
        "docker",
        "run",
        "--rm",
        "-u",
        f"{os.getuid()}:{os.getgid()}",
        *mounts,
        "-w",
        "/work",
        image,
    ]

    return cmd, ref_in, cram_in


def check_docker() -> None:
    try:
        subprocess.run(
            ["docker", "version", "--format", "{{.Server.Version}}"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "Docker is required to recover indels from CRAM files, but Docker "
            "is not available on this server."
        ) from exc


def recover_indels_from_cram(
    project: Path,
    cram: Path,
    reference: Path,
    image: str,
    threads: int,
) -> Dict[str, Optional[float]]:
    """
    Run only samtools stats on an existing CRAM.

    No alignment, clipping, or coverage computation is performed.
    """
    base, ref_in, cram_in = docker_base(project, reference, cram, image)

    cmd = base + [
        "samtools",
        "stats",
        "-@",
        str(threads),
        "--reference",
        ref_in,
        "--remove-overlaps",
        cram_in,
    ]

    proc = subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    result = parse_id_records(proc.stdout.splitlines())

    # Also extract bases mapped (cigar) from the same samtools output as a
    # denominator fallback if the input table lacks it.
    for line in proc.stdout.splitlines():
        if not line.startswith("SN\t"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        label = fields[1].strip().rstrip(":").lower()
        if label == "bases mapped (cigar)":
            result["_cigar_bases_from_stats"] = num(fields[2])
            break

    return result


def ensure_indel_columns(fieldnames: List[str]) -> List[str]:
    out = list(fieldnames)
    for col in INDEL_COLUMNS:
        if col not in out:
            out.append(col)
    return out


def read_table(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in {path}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)

    required = {"sample", "read_technology", "aligner", "bases_mapped_cigar"}
    missing = required - set(fieldnames)
    if missing:
        raise ValueError(
            "Input table is missing required columns: " + ", ".join(sorted(missing))
        )

    return rows, ensure_indel_columns(fieldnames)


def row_needs_recovery(row: Dict[str, str]) -> bool:
    return any(is_missing(row.get(col)) for col in INDEL_COLUMNS)


def fill_indels(
    row: Dict[str, str],
    raw: Dict[str, Optional[float]],
) -> bool:
    """
    Fill only missing indel values. Existing values are never overwritten.
    Returns True if at least one value was added.
    """
    if not raw:
        return False

    cigar_bases = num(row.get("bases_mapped_cigar"))
    if cigar_bases in (None, 0):
        cigar_bases = num(raw.get("_cigar_bases_from_stats"))

    values = {
        "insertion_events": raw.get("insertion_events"),
        "deletion_events": raw.get("deletion_events"),
        "inserted_bases": raw.get("inserted_bases"),
        "deleted_bases": raw.get("deleted_bases"),
    }

    values["insertion_events_per_100kb"] = per100k(
        values["insertion_events"], cigar_bases
    )
    values["deletion_events_per_100kb"] = per100k(
        values["deletion_events"], cigar_bases
    )

    changed = False
    for col, value in values.items():
        if is_missing(row.get(col)) and value is not None:
            row[col] = fmt(value)
            changed = True

    return changed


def write_table(path: Path, rows: List[Dict[str, str]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, records: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "read_technology",
        "aligner",
        "mapper_tag",
        "status",
        "source_type",
        "source_file",
    ]

    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Fill only missing insertion/deletion metrics in an existing "
            "alignment summary TSV."
        )
    )

    p.add_argument(
        "--project",
        type=Path,
        required=True,
        help="Root of lrs_benchmarking on the server.",
    )

    p.add_argument(
        "--input-table",
        type=Path,
        required=True,
        help="Existing trusted alignment summary TSV to update.",
    )

    p.add_argument(
        "--cram-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing input CRAM files (searched recursively). "
            "Defaults to PROJECT/cram when that directory exists; otherwise "
            "defaults to PROJECT for backward compatibility."
        ),
    )

    p.add_argument(
        "--stats-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing full .cram.stats files (searched "
            "recursively). Defaults to "
            "PROJECT/alignment_analysis/statistics_cram_files when present."
        ),
    )

    p.add_argument(
        "--reference-fasta",
        type=Path,
        default=None,
        help=(
            "The same GRCh38 FASTA used for the CRAM alignments. Required "
            "only when recovery must fall back to processing CRAM files."
        ),
    )

    p.add_argument(
        "--out",
        type=Path,
        required=True,
        help="New TSV. The input file is never overwritten.",
    )

    p.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional provenance/recovery report TSV.",
    )

    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads used by samtools stats for rows that require CRAM recovery.",
    )

    p.add_argument(
        "--samtools-image",
        default=SAMTOOLS_IMAGE_DEFAULT,
        help="Pinned samtools Docker image.",
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    project = args.project.expanduser().resolve()
    input_table = args.input_table.expanduser().resolve()
    reference = (
        args.reference_fasta.expanduser().resolve()
        if args.reference_fasta is not None
        else None
    )
    output = args.out.expanduser().resolve()
    if args.stats_dir is not None:
        stats_dir = args.stats_dir.expanduser().resolve()
    else:
        conventional_stats_dir = (
            project / "alignment_analysis" / "statistics_cram_files"
        )
        stats_dir = conventional_stats_dir if conventional_stats_dir.is_dir() else None
    if args.cram_dir is not None:
        cram_dir = args.cram_dir.expanduser().resolve()
    else:
        conventional_cram_dir = project / "cram"
        cram_dir = conventional_cram_dir if conventional_cram_dir.is_dir() else project

    report = (
        args.report.expanduser().resolve()
        if args.report
        else output.with_name(output.stem + ".indel_recovery_report.tsv")
    )

    if not project.exists():
        print(f"ERROR: project directory not found: {project}", file=sys.stderr)
        return 1

    if not input_table.is_file():
        print(f"ERROR: input TSV not found: {input_table}", file=sys.stderr)
        return 1

    if reference is not None and not reference.is_file():
        print(f"ERROR: reference FASTA not found: {reference}", file=sys.stderr)
        return 1

    if stats_dir is not None and not stats_dir.is_dir():
        print(f"ERROR: statistics directory not found: {stats_dir}", file=sys.stderr)
        return 1

    if not cram_dir.is_dir():
        print(f"ERROR: CRAM directory not found: {cram_dir}", file=sys.stderr)
        return 1

    try:
        rows, fieldnames = read_table(input_table)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    missing_before = sum(row_needs_recovery(row) for row in rows)

    print(f"Rows in input table: {len(rows)}")
    print(f"Rows with at least one missing indel metric: {missing_before}")

    if stats_dir is not None:
        print(f"Indexing samtools stats files under {stats_dir}...", file=sys.stderr)
    else:
        print("Indexing existing full samtools stats files...", file=sys.stderr)
    full_stats = discover_full_stats(project, stats_dir)

    print(f"Indexing CRAM files under {cram_dir}...", file=sys.stderr)
    crams = discover_crams(cram_dir)

    # Only require Docker if at least one row may need CRAM recovery.
    need_docker = False
    for row in rows:
        if not row_needs_recovery(row):
            continue
        ident = row_identity(row)
        if ident is not None and ident not in full_stats and ident in crams:
            need_docker = True
            break

    if need_docker:
        if reference is None:
            print(
                "ERROR: --reference-fasta is required because one or more "
                "rows need CRAM fallback processing.",
                file=sys.stderr,
            )
            return 1
        try:
            check_docker()
        except RuntimeError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1

    records: List[Dict[str, str]] = []
    recovered_rows = 0

    for i, row in enumerate(rows, start=1):
        ident = row_identity(row)

        base_record = {
            "sample": str(row.get("sample", "")),
            "read_technology": str(row.get("read_technology", "")),
            "aligner": str(row.get("aligner", "")),
            "mapper_tag": infer_mapper_tag(row) or "",
            "status": "",
            "source_type": "",
            "source_file": "",
        }

        if not row_needs_recovery(row):
            base_record["status"] = "already_complete"
            records.append(base_record)
            continue

        if ident is None:
            base_record["status"] = "unrecognized_row_identity"
            records.append(base_record)
            warn(
                f"Could not identify row {i}: "
                f"{row.get('sample')} {row.get('read_technology')} {row.get('aligner')}"
            )
            continue

        sample, tech, coverage, tag = ident
        print(
            f"[{i}/{len(rows)}] Recovering {sample} {tech} {tag} ...",
            file=sys.stderr,
        )

        recovered = False

        # 1) Fast path: parse an existing complete stats file.
        stats_path = full_stats.get(ident)
        if stats_path is not None:
            try:
                raw = parse_full_stats(stats_path)
                recovered = fill_indels(row, raw)
            except OSError as exc:
                warn(f"Could not read {stats_path}: {exc}")
                raw = {}

            if recovered:
                recovered_rows += 1
                base_record["status"] = "recovered"
                base_record["source_type"] = "existing_full_samtools_stats"
                base_record["source_file"] = str(stats_path)
                records.append(base_record)
                continue

        # 2) Slow path: run only samtools stats on the existing CRAM.
        cram = crams.get(ident)
        if cram is not None:
            try:
                raw = recover_indels_from_cram(
                    project=project,
                    cram=cram,
                    reference=reference,
                    image=args.samtools_image,
                    threads=args.threads,
                )
                recovered = fill_indels(row, raw)
            except subprocess.CalledProcessError as exc:
                stderr = exc.stderr.strip() if exc.stderr else str(exc)
                warn(f"samtools stats failed for {cram}: {stderr}")
                recovered = False
            except RuntimeError as exc:
                warn(str(exc))
                recovered = False

            if recovered:
                recovered_rows += 1
                base_record["status"] = "recovered"
                base_record["source_type"] = "samtools_stats_from_existing_cram"
                base_record["source_file"] = str(cram)
            else:
                base_record["status"] = "no_ID_records_or_missing_denominator"
                base_record["source_type"] = "samtools_stats_from_existing_cram"
                base_record["source_file"] = str(cram)

            records.append(base_record)
            continue

        base_record["status"] = "missing_full_stats_and_cram"
        records.append(base_record)

    missing_after = sum(row_needs_recovery(row) for row in rows)

    try:
        write_table(output, rows, fieldnames)
        write_report(report, records)
    except OSError as exc:
        print(f"ERROR writing output: {exc}", file=sys.stderr)
        return 1

    print()
    print("========================================")
    print("INDEL RECOVERY SUMMARY")
    print("========================================")
    print(f"Input rows:                         {len(rows)}")
    print(f"Rows missing indels before:         {missing_before}")
    print(f"Rows recovered in this run:         {recovered_rows}")
    print(f"Rows still missing indels after:    {missing_after}")
    print()
    print(f"Updated table: {output}")
    print(f"Recovery report: {report}")
    print()
    print(
        "Only missing indel columns were filled. "
        "All other columns and existing values were preserved."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
