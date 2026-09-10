#!/usr/bin/env python3
"""
Build a comparable long-read alignment summary for:
  minimap2, pbmm2, VACMap and VG Giraffe
on ONT and PacBio HiFi data.

Core metrics are parsed from existing samtools stats files.
The repository's committed samtools_stats_30x_Christian/*.cram.stats.SN.txt
files are supported as a fallback when the full server-side *.cram.stats
reports are not present.
With --reference-fasta, CRAM-derived soft clipping and coverage metrics are
calculated automatically using a pinned samtools Docker image.

Runtime/RAM/threads and original aligner version/image are loaded from recorded
run metadata when available. They cannot be reconstructed after a mapper run
if they were never recorded.

Missing values are NA, never zero.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

NA = "NA"
SAMTOOLS_IMAGE_DEFAULT = "quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
EXCLUDE_FLAGS = 3844  # UNMAP|SECONDARY|QCFAIL|DUP|SUPPLEMENTARY

SUPPORTED_SUFFIXES = (
    ".samtools_stats.txt",
    ".cram.stats.SN.txt",
    ".cram.stats",
    ".stats.txt",
    ".stats",
)

METHOD_SPECS = {
    "mm2-ont":       ("mm2-ont",    "minimap2",   "map-ont"),
    "mm2-pb":        ("mm2-pb",     "minimap2",   "map-hifi"),
    "pbmm2-ont":     ("pbmm2-ont",  "pbmm2",      "SUBREAD"),
    "pbmm2-subread": ("pbmm2-ont",  "pbmm2",      "SUBREAD"),   # legacy
    "pbmm2-pb":      ("pbmm2-pb",   "pbmm2",      "CCS/HIFI"),
    "pbmm2-ccs":     ("pbmm2-pb",   "pbmm2",      "CCS/HIFI"),  # legacy
    "vacmap-ont":    ("vacmap-ont", "VACMap",     "vacmap-ont"),
    "vacmap-pb":     ("vacmap-pb",  "VACMap",     "vacmap-pb"),
    "vg-ont":        ("vg-ont",     "VG Giraffe", "vg-ont"),
    "vg-pb":         ("vg-pb",      "VG Giraffe", "vg-pb"),
}
METHOD_TAGS = sorted(METHOD_SPECS, key=len, reverse=True)

RUN_METADATA_COLUMNS = [
    "runtime_seconds",
    "peak_ram_mb",
    "threads",
    "aligner_version",
    "aligner_docker_image",
]

OUTPUT_COLUMNS = [
    # identification
    "sample", "read_technology", "coverage", "dataset", "reference",
    "aligner", "preset", "mapper_tag", "configuration",
    "statistics_file", "cram_file",

    # mapping
    "raw_total_sequences", "reads_mapped", "reads_unmapped",
    "mapped_reads_percent",

    # mapping quality
    "reads_mq0", "reads_mq0_percent", "mapq_mean", "mapq_median",

    # alignment structure
    "secondary_alignments", "secondary_alignments_per_100_mapped_reads",
    "supplementary_alignments",
    "supplementary_alignments_per_100_mapped_reads",

    # bases / disagreement
    "total_length", "bases_mapped", "bases_mapped_cigar",
    "mapped_bases_percent", "mismatches", "error_rate", "error_percent",

    # indels
    "insertion_events", "deletion_events", "inserted_bases", "deleted_bases",
    "insertion_events_per_100kb", "deletion_events_per_100kb",

    # clipping
    "soft_clipped_bases", "soft_clipped_percent",

    # read lengths
    "average_length", "maximum_length",

    # coverage
    "mean_coverage", "median_coverage",
    "breadth_1x_percent", "breadth_10x_percent",
    "breadth_20x_percent", "breadth_30x_percent",

    # performance / reproducibility
    "runtime_seconds", "peak_ram_mb", "threads",
    "aligner_version", "samtools_version",
    "aligner_docker_image", "samtools_docker_image",

    # completeness
    "metrics_complete", "missing_metrics",
]

DATASET_RE = re.compile(
    r"(?P<sample>HG00[234])[._](?P<technology>ont|pb)[._](?P<coverage>1k|30x)(?=[._])",
    re.IGNORECASE,
)

CORE_REQUIRED = [
    "raw_total_sequences", "reads_mapped", "reads_unmapped",
    "mapped_reads_percent", "reads_mq0", "reads_mq0_percent",
    "mapq_mean", "mapq_median", "secondary_alignments",
    "supplementary_alignments", "total_length", "bases_mapped",
    "bases_mapped_cigar", "mapped_bases_percent", "mismatches",
    "error_rate", "error_percent", "insertion_events", "deletion_events",
    "inserted_bases", "deleted_bases", "insertion_events_per_100kb",
    "deletion_events_per_100kb", "average_length", "maximum_length",
]

FULL_REQUIRED = CORE_REQUIRED + [
    "soft_clipped_bases", "soft_clipped_percent",
    "mean_coverage", "median_coverage",
    "breadth_1x_percent", "breadth_10x_percent",
    "breadth_20x_percent", "breadth_30x_percent",
    "runtime_seconds", "peak_ram_mb", "threads",
    "aligner_version", "samtools_version",
    "aligner_docker_image", "samtools_docker_image",
]


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def num(x) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip().replace(",", "")
    if not s or s.upper() == NA:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def fmt(x, digits: int = 6) -> str:
    if x is None:
        return NA
    if isinstance(x, str):
        return x if x.strip() else NA
    x = float(x)
    if abs(x - round(x)) < 1e-12:
        return str(int(round(x)))
    return f"{x:.{digits}f}".rstrip("0").rstrip(".")


def pct(a, b) -> Optional[float]:
    a, b = num(a), num(b)
    return None if a is None or b in (None, 0) else a / b * 100.0


def per100k(a, b) -> Optional[float]:
    a, b = num(a), num(b)
    return None if a is None or b in (None, 0) else a / b * 100000.0


def rel(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def identity_from_name(path: Path) -> Optional[Dict[str, str]]:
    m = DATASET_RE.search(path.name)
    if not m:
        return None

    sample = m.group("sample").upper()
    tech_token = m.group("technology").lower()
    coverage = m.group("coverage").lower()
    technology = "ONT" if tech_token == "ont" else "PacBio"

    found = None
    lower = path.name.lower()
    for tag in METHOD_TAGS:
        if re.search(rf"(?:^|\.){re.escape(tag)}(?=\.|$)", lower):
            found = tag
            break
    if not found:
        return None

    canonical, aligner, preset = METHOD_SPECS[found]
    expected = "-ont" if tech_token == "ont" else "-pb"
    if not canonical.endswith(expected):
        return None

    return {
        "sample": sample,
        "read_technology": technology,
        "coverage": coverage,
        "dataset": f"{sample}.{tech_token}.{coverage}",
        "reference": "hg38",
        "aligner": aligner,
        "preset": preset,
        "mapper_tag": canonical,
        "configuration": canonical,
    }


def row_key(d: Dict[str, str]) -> Tuple[str, str, str, str, str]:
    return (
        d["sample"], d["read_technology"], d["coverage"],
        d["reference"], d["mapper_tag"]
    )


def parse_stats(path: Path) -> Dict[str, object]:
    sn: Dict[str, float] = {}
    mapq: Dict[int, int] = defaultdict(int)
    ins_events = del_events = ins_bases = del_bases = 0
    saw_id = False

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue

            if f[0] == "SN" and len(f) >= 3:
                v = num(f[2])
                if v is not None:
                    sn[f[1].strip().rstrip(":").lower()] = v

            elif f[0] == "MAPQ" and len(f) >= 3:
                q, c = num(f[1]), num(f[2])
                if q is not None and c is not None:
                    mapq[int(q)] += int(c)

            elif f[0] == "ID" and len(f) >= 4:
                size, ni, nd = num(f[1]), num(f[2]), num(f[3])
                if size is not None and ni is not None and nd is not None:
                    saw_id = True
                    size, ni, nd = int(size), int(ni), int(nd)
                    ins_events += ni
                    del_events += nd
                    ins_bases += size * ni
                    del_bases += size * nd

    def S(name: str):
        return sn.get(name)

    raw = S("raw total sequences")
    mapped = S("reads mapped")
    unmapped = S("reads unmapped")
    if unmapped is None and raw is not None and mapped is not None:
        unmapped = max(raw - mapped, 0)

    mq0 = S("reads mq0")
    secondary = S("non-primary alignments")
    supplementary = S("supplementary alignments")
    total_len = S("total length")
    bases_mapped = S("bases mapped")
    cigar_bases = S("bases mapped (cigar)")
    mismatches = S("mismatches")
    error_rate = S("error rate")

    mean_mapq = median_mapq = None
    valid_q = sorted((q, c) for q, c in mapq.items() if q != 255 and c > 0)
    n = sum(c for _, c in valid_q)
    if n:
        mean_mapq = sum(q * c for q, c in valid_q) / n
        target1 = (n - 1) // 2
        target2 = n // 2

        def q_at(rank: int) -> int:
            acc = 0
            for q, c in valid_q:
                acc += c
                if acc > rank:
                    return q
            return valid_q[-1][0]

        median_mapq = (q_at(target1) + q_at(target2)) / 2.0

    return {
        "raw_total_sequences": raw,
        "reads_mapped": mapped,
        "reads_unmapped": unmapped,
        "mapped_reads_percent": pct(mapped, raw),
        "reads_mq0": mq0,
        "reads_mq0_percent": pct(mq0, mapped),
        "mapq_mean": mean_mapq,
        "mapq_median": median_mapq,
        "secondary_alignments": secondary,
        "secondary_alignments_per_100_mapped_reads": pct(secondary, mapped),
        "supplementary_alignments": supplementary,
        "supplementary_alignments_per_100_mapped_reads":
            pct(supplementary, mapped),
        "total_length": total_len,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": cigar_bases,
        "mapped_bases_percent": pct(bases_mapped, total_len),
        "mismatches": mismatches,
        "error_rate": error_rate,
        "error_percent": error_rate * 100 if error_rate is not None else None,
        "insertion_events": ins_events if saw_id else None,
        "deletion_events": del_events if saw_id else None,
        "inserted_bases": ins_bases if saw_id else None,
        "deleted_bases": del_bases if saw_id else None,
        "insertion_events_per_100kb":
            per100k(ins_events, cigar_bases) if saw_id else None,
        "deletion_events_per_100kb":
            per100k(del_events, cigar_bases) if saw_id else None,
        "average_length": S("average length"),
        "maximum_length": S("maximum length"),
    }


def discover_stats(project: Path) -> List[Path]:
    """Find full samtools reports and the committed 30x SN extracts."""
    found = set()
    roots = (
        project / "alignment_analysis" / "tables",
        project / "cram",
        project / "samtools_stats_30x_Christian",
    )
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if p.is_file() and any(p.name.endswith(s) for s in SUPPORTED_SUFFIXES):
                found.add(p.resolve())
    return sorted(found)


def stats_priority(p: Path, project: Path) -> Tuple[int, str]:
    """
    Prefer the complete server-side samtools report because it contains
    SN + MAPQ + ID sections. The committed *.SN.txt exports are fallback only.
    """
    in_cram = False
    try:
        p.resolve().relative_to((project / "cram").resolve())
        in_cram = True
    except ValueError:
        pass

    if in_cram and p.name.endswith(".cram.stats"):
        rank = 0
    elif p.name.endswith(".cram.stats"):
        rank = 1
    elif p.name.endswith(".samtools_stats.txt"):
        rank = 2
    elif p.name.endswith(".cram.stats.SN.txt"):
        rank = 5
    elif p.name.endswith(".stats.txt"):
        rank = 3
    else:
        rank = 4
    return rank, str(p)


def discover_cram(project: Path, ident: Dict[str, str]) -> Optional[Path]:
    candidates = []
    cram_root = project / "cram"
    if not cram_root.exists():
        return None
    for p in cram_root.rglob("*.cram"):
        pi = identity_from_name(p)
        if pi and row_key(pi) == row_key(ident):
            candidates.append(p)
    if not candidates:
        return None
    candidates.sort(key=lambda p: (".hg38." not in p.name.lower(), str(p)))
    if len(candidates) > 1:
        warn(f"Multiple CRAMs match {row_key(ident)}; using {rel(candidates[0], project)}")
    return candidates[0]


def docker_base(project: Path, reference: Path, image: str) -> Tuple[List[str], str, str]:
    project = project.resolve()
    reference = reference.resolve()

    try:
        ref_rel = reference.relative_to(project)
        mounts = ["-v", f"{project}:/work:ro"]
        ref_in = f"/work/{ref_rel}"
    except ValueError:
        mounts = [
            "-v", f"{project}:/work:ro",
            "-v", f"{reference.parent}:/reference:ro",
        ]
        ref_in = f"/reference/{reference.name}"

    cmd = [
        "docker", "run", "--rm",
        "-u", f"{os.getuid()}:{os.getgid()}",
        *mounts,
        "-w", "/work",
        image,
    ]
    return cmd, "/work", ref_in


def check_docker() -> None:
    try:
        subprocess.run(
            ["docker", "version", "--format", "{{.Server.Version}}"],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "Docker is required for CRAM-derived full metrics but is not available."
        ) from exc


def docker_samtools_version(project: Path, reference: Path, image: str) -> str:
    base, _, _ = docker_base(project, reference, image)
    proc = subprocess.run(
        base + ["samtools", "--version"],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    first = proc.stdout.splitlines()[0].strip()
    return first.replace("samtools ", "", 1)


def cram_path_in_container(cram: Path, project: Path) -> str:
    return f"/work/{cram.resolve().relative_to(project.resolve())}"


def compute_soft_clipping(
    project: Path, cram: Path, reference: Path, image: str, threads: int,
    total_length
) -> Dict[str, object]:
    base, _, ref_in = docker_base(project, reference, image)
    cram_in = cram_path_in_container(cram, project)

    cmd = base + [
        "samtools", "view",
        "-@", str(threads),
        "-T", ref_in,
        "-F", str(EXCLUDE_FLAGS),
        cram_in,
    ]

    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    assert proc.stdout is not None
    soft = 0

    cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")
    for line in proc.stdout:
        f = line.split("\t", 6)
        if len(f) < 6:
            continue
        cigar = f[5]
        if cigar == "*":
            continue
        soft += sum(int(n) for n, op in cigar_re.findall(cigar) if op == "S")

    stderr = proc.stderr.read() if proc.stderr else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"samtools view failed for {cram.name}: {stderr.strip()}")

    return {
        "soft_clipped_bases": soft,
        "soft_clipped_percent": pct(soft, total_length),
    }


def compute_coverage(
    project: Path, cram: Path, reference: Path, image: str, threads: int
) -> Dict[str, object]:
    """
    Produce a depth histogram inside Docker, not one Python object per genome base.

    samtools view -T makes CRAM decoding independent of host reference lookup.
    samtools depth -aa includes zero-depth positions.
    """
    base, _, ref_in = docker_base(project, reference, image)
    cram_in = cram_path_in_container(cram, project)

    awk = (
        r'{c[$3]++; n++; s+=$3} '
        r'END {for (d in c) print d "\t" c[d]; print "SUM\t" s; print "N\t" n}'
    )
    shell = (
        f"samtools view -@ {threads} -T {shlex.quote(ref_in)} -b "
        f"{shlex.quote(cram_in)} | "
        f"samtools depth -aa -G {EXCLUDE_FLAGS} - | "
        f"awk {shlex.quote(awk)}"
    )

    proc = subprocess.run(
        base + ["sh", "-c", shell],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    hist: Dict[int, int] = {}
    total_sum = total_n = None

    for line in proc.stdout.splitlines():
        a, b = line.split("\t", 1)
        if a == "SUM":
            total_sum = int(float(b))
        elif a == "N":
            total_n = int(float(b))
        else:
            hist[int(a)] = int(float(b))

    if not total_n:
        return {
            "mean_coverage": None, "median_coverage": None,
            "breadth_1x_percent": None, "breadth_10x_percent": None,
            "breadth_20x_percent": None, "breadth_30x_percent": None,
        }

    mean_cov = total_sum / total_n if total_sum is not None else None

    target1 = (total_n - 1) // 2
    target2 = total_n // 2
    cumulative = 0
    med1 = med2 = None
    for depth in sorted(hist):
        prev = cumulative
        cumulative += hist[depth]
        if med1 is None and prev <= target1 < cumulative:
            med1 = depth
        if med2 is None and prev <= target2 < cumulative:
            med2 = depth
            break
    median_cov = (med1 + med2) / 2.0 if med1 is not None and med2 is not None else None

    def breadth(threshold: int) -> float:
        covered = sum(count for depth, count in hist.items() if depth >= threshold)
        return covered / total_n * 100.0

    return {
        "mean_coverage": mean_cov,
        "median_coverage": median_cov,
        "breadth_1x_percent": breadth(1),
        "breadth_10x_percent": breadth(10),
        "breadth_20x_percent": breadth(20),
        "breadth_30x_percent": breadth(30),
    }


def load_run_metadata(path: Path) -> Dict[Tuple[str, str, str, str, str], Dict[str, str]]:
    if not path.exists():
        return {}

    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample", "read_technology", "coverage", "reference", "mapper_tag"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain: " + ", ".join(sorted(required))
            )

        data = {}
        for row in reader:
            tech_raw = row["read_technology"].strip().lower()
            tech = "ONT" if tech_raw == "ont" else "PacBio" if tech_raw in {"pb", "pacbio"} else row["read_technology"].strip()
            tag_raw = row["mapper_tag"].strip().lower()
            canonical = METHOD_SPECS.get(tag_raw, (tag_raw, "", ""))[0]
            key = (
                row["sample"].strip().upper(),
                tech,
                row["coverage"].strip().lower(),
                row["reference"].strip().lower() or "hg38",
                canonical,
            )
            data[key] = row
        return data


def discover_snakemake_benchmark(project: Path, ident: Dict[str, str]) -> Dict[str, str]:
    """
    Best-effort reader for standard Snakemake benchmark TSVs.

    It intentionally does not guess a benchmark if several equally plausible
    files exist.
    """
    roots = [
        project / "benchmarks",
        project / "benchmark",
        project / "alignment_analysis" / "benchmarks",
        project / "alignment_analysis" / "benchmark",
    ]
    tokens = [
        ident["sample"].lower(),
        ident["coverage"].lower(),
        ident["mapper_tag"].lower(),
    ]
    tech_token = "ont" if ident["read_technology"] == "ONT" else "pb"

    candidates = []
    for root in roots:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            n = p.name.lower()
            if all(t in n for t in tokens) and tech_token in n:
                candidates.append(p)

    if len(candidates) != 1:
        return {}

    p = candidates[0]
    try:
        with p.open(encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            row = next(reader, None)
            if not row:
                return {}
    except Exception:
        return {}

    out = {}

    # Standard Snakemake benchmark columns normally include "s" and "max_rss".
    sec = num(row.get("s"))
    if sec is not None:
        out["runtime_seconds"] = fmt(sec)

    rss = num(row.get("max_rss"))
    if rss is not None:
        # Snakemake benchmark max_rss is MB.
        out["peak_ram_mb"] = fmt(rss)

    return out


def complete_status(row: Dict[str, object]) -> Tuple[str, str]:
    missing = [k for k in FULL_REQUIRED if str(row.get(k, NA)).strip().upper() == NA]
    return ("YES" if not missing else "NO", ",".join(missing) if missing else "")


def build(args) -> List[Dict[str, str]]:
    project = args.project.resolve()
    metadata_path = args.run_metadata.resolve()
    run_metadata = load_run_metadata(metadata_path)

    selected = {}
    for p in discover_stats(project):
        ident = identity_from_name(p)
        if not ident:
            continue
        key = row_key(ident)
        if key not in selected or stats_priority(p, project) < stats_priority(selected[key][0], project):
            selected[key] = (p, ident)

    if not selected:
        raise RuntimeError("No recognized samtools stats files found.")

    samtools_version = NA
    if args.reference_fasta:
        check_docker()
        samtools_version = docker_samtools_version(
            project, args.reference_fasta, args.samtools_image
        )

    rows = []
    for key, (stats_file, ident) in selected.items():
        print(f"Processing {ident['dataset']} {ident['mapper_tag']} ...", file=sys.stderr)

        row: Dict[str, object] = {
            **ident,
            "statistics_file": rel(stats_file, project),
            **parse_stats(stats_file),
        }

        cram = discover_cram(project, ident)
        row["cram_file"] = rel(cram, project) if cram else NA

        # CRAM-derived metrics
        for k in [
            "soft_clipped_bases", "soft_clipped_percent",
            "mean_coverage", "median_coverage",
            "breadth_1x_percent", "breadth_10x_percent",
            "breadth_20x_percent", "breadth_30x_percent",
        ]:
            row[k] = NA

        if args.reference_fasta and cram:
            try:
                row.update(compute_soft_clipping(
                    project, cram, args.reference_fasta, args.samtools_image,
                    args.threads, row.get("total_length")
                ))
                row.update(compute_coverage(
                    project, cram, args.reference_fasta, args.samtools_image,
                    args.threads
                ))
            except subprocess.CalledProcessError as exc:
                warn(f"Coverage/clipping failed for {cram.name}: {exc.stderr.strip() if exc.stderr else exc}")
            except RuntimeError as exc:
                warn(str(exc))

        # Run provenance/performance.
        for k in RUN_METADATA_COLUMNS:
            row[k] = NA

        auto_bench = discover_snakemake_benchmark(project, ident)
        for k, v in auto_bench.items():
            row[k] = v

        explicit = run_metadata.get(key, {})
        for k in RUN_METADATA_COLUMNS:
            v = explicit.get(k, "")
            if str(v).strip():
                row[k] = str(v).strip()

        row["samtools_version"] = samtools_version
        row["samtools_docker_image"] = (
            args.samtools_image if args.reference_fasta else NA
        )

        # Format all numeric objects.
        clean = {}
        for col in OUTPUT_COLUMNS:
            v = row.get(col, NA)
            clean[col] = fmt(v) if isinstance(v, (int, float)) or v is None else (str(v) if str(v).strip() else NA)

        complete, missing = complete_status(clean)
        clean["metrics_complete"] = complete
        clean["missing_metrics"] = missing

        rows.append(clean)

    order_aligner = {"minimap2": 0, "pbmm2": 1, "VACMap": 2, "VG Giraffe": 3}
    order_tech = {"ONT": 0, "PacBio": 1}
    rows.sort(key=lambda r: (
        r["sample"],
        order_tech.get(r["read_technology"], 9),
        r["coverage"],
        order_aligner.get(r["aligner"], 9),
    ))
    return rows


def parse_args():
    here = Path(__file__).resolve()
    default_project = here.parents[3] if len(here.parents) >= 4 else Path.cwd()

    p = argparse.ArgumentParser(
        description="Build full long-read alignment benchmark summary."
    )
    p.add_argument("--project", type=Path, default=default_project)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument(
        "--reference-fasta", type=Path, default=None,
        help="hg38 FASTA. Required to automatically compute CRAM coverage/clipping."
    )
    p.add_argument(
        "--run-metadata", type=Path, default=None,
        help="Optional alignment_run_metadata.tsv with mapper runtime/provenance."
    )
    p.add_argument("--threads", type=int, default=8)
    p.add_argument(
        "--samtools-image",
        default=SAMTOOLS_IMAGE_DEFAULT,
        help="Pinned Docker image used only for CRAM-derived metrics."
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.project = args.project.expanduser().resolve()
    args.out = (
        args.out.expanduser().resolve()
        if args.out
        else args.project / "alignment_analysis" / "tables" / "alignment_summary.tsv"
    )
    args.run_metadata = (
        args.run_metadata.expanduser().resolve()
        if args.run_metadata
        else args.project / "alignment_analysis" / "tables" / "alignment_run_metadata.tsv"
    )
    if args.reference_fasta:
        args.reference_fasta = args.reference_fasta.expanduser().resolve()
        if not args.reference_fasta.exists():
            print(f"ERROR: reference FASTA not found: {args.reference_fasta}", file=sys.stderr)
            return 1

    try:
        rows = build(args)
    except (RuntimeError, ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    complete = sum(r["metrics_complete"] == "YES" for r in rows)
    print(f"Wrote {len(rows)} rows to {args.out}")
    print(f"Rows with every requested metric present: {complete}/{len(rows)}")
    if complete != len(rows):
        print(
            "Any remaining NA values are listed in the missing_metrics column. "
            "Mapper runtime/RAM/version/image cannot be reconstructed if the "
            "original mapping workflow never recorded them."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
