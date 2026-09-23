#!/usr/bin/env bash
set -euo pipefail

echo "============================================================"
echo " FULL ALIGNER PIPELINE SETUP: minimap2 + pbmm2 + VACmap + vg"
echo "============================================================"

WORKDIR="$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"
PROJECT_ROOT="$HOME/lrs_benchmarking"

cd "$WORKDIR"

echo "[1/10] Creating folders..."
mkdir -p envs scripts data/samples data/reference data/graph results/logs results/bam results/stats results/tables results/figures tools

echo "[2/10] Backing up old files..."
STAMP="$(date +%Y%m%d_%H%M%S)"
for f in snakemake_aligners.smk config.yaml; do
    if [ -f "$f" ]; then
        cp "$f" "${f}.backup_${STAMP}"
    fi
done

echo "[3/10] Linking FASTQ files from main project..."
MISSING_INPUTS=0

for sample in HG002 HG003 HG004; do
    for tech in ont pb; do
        SRC="$PROJECT_ROOT/samples_try/${sample}.${tech}.1k.fastq.gz"
        DEST="data/samples/${sample}.${tech}.1k.fastq.gz"

        if [ -f "$SRC" ]; then
            ln -sf "$SRC" "$DEST"
            echo "Linked: $DEST"
        else
            echo "MISSING FASTQ: $SRC"
            MISSING_INPUTS=1
        fi
    done
done

echo "[4/10] Linking reference genome..."
KNOWN_REF="$PROJECT_ROOT/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"

if [ -f "$KNOWN_REF" ]; then
    ln -sf "$KNOWN_REF" data/reference/genome.fasta
    echo "Linked reference: data/reference/genome.fasta"
else
    FIRST_REF="$(find "$PROJECT_ROOT/reference" -maxdepth 1 -type f -name "*.fasta" | head -n 1 || true)"
    if [ -n "$FIRST_REF" ]; then
        ln -sf "$FIRST_REF" data/reference/genome.fasta
        echo "Linked fallback reference: data/reference/genome.fasta"
    else
        echo "MISSING REFERENCE FASTA in $PROJECT_ROOT/reference"
        MISSING_INPUTS=1
    fi
fi

if [ "$MISSING_INPUTS" -eq 1 ]; then
    echo
    echo "STOP: Some FASTQ/reference input files are missing."
    echo "Fix those paths first, then run this script again."
    exit 1
fi

echo "[5/10] Writing config.yaml with all four aligners active..."
cat > config.yaml <<'YAML'
samples_dir: data/samples
reference: data/reference/genome.fasta

samples:
  - HG002
  - HG003
  - HG004

technologies:
  - ont
  - pb

active_aligners:
  - minimap2
  - pbmm2
  - vacmap
  - vg_giraffe

vacmap_bin: "conda run -n vacmap_env vacmap"

vg_gbz: data/graph/graph.gbz
vg_min: data/graph/graph.min
vg_dist: data/graph/graph.dist
vg_zipcodes: ""
YAML

echo "[6/10] Writing Conda/Mamba environment files..."
cat > envs/alignment.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - python=3.11
  - minimap2
  - pbmm2
  - samtools
  - pandas
  - numpy
  - scipy
  - matplotlib
  - pip
YAML

cat > envs/vg.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - vg
  - samtools
YAML

echo "[7/10] Writing Snakefile..."
cat > snakemake_aligners.smk <<'SNAKE'
#####################################################################################
# ALIGNMENT PIPELINE: MINIMAP2, PBMM2, VACMAP, AND VG GIRAFFE
#####################################################################################

from pathlib import Path

configfile: "config.yaml"


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

SAMPLES_DIR = config.get("samples_dir", "data/samples")
REFERENCE = config.get("reference", "data/reference/genome.fasta")

SAMPLES = config.get("samples", ["HG002", "HG003", "HG004"])
TECHNOLOGIES = config.get("technologies", ["ont", "pb"])
ACTIVE_ALIGNERS = config.get("active_aligners", ["minimap2", "pbmm2", "vacmap", "vg_giraffe"])

VACMAP_BIN = config.get("vacmap_bin", "conda run -n vacmap_env vacmap")

VG_GBZ = config.get("vg_gbz", "data/graph/graph.gbz")
VG_MIN = config.get("vg_min", "data/graph/graph.min")
VG_DIST = config.get("vg_dist", "data/graph/graph.dist")
VG_ZIPCODES = config.get("vg_zipcodes", "")

PREFLIGHT_OK = "results/logs/preflight.ok"

ENV_ALIGNMENT = "envs/alignment.yaml"
ENV_VG = "envs/vg.yaml"


wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES),
    aligner="|".join(ACTIVE_ALIGNERS)


# -----------------------------------------------------------------------------
# Preset / mode functions
# -----------------------------------------------------------------------------

def minimap2_preset(wildcards):
    if wildcards.technology == "ont":
        return "map-ont"
    if wildcards.technology == "pb":
        return "map-hifi"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def pbmm2_preset(wildcards):
    if wildcards.technology == "ont":
        return "SUBREAD"
    if wildcards.technology == "pb":
        return "HIFI"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vacmap_mode(wildcards):
    if wildcards.technology == "ont":
        return "H"
    if wildcards.technology == "pb":
        return "L"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vg_giraffe_preset(wildcards):
    if wildcards.technology == "ont":
        return "r10"
    if wildcards.technology == "pb":
        return "hifi"
    raise ValueError(f"Unknown technology: {wildcards.technology}")


def vg_zipcodes_arg(wildcards):
    if VG_ZIPCODES:
        return f"-z {VG_ZIPCODES}"
    return ""


# -----------------------------------------------------------------------------
# Final target
# -----------------------------------------------------------------------------

rule all:
    input:
        PREFLIGHT_OK,
        "results/tables/alignment_summary.tsv",
        "results/tables/quality_report.tsv",
        "results/tables/plot_data_clean.tsv",
        "results/tables/plot_data_pivot.tsv",
        "results/tables/plot_arrays.npz",
        "results/tables/paired_alignment_tests.tsv",
        "results/figures/final_alignment_error_rate.png",
        "results/figures/final_alignment_error_rate.pdf",
        "results/figures/mapped_bases.png",
        "results/figures/mapped_bases.pdf",
        "results/figures/mapped_bases_cigar.png",
        "results/figures/mapped_bases_cigar.pdf"


# -----------------------------------------------------------------------------
# Preflight check
# -----------------------------------------------------------------------------

rule preflight:
    output:
        PREFLIGHT_OK
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)

        errors = []

        for required_file in [ENV_ALIGNMENT, ENV_VG, REFERENCE]:
            if not Path(required_file).exists():
                errors.append(f"Missing required file: {required_file}")

        for sample in SAMPLES:
            for technology in TECHNOLOGIES:
                fastq = Path(f"{SAMPLES_DIR}/{sample}.{technology}.1k.fastq.gz")
                if not fastq.exists():
                    errors.append(f"Missing FASTQ file: {fastq}")

        if "vg_giraffe" in ACTIVE_ALIGNERS:
            for graph_file in [VG_GBZ, VG_MIN, VG_DIST]:
                if not Path(graph_file).exists():
                    errors.append(f"Missing vg Giraffe graph/index file: {graph_file}")

            if VG_ZIPCODES and not Path(VG_ZIPCODES).exists():
                errors.append(f"Missing vg Giraffe zipcodes file: {VG_ZIPCODES}")

        required_scripts = [
            "scripts/alignment_summary.py",
            "scripts/quality_check.py",
            "scripts/prepare_plot_data.py",
            "scripts/paired_t_tests.py",
            "scripts/alignment_metrics.py",
        ]

        for script in required_scripts:
            if not Path(script).exists():
                errors.append(f"Missing Python script: {script}")

        if errors:
            message = "\n".join(errors)
            raise RuntimeError(
                "\nPREFLIGHT CHECK FAILED\n"
                "Fix these problems before running the workflow:\n\n"
                f"{message}\n"
            )

        Path(output[0]).write_text("Preflight check passed.\n")


# -----------------------------------------------------------------------------
# Rule 1: minimap2 alignment
# -----------------------------------------------------------------------------

rule align_minimap2:
    input:
        preflight=PREFLIGHT_OK,
        fastq=f"{SAMPLES_DIR}/{{sample}}.{{technology}}.1k.fastq.gz",
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.minimap2.sorted.bam",
        bai="results/bam/{sample}.{technology}.minimap2.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        preset=minimap2_preset
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fastq} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 2: pbmm2 alignment
# -----------------------------------------------------------------------------

rule align_pbmm2:
    input:
        preflight=PREFLIGHT_OK,
        fastq=f"{SAMPLES_DIR}/{{sample}}.{{technology}}.1k.fastq.gz",
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.pbmm2.sorted.bam",
        bai="results/bam/{sample}.{technology}.pbmm2.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        preset=pbmm2_preset
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        pbmm2 align \
            --preset {params.preset} \
            --sort \
            -j {threads} \
            {input.ref} \
            {input.fastq} \
            {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 3: VACmap alignment
# -----------------------------------------------------------------------------

rule align_vacmap:
    input:
        preflight=PREFLIGHT_OK,
        fastq=f"{SAMPLES_DIR}/{{sample}}.{{technology}}.1k.fastq.gz",
        ref=REFERENCE
    output:
        bam="results/bam/{sample}.{technology}.vacmap.sorted.bam",
        bai="results/bam/{sample}.{technology}.vacmap.sorted.bam.bai"
    conda:
        ENV_ALIGNMENT
    params:
        mode=vacmap_mode,
        vacmap_bin=VACMAP_BIN
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        {params.vacmap_bin} \
            -ref {input.ref} \
            -read {input.fastq} \
            -mode {params.mode} \
            -t {threads} \
            --force \
            -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 4: vg Giraffe alignment
# -----------------------------------------------------------------------------

rule align_vg_giraffe:
    input:
        preflight=PREFLIGHT_OK,
        fastq=f"{SAMPLES_DIR}/{{sample}}.{{technology}}.1k.fastq.gz",
        gbz=VG_GBZ,
        minimizer=VG_MIN,
        dist=VG_DIST
    output:
        bam="results/bam/{sample}.{technology}.vg_giraffe.sorted.bam",
        bai="results/bam/{sample}.{technology}.vg_giraffe.sorted.bam.bai"
    conda:
        ENV_VG
    params:
        preset=vg_giraffe_preset,
        zipcodes=vg_zipcodes_arg
    threads: 4
    shell:
        r"""
        mkdir -p results/bam

        vg giraffe \
            -Z {input.gbz} \
            -m {input.minimizer} \
            -d {input.dist} \
            {params.zipcodes} \
            -f {input.fastq} \
            -t {threads} \
            -b {params.preset} \
            -o BAM \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """


# -----------------------------------------------------------------------------
# Rule 5: Samtools quality statistics
# -----------------------------------------------------------------------------

rule samtools_quality:
    input:
        preflight=PREFLIGHT_OK,
        bam="results/bam/{sample}.{technology}.{aligner}.sorted.bam",
        bai="results/bam/{sample}.{technology}.{aligner}.sorted.bam.bai"
    output:
        flagstat="results/stats/{sample}.{technology}.{aligner}.flagstat.txt",
        idxstat="results/stats/{sample}.{technology}.{aligner}.idxstat.txt",
        stats="results/stats/{sample}.{technology}.{aligner}.stats.txt"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/stats

        samtools flagstat {input.bam} > {output.flagstat}
        samtools idxstats {input.bam} > {output.idxstat}
        samtools stats {input.bam} > {output.stats}
        """


# -----------------------------------------------------------------------------
# Rule 6: Build alignment summary table
# -----------------------------------------------------------------------------

rule alignment_summary:
    input:
        stats=expand(
            "results/stats/{sample}.{technology}.{aligner}.stats.txt",
            sample=SAMPLES,
            technology=TECHNOLOGIES,
            aligner=ACTIVE_ALIGNERS,
        )
    output:
        summary="results/tables/alignment_summary.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/tables

        python scripts/alignment_summary.py \
            --stats {input.stats} \
            --out {output.summary}
        """


# -----------------------------------------------------------------------------
# Rule 7: Quality check of the alignment summary
# -----------------------------------------------------------------------------

rule quality_check:
    input:
        summary="results/tables/alignment_summary.tsv"
    output:
        report="results/tables/quality_report.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/quality_check.py \
            --input {input.summary} \
            --out {output.report}
        """


# -----------------------------------------------------------------------------
# Rule 8: Prepare data for plotting
# -----------------------------------------------------------------------------

rule prepare_plot_data:
    input:
        summary="results/tables/alignment_summary.tsv",
        quality_report="results/tables/quality_report.tsv"
    output:
        clean="results/tables/plot_data_clean.tsv",
        pivot="results/tables/plot_data_pivot.tsv",
        arrays="results/tables/plot_arrays.npz"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/prepare_plot_data.py \
            --input {input.summary} \
            --quality-report {input.quality_report} \
            --cleaned-data {output.clean} \
            --pivoted-data {output.pivot} \
            --arrays {output.arrays}
        """


# -----------------------------------------------------------------------------
# Rule 9: Paired t-tests
# -----------------------------------------------------------------------------

rule paired_t_tests:
    input:
        pivot="results/tables/plot_data_pivot.tsv"
    output:
        tests="results/tables/paired_alignment_tests.tsv"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        python scripts/paired_t_tests.py \
            --input {input.pivot} \
            --out {output.tests}
        """


# -----------------------------------------------------------------------------
# Rule 10: Final alignment figures
# -----------------------------------------------------------------------------

rule alignment_figures:
    input:
        clean="results/tables/plot_data_clean.tsv",
        pivot="results/tables/plot_data_pivot.tsv",
        arrays="results/tables/plot_arrays.npz",
        tests="results/tables/paired_alignment_tests.tsv"
    output:
        error_png="results/figures/final_alignment_error_rate.png",
        error_pdf="results/figures/final_alignment_error_rate.pdf",
        mapped_png="results/figures/mapped_bases.png",
        mapped_pdf="results/figures/mapped_bases.pdf",
        cigar_png="results/figures/mapped_bases_cigar.png",
        cigar_pdf="results/figures/mapped_bases_cigar.pdf"
    conda:
        ENV_ALIGNMENT
    shell:
        r"""
        mkdir -p results/figures

        python scripts/alignment_metrics.py \
            --cleaned-data {input.clean} \
            --pivoted-data {input.pivot} \
            --arrays {input.arrays} \
            --tests {input.tests} \
            --error-png {output.error_png} \
            --error-pdf {output.error_pdf} \
            --mapped-png {output.mapped_png} \
            --mapped-pdf {output.mapped_pdf} \
            --cigar-png {output.cigar_png} \
            --cigar-pdf {output.cigar_pdf}
        """
SNAKE

echo "[8/10] Writing Python scripts..."

cat > scripts/alignment_summary.py <<'PY'
#!/usr/bin/env python3

import argparse
from pathlib import Path
import math
import pandas as pd


def parse_number(value):
    value = str(value).strip().replace(",", "")
    try:
        if "." in value or "e" in value.lower():
            return float(value)
        return int(value)
    except ValueError:
        return math.nan


def parse_samtools_stats(path):
    data = {}

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            key = parts[1].strip().rstrip(":").lower()
            value = parse_number(parts[2])
            data[key] = value

    name = Path(path).name.replace(".stats.txt", "")
    parts = name.split(".", 2)

    if len(parts) == 3:
        sample, technology, aligner = parts
    else:
        sample, technology, aligner = "unknown", "unknown", name

    raw_total_sequences = data.get("raw total sequences", math.nan)
    reads_mapped = data.get("reads mapped", math.nan)
    reads_unmapped = data.get("reads unmapped", math.nan)
    total_length = data.get("total length", math.nan)
    bases_mapped = data.get("bases mapped", math.nan)
    bases_mapped_cigar = data.get("bases mapped (cigar)", math.nan)
    mismatches = data.get("mismatches", math.nan)
    error_rate = data.get("error rate", math.nan)
    insertions = data.get("insertions", math.nan)
    deletions = data.get("deletions", math.nan)
    average_length = data.get("average length", math.nan)
    maximum_length = data.get("maximum length", math.nan)
    non_primary_alignments = data.get("non-primary alignments", math.nan)

    mapped_reads_percent = (
        reads_mapped / raw_total_sequences * 100
        if raw_total_sequences and not math.isnan(raw_total_sequences)
        else math.nan
    )

    mapped_bases_percent = (
        bases_mapped_cigar / total_length * 100
        if total_length and not math.isnan(total_length)
        else math.nan
    )

    error_percent = error_rate * 100 if not math.isnan(error_rate) else math.nan

    return {
        "sample": sample,
        "read_technology": technology,
        "aligner": aligner,
        "statistics_file": str(path),
        "raw_total_sequences": raw_total_sequences,
        "reads_mapped": reads_mapped,
        "reads_unmapped": reads_unmapped,
        "mapped_reads_percent": mapped_reads_percent,
        "total_length": total_length,
        "bases_mapped": bases_mapped,
        "bases_mapped_cigar": bases_mapped_cigar,
        "mapped_bases_percent": mapped_bases_percent,
        "mismatches": mismatches,
        "error_rate": error_rate,
        "error_percent": error_percent,
        "insertions": insertions,
        "deletions": deletions,
        "average_length": average_length,
        "maximum_length": maximum_length,
        "non_primary_alignments": non_primary_alignments,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = [parse_samtools_stats(path) for path in args.stats]
    table = pd.DataFrame(rows)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
PY

cat > scripts/quality_check.py <<'PY'
#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def add(rows, check, status, detail):
    rows.append({"check": check, "status": status, "detail": detail})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    rows = []

    add(rows, "rows_exist", "PASS" if len(df) > 0 else "FAIL", f"rows={len(df)}")

    required = ["sample", "read_technology", "aligner", "statistics_file"]
    missing_required = [column for column in required if column not in df.columns]
    add(
        rows,
        "required_columns",
        "PASS" if not missing_required else "FAIL",
        "missing=" + ",".join(missing_required) if missing_required else "all required columns present",
    )

    if not missing_required:
        duplicates = df.duplicated(subset=["sample", "read_technology", "aligner"]).sum()
        add(rows, "duplicate_conditions", "PASS" if duplicates == 0 else "FAIL", f"duplicates={duplicates}")

    missing_values = int(df.isna().sum().sum())
    add(rows, "missing_values", "PASS" if missing_values == 0 else "WARN", f"missing_values={missing_values}")

    numeric_columns = [
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "total_length",
        "bases_mapped",
        "bases_mapped_cigar",
        "error_rate",
        "error_percent",
    ]

    existing_numeric = [column for column in numeric_columns if column in df.columns]
    negative_count = 0

    for column in existing_numeric:
        values = pd.to_numeric(df[column], errors="coerce")
        negative_count += int((values < 0).sum())

    add(rows, "negative_numeric_values", "PASS" if negative_count == 0 else "FAIL", f"negative_values={negative_count}")

    report = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
PY

cat > scripts/prepare_plot_data.py <<'PY'
#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--quality-report", required=True)
    parser.add_argument("--cleaned-data", required=True)
    parser.add_argument("--pivoted-data", required=True)
    parser.add_argument("--arrays", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    numeric_columns = [
        "raw_total_sequences",
        "reads_mapped",
        "reads_unmapped",
        "mapped_reads_percent",
        "total_length",
        "bases_mapped",
        "bases_mapped_cigar",
        "mapped_bases_percent",
        "mismatches",
        "error_rate",
        "error_percent",
        "insertions",
        "deletions",
        "average_length",
        "maximum_length",
        "non_primary_alignments",
    ]

    for column in numeric_columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    Path(args.cleaned_data).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.cleaned_data, sep="\t", index=False)

    pivot = df.copy()
    pivot.to_csv(args.pivoted_data, sep="\t", index=False)

    np.savez(
        args.arrays,
        error_percent=df["error_percent"].to_numpy() if "error_percent" in df.columns else np.array([]),
        bases_mapped=df["bases_mapped"].to_numpy() if "bases_mapped" in df.columns else np.array([]),
        bases_mapped_cigar=df["bases_mapped_cigar"].to_numpy() if "bases_mapped_cigar" in df.columns else np.array([]),
    )


if __name__ == "__main__":
    main()
PY

cat > scripts/paired_t_tests.py <<'PY'
#!/usr/bin/env python3

import argparse
from pathlib import Path
import itertools
import pandas as pd
from scipy.stats import ttest_rel


def safe_ttest(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan")
    result = ttest_rel(a, b, nan_policy="omit")
    return result.statistic, result.pvalue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    metrics = [
        "error_percent",
        "mapped_reads_percent",
        "mapped_bases_percent",
        "bases_mapped",
        "bases_mapped_cigar",
    ]

    rows = []

    for technology in sorted(df["read_technology"].dropna().unique()):
        tech_df = df[df["read_technology"] == technology]

        aligners = sorted(tech_df["aligner"].dropna().unique())

        for metric in metrics:
            if metric not in tech_df.columns:
                continue

            metric_df = tech_df[["sample", "aligner", metric]].copy()
            metric_df[metric] = pd.to_numeric(metric_df[metric], errors="coerce")

            wide = metric_df.pivot_table(index="sample", columns="aligner", values=metric, aggfunc="mean")

            for aligner_a, aligner_b in itertools.combinations(aligners, 2):
                if aligner_a not in wide.columns or aligner_b not in wide.columns:
                    continue

                paired = wide[[aligner_a, aligner_b]].dropna()
                stat, pvalue = safe_ttest(paired[aligner_a], paired[aligner_b])

                rows.append(
                    {
                        "read_technology": technology,
                        "metric": metric,
                        "aligner_a": aligner_a,
                        "aligner_b": aligner_b,
                        "n_pairs": len(paired),
                        "t_statistic": stat,
                        "p_value": pvalue,
                    }
                )

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
PY

cat > scripts/alignment_metrics.py <<'PY'
#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


ALIGNER_ORDER = ["minimap2", "pbmm2", "vacmap", "vg_giraffe"]
TECH_ORDER = ["ont", "pb"]


def clean_label(value):
    mapping = {
        "ont": "ONT",
        "pb": "PacBio HiFi",
        "minimap2": "minimap2",
        "pbmm2": "pbmm2",
        "vacmap": "VACmap",
        "vg_giraffe": "vg Giraffe",
    }
    return mapping.get(value, value)


def plot_metric(df, metric, ylabel, title, png, pdf):
    df = df.copy()

    if metric not in df.columns:
        raise ValueError(f"Missing metric column: {metric}")

    df[metric] = pd.to_numeric(df[metric], errors="coerce")
    df = df.dropna(subset=[metric])

    if df.empty:
        raise ValueError(f"No valid numeric data for metric: {metric}")

    aligners = [a for a in ALIGNER_ORDER if a in set(df["aligner"])]
    technologies = [t for t in TECH_ORDER if t in set(df["read_technology"])]

    fig, ax = plt.subplots(figsize=(10, 6))

    x_positions = list(range(len(aligners)))
    offsets = {}

    if len(technologies) == 1:
        offsets[technologies[0]] = 0
    else:
        offsets = {technologies[0]: -0.15, technologies[1]: 0.15}

    for tech in technologies:
        tech_df = df[df["read_technology"] == tech]

        means = []
        xs = []

        for i, aligner in enumerate(aligners):
            values = tech_df[tech_df["aligner"] == aligner][metric]
            means.append(values.mean())
            xs.append(i + offsets.get(tech, 0))

            sample_values = tech_df[tech_df["aligner"] == aligner]
            ax.scatter(
                [i + offsets.get(tech, 0)] * len(sample_values),
                sample_values[metric],
                alpha=0.75,
                s=45,
            )

        ax.plot(xs, means, marker="o", linewidth=2, label=clean_label(tech))

    ax.set_xticks(x_positions)
    ax.set_xticklabels([clean_label(a) for a in aligners], rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Alignment configuration")
    ax.set_title(title)
    ax.legend(title="Read technology")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    Path(png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cleaned-data", required=True)
    parser.add_argument("--pivoted-data", required=True)
    parser.add_argument("--arrays", required=True)
    parser.add_argument("--tests", required=True)
    parser.add_argument("--error-png", required=True)
    parser.add_argument("--error-pdf", required=True)
    parser.add_argument("--mapped-png", required=True)
    parser.add_argument("--mapped-pdf", required=True)
    parser.add_argument("--cigar-png", required=True)
    parser.add_argument("--cigar-pdf", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.cleaned_data, sep="\t")

    plot_metric(
        df=df,
        metric="error_percent",
        ylabel="Alignment error rate (%)",
        title="Alignment error rate across aligners",
        png=args.error_png,
        pdf=args.error_pdf,
    )

    plot_metric(
        df=df,
        metric="bases_mapped",
        ylabel="Mapped bases",
        title="Mapped bases across aligners",
        png=args.mapped_png,
        pdf=args.mapped_pdf,
    )

    plot_metric(
        df=df,
        metric="bases_mapped_cigar",
        ylabel="Mapped bases from CIGAR",
        title="CIGAR mapped bases across aligners",
        png=args.cigar_png,
        pdf=args.cigar_pdf,
    )


if __name__ == "__main__":
    main()
PY

chmod +x scripts/*.py

echo "[9/10] Installing or updating VACmap..."
if [ -d "tools/VACmap/.git" ]; then
    echo "VACmap folder already exists. Updating..."
    git -C tools/VACmap pull || true
else
    git clone https://github.com/micahvista/VACmap.git tools/VACmap
fi

if conda env list | awk '{print $1}' | grep -qx "vacmap_env"; then
    echo "Conda environment vacmap_env already exists."
else
    if command -v mamba >/dev/null 2>&1; then
        mamba env create --name vacmap_env --file tools/VACmap/VACmap_environment.yml
    else
        conda env create --name vacmap_env --file tools/VACmap/VACmap_environment.yml
    fi
fi

echo "Installing VACmap Python package into vacmap_env..."
(
    cd tools/VACmap
    conda run -n vacmap_env python setup.py install
)

echo "Testing VACmap command..."
conda run -n vacmap_env vacmap -h >/dev/null 2>&1 || true

echo "[10/10] Checking vg Giraffe graph files..."
GRAPH_OK=1

for graph_file in data/graph/graph.gbz data/graph/graph.min data/graph/graph.dist; do
    if [ ! -f "$graph_file" ]; then
        echo "MISSING VG GRAPH FILE: $graph_file"
        GRAPH_OK=0
    else
        echo "Found: $graph_file"
    fi
done

echo
echo "============================================================"
echo " SETUP FINISHED"
echo "============================================================"

if [ "$GRAPH_OK" -eq 0 ]; then
    echo
    echo "STOP: Everything was created, but vg Giraffe cannot run yet."
    echo
    echo "You need these three files:"
    echo "  data/graph/graph.gbz"
    echo "  data/graph/graph.min"
    echo "  data/graph/graph.dist"
    echo
    echo "After placing them there, run:"
    echo
    echo "  cd \"$WORKDIR\""
    echo "  snakemake -s snakemake_aligners.smk --sdm conda --cores 4 -np"
    echo "  snakemake -s snakemake_aligners.smk --sdm conda --cores 4 -p"
    echo
    exit 0
fi

echo "Running dry-run..."
snakemake -s snakemake_aligners.smk --sdm conda --cores 4 -np

echo
echo "Dry-run passed. Running real workflow now..."
snakemake -s snakemake_aligners.smk --sdm conda --cores 4 -p

echo
echo "============================================================"
echo " PIPELINE COMPLETE"
echo "============================================================"
echo "Main outputs:"
echo "  results/tables/alignment_summary.tsv"
echo "  results/tables/quality_report.tsv"
echo "  results/tables/paired_alignment_tests.tsv"
echo "  results/figures/final_alignment_error_rate.png"
echo "  results/figures/mapped_bases.png"
echo "  results/figures/mapped_bases_cigar.png"
