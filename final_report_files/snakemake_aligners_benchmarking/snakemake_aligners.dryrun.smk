#####################################################################################
# ALIGNMENT PIPELINE: MINIMAP2, PBMM2, VACMAP, AND VG GIRAFFE
#####################################################################################

from pathlib import Path

configfile: "config.dryrun.yaml"


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
