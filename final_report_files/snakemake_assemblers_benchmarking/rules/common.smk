#####################################################################################
# SHARED ASSEMBLER WORKFLOW CONFIGURATION
#
# Responsibilities:
#   - Read sample metadata from samples.tsv.
#   - Keep machine-specific paths OUT of samples.tsv.
#   - Resolve source FASTQs from input_root automatically.
#   - Support local pre-extracted Chr21 testing and server WGS mode.
#   - Define shared assembler parameters and expected workflow outputs.
#####################################################################################

from pathlib import Path
import csv
import hashlib
import re
import shlex
import shutil
import subprocess


# -----------------------------------------------------------------------------
# General workflow configuration
# -----------------------------------------------------------------------------

PROJECT_DIR = Path(workflow.basedir).resolve()

SAMPLE_SHEET = config.get("sample_sheet", "samples.tsv")
ACTIVE_ASSEMBLERS = config.get("active_assemblers", ["flye"])
GENOME_SIZE = config.get("genome_size", 46709983)
DEFAULT_THREADS = int(config.get("threads", 4))

CHROMOSOME = config.get("chromosome", "chr21")
CHROMOSOME_LENGTH = int(
    config.get("chromosome_length", GENOME_SIZE)
)
MIN_MAPQ = int(
    config.get("min_mapq", 20)
)
TARGET_COVERAGE = float(
    config.get("target_coverage", 30)
)
TARGET_BASES = int(
    config.get(
        "target_bases",
        round(CHROMOSOME_LENGTH * TARGET_COVERAGE)
    )
)

EXPECTED_TARGET_BASES = round(
    CHROMOSOME_LENGTH * TARGET_COVERAGE
)

if TARGET_BASES != EXPECTED_TARGET_BASES:
    raise ValueError(
        "target_bases is inconsistent with chromosome length "
        "and target coverage: "
        f"configured={TARGET_BASES}, "
        f"expected={EXPECTED_TARGET_BASES}"
    )

PROGRESS_DIR = config.get(
    "progress_dir",
    "results/progress"
)

LOG_DIR = config.get(
    "log_dir",
    "results/logs"
)

RESOURCE_CONFIG = (
    config.get("resources", {})
    or {}
)

MAPPING_THREADS = int(
    RESOURCE_CONFIG.get(
        "mapping_threads",
        DEFAULT_THREADS
    )
)

MAPPING_MEM_MB = int(
    RESOURCE_CONFIG.get(
        "mapping_mem_mb",
        8000
    )
)

# Keep the total CPU usage of the mapping pipeline near the Snakemake
# threads allocation. minimap2 and both samtools sort processes can overlap.
MAP_THREADS = max(
    1,
    MAPPING_THREADS - 4
)

SORT_THREADS = max(
    1,
    min(2, MAPPING_THREADS // 4)
)

RAW_REFERENCE = config.get("reference")

REFERENCE = None

if RAW_REFERENCE:
    REFERENCE = Path(
        RAW_REFERENCE
    ).expanduser()

    if not REFERENCE.is_absolute():
        REFERENCE = (
            PROJECT_DIR
            / REFERENCE
        )

    REFERENCE = REFERENCE.resolve(
        strict=False
    )

INPUT_MODE = config.get("chr21_fastq")

if INPUT_MODE != "chr21_fastq":
    raise ValueError(
        "Missing or invalid input_mode. "
        "Expected input_mode to be 'chr21_fastq'."
    )


# -----------------------------------------------------------------------------
# Input root
# -----------------------------------------------------------------------------

RAW_INPUT_ROOT = config.get("/home/nicolas/lrs_benchmarking/final_report_files/snakemake_assemblers_benchmarking/benchmark_chr21_real/inputs/final_normalized")

if not RAW_INPUT_ROOT:
    raise ValueError(
        "Missing input_root. "
        "Set input_root in config/local.yaml or config/server.yaml."
    )

INPUT_ROOT = Path(RAW_INPUT_ROOT).expanduser()

if not INPUT_ROOT.is_absolute():
    INPUT_ROOT = PROJECT_DIR / INPUT_ROOT

INPUT_ROOT = INPUT_ROOT.resolve(strict=False)


# -----------------------------------------------------------------------------
# Input/output naming
# -----------------------------------------------------------------------------

CHR21_OUTPUT_DIR = config.get(
    "chr21_output_dir",
    "results/chr21"
)

# Separate validation markers prevent a successful local validation from being
# accidentally reused during a server run, or vice versa.
VALIDATION_OK = (
    f"results/logs/validate_inputs.{INPUT_MODE}.ok"
)


# -----------------------------------------------------------------------------
# Assembler Docker containers
# -----------------------------------------------------------------------------

CONTAINERS = config.get("containers", {}) or {}
ASSEMBLER_CONTAINER_KEYS = {
    "flye": "flye2",
    "goldrush": "goldrush",
    "ntlink": "ntlink",
    "verkko": "verkko2",
}

missing_container_keys = [
    key for assembler, key in ASSEMBLER_CONTAINER_KEYS.items()
    if assembler in ACTIVE_ASSEMBLERS and key not in CONTAINERS
]

if missing_container_keys:
    raise ValueError(
        "Missing container definitions in config/base.yaml: "
        + ", ".join(sorted(missing_container_keys))
    )

DOCKER_IMAGES = {}

for key, image in CONTAINERS.items():
    image = str(image).strip()

    if image.startswith("docker://"):
        image = image[len("docker://"):]

    if not image:
        raise ValueError(f"Empty Docker image definition: {key}")

    DOCKER_IMAGES[key] = image


def is_within(path, parent):
    """Return True when path is inside parent after path normalization."""
    path = Path(path).resolve(strict=False)
    parent = Path(parent).resolve(strict=False)

    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def docker_mount_args():
    """Mount the project writable and external input data read-only."""
    mounts = ["-v " + shlex.quote(f"{PROJECT_DIR}:{PROJECT_DIR}")]

    if not is_within(INPUT_ROOT, PROJECT_DIR):
        mounts.append("-v " + shlex.quote(f"{INPUT_ROOT}:{INPUT_ROOT}:ro"))

    return " ".join(mounts)


PROJECT_DIR_STR = str(PROJECT_DIR)
DOCKER_MOUNTS = docker_mount_args()

# Assessment remains Conda-based for now; assembler execution does not.
ENV_ASSESSMENT = str(PROJECT_DIR / "envs" / "assessment.yaml")


# -----------------------------------------------------------------------------
# Supported assemblers
# -----------------------------------------------------------------------------

SUPPORTED_ASSEMBLERS = {
    "flye",
    "goldrush",
    "ntlink",
    "verkko",
}

UNKNOWN_ASSEMBLERS = (
    set(ACTIVE_ASSEMBLERS)
    - SUPPORTED_ASSEMBLERS
)

if UNKNOWN_ASSEMBLERS:
    raise ValueError(
        "Unsupported assembler names in configuration: "
        f"{sorted(UNKNOWN_ASSEMBLERS)}"
    )


# -----------------------------------------------------------------------------
# Read portable sample table
# -----------------------------------------------------------------------------

def read_sample_sheet(path):
    """
    Read sample metadata.

    samples.tsv deliberately contains NO file paths.

    Required columns:
        sample
        technology

    Example:
        HG002   ont
        HG002   pb
    """

    path = Path(path)

    if not path.is_absolute():
        path = PROJECT_DIR / path

    if not path.exists():
        raise FileNotFoundError(
            f"Missing sample sheet: {path}"
        )

    rows = []
    seen_pairs = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        required_columns = {
            "sample",
            "technology",
        }

        observed_columns = set(
            reader.fieldnames or []
        )

        missing_columns = (
            required_columns
            - observed_columns
        )

        if missing_columns:
            raise ValueError(
                "Sample sheet is missing columns: "
                f"{sorted(missing_columns)}"
            )

        for row in reader:
            sample = row["sample"].strip()
            technology = row["technology"].strip()

            if not sample or not technology:
                raise ValueError(
                    f"Invalid empty value in sample sheet row: {row}"
                )

            if technology not in {"ont", "pb"}:
                raise ValueError(
                    "Unsupported sequencing technology: "
                    f"{technology}"
                )

            pair = (
                sample,
                technology,
            )

            if pair in seen_pairs:
                raise ValueError(
                    "Duplicate sample/technology combination: "
                    f"{sample} {technology}"
                )

            seen_pairs.add(pair)

            rows.append(
                {
                    "sample": sample,
                    "technology": technology,
                }
            )

    if not rows:
        raise ValueError(
            f"Sample sheet is empty: {path}"
        )

    return rows


SAMPLE_ROWS = read_sample_sheet(
    SAMPLE_SHEET
)

SAMPLES = list(
    dict.fromkeys(
        row["sample"]
        for row in SAMPLE_ROWS
    )
)

TECHNOLOGIES = list(
    dict.fromkeys(
        row["technology"]
        for row in SAMPLE_ROWS
    )
)

SAMPLE_TECH_PAIRS = list(
    dict.fromkeys(
        (
            row["sample"],
            row["technology"],
        )
        for row in SAMPLE_ROWS
    )
)

SAMPLE_TECH_SET = set(
    SAMPLE_TECH_PAIRS
)


# -----------------------------------------------------------------------------
# Source input path construction
# -----------------------------------------------------------------------------

def _check_sample_technology(
    sample,
    technology
):
    key = (
        sample,
        technology,
    )

    if key not in SAMPLE_TECH_SET:
        raise ValueError(
            "No dataset defined for "
            f"sample={sample}, "
            f"technology={technology}"
        )



def source_fastq_for_values(
    sample,
    technology
):
    """
    Return the final normalized chromosome-21 FASTQ.

    The assembler workflow starts directly from already extracted
    and coverage-normalized chromosome-21 reads.

    Examples:
        HG002 + ont
        -> HG002/HG002.ont.chr21.primary.mapq20.30x.fastq.gz

        HG002 + pb
        -> HG002/HG002.hifi.chr21.primary.mapq20.30x.fastq.gz
    """

    _check_sample_technology(
        sample,
        technology
    )

    if INPUT_MODE != "chr21_fastq":
        raise ValueError(
            f"Unsupported input_mode: {INPUT_MODE}. "
            "Expected 'chr21_fastq'."
        )

    filename_technology = {
        "ont": "ont",
        "pb": "hifi",
    }[technology]

    filename = (
        f"{sample}."
        f"{filename_technology}."
        "chr21.primary.mapq20.30x.fastq.gz"
    )

    return str(
        INPUT_ROOT
        / sample
        / filename
    )


def source_fastq_for_values(
    sample,
    technology
):
    """
    Return the ORIGINAL source FASTQ.

    Local mode:
        input_root/HG002.ont.fastq.gz

    Server mode:
        input_root/HG002.ont.30x.fastq.gz

    The filename is derived automatically.
    It is never stored in samples.tsv.
    """

    _check_sample_technology(
        sample,
        technology
    )

    if INPUT_MODE != "chr21_fastq":
        raise ValueError(
        f"Unsupported input_mode: {INPUT_MODE}. "
        "Expected 'chr21_fastq'."
    )


    filename_technology = {"ont": "ont","pb": "hifi",}[technology]

    filename = (f"{sample}." f"{filename_technology}.""chr21.primary.mapq20.30x.fastq.gz")

    return str(INPUT_ROOT
        / sample
        / filename)


def source_fastq_for(wildcards):
    """
    Return the source FASTQ for a Snakemake wildcard object.

    This helper will be used by the future Chr21 extraction rule.
    """

    return source_fastq_for_values(
        wildcards.sample,
        wildcards.technology,
    )


# -----------------------------------------------------------------------------
# Canonical chromosome-21 paths
# -----------------------------------------------------------------------------

def generated_chr21_fastq_for_values(
    sample,
    technology
):
    """
    Canonical Chr21 FASTQ generated by server preprocessing.
    """

    _check_sample_technology(
        sample,
        technology
    )

    return (
        f"{CHR21_OUTPUT_DIR}/"
        f"{sample}.{technology}.fastq.gz"
    )


# -----------------------------------------------------------------------------
# Shared assembler input helpers
# -----------------------------------------------------------------------------

def assembler_fastq_for_values(
    sample,
    technology
):
    """
    Return the final normalized Chr21 FASTQ consumed by assemblers.

    Chr21 FASTQ is now the production workflow entry point.
    """

    _check_sample_technology(
        sample,
        technology
    )

    return source_fastq_for_values(
        sample,
        technology
    )


# Preserve the existing helper API used by Flye, GoldRush and ntLink.
def fastq_for(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        wildcards.technology,
    )


def ont_fastq_for_sample(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        "ont",
    )


def pb_fastq_for_sample(wildcards):
    return assembler_fastq_for_values(
        wildcards.sample,
        "pb",
    )


# -----------------------------------------------------------------------------
# Technology-specific assembler parameters
# -----------------------------------------------------------------------------

def flye_read_option(wildcards):
    if wildcards.technology == "ont":
        return "--nano-raw"

    if wildcards.technology == "pb":
        return "--pacbio-hifi"

    raise ValueError(
        "Unsupported sequencing technology: "
        f"{wildcards.technology}"
    )


# -----------------------------------------------------------------------------
# Expected workflow outputs
# -----------------------------------------------------------------------------

EXPECTED_FINAL_OUTPUTS = []


# Flye
if "flye" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            f"results/flye/{sample}.{technology}.fasta"
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# GoldRush
if "goldrush" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            f"results/goldrush/{sample}.{technology}.fasta"
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# ntLink
if "ntlink" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            f"results/ntlink/{sample}.{technology}.fasta"
            for sample, technology
            in SAMPLE_TECH_PAIRS
        ]
    )


# Verkko
if "verkko" in ACTIVE_ASSEMBLERS:
    EXPECTED_FINAL_OUTPUTS.extend(
        [
            f"results/verkko/{sample}.fasta"
            for sample
            in SAMPLES
        ]
    )


# -----------------------------------------------------------------------------
# Wildcard constraints
# -----------------------------------------------------------------------------

wildcard_constraints:
    sample="|".join(SAMPLES),
    technology="|".join(TECHNOLOGIES)


#####################################################################################
# ASSEMBLER RESOURCE SETTINGS
#
# Values come from config/local.yaml or config/server.yaml.
# Snakemake treats mem_mb as total memory requested by one job.
#####################################################################################

FLYE_THREADS = int(
    RESOURCE_CONFIG.get(
        "flye_threads",
        DEFAULT_THREADS
    )
)

GOLDRUSH_THREADS = int(
    RESOURCE_CONFIG.get(
        "goldrush_threads",
        DEFAULT_THREADS
    )
)

NTLINK_THREADS = int(
    RESOURCE_CONFIG.get(
        "ntlink_threads",
        DEFAULT_THREADS
    )
)

VERKKO_THREADS = int(
    RESOURCE_CONFIG.get(
        "verkko_threads",
        DEFAULT_THREADS
    )
)

FLYE_MEM_MB = int(
    RESOURCE_CONFIG.get(
        "flye_mem_mb",
        12000
    )
)

GOLDRUSH_MEM_MB = int(
    RESOURCE_CONFIG.get(
        "goldrush_mem_mb",
        12000
    )
)

GOLDRUSH_SHM_SIZE = str(RESOURCE_CONFIG.get("goldrush_shm_size", "2g"))

NTLINK_MEM_MB = int(
    RESOURCE_CONFIG.get(
        "ntlink_mem_mb",
        8000
    )
)

VERKKO_MEM_MB = int(
    RESOURCE_CONFIG.get(
        "verkko_mem_mb",
        12000
    )
)

VERKKO_LOCAL_MEMORY_GB = int(
    RESOURCE_CONFIG.get(
        "verkko_local_memory_gb",
        max(1, VERKKO_MEM_MB // 1000)
    )
)
