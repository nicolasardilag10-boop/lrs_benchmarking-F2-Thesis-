# ************************************************************************************************
#
# assembly_quality_quast.smk
#
# Whole-genome assembly quality assessment with QUAST-LG
#
# Evaluates existing Flye, GoldRush, Verkko and ntLink assembly outputs.
#
# Constitution: Articles I, II, III, V, VI and VII
#
# ************************************************************************************************

import os

# ************************************************************************************************
# Repository root
# ************************************************************************************************

CWD = os.getcwd()  # From everywhere the smk is run, go to the repository root

print("Current working directory: " + CWD)

# ************************************************************************************************
# QUAST configuration
# ************************************************************************************************

QUAST_VERSION = "5.3.0"

DOCKER_QUAST = (
    "quay.io/biocontainers/"
    "quast:5.3.0--py313pl5321h5ca1c30_2")

print("QUAST version:", QUAST_VERSION)
print("QUAST Docker image: " + DOCKER_QUAST)

# ************************************************************************************************
# Reference configuration
# ************************************************************************************************

DEFAULT_REFERENCE = os.path.expanduser(
    "~/smb/Analyses/Reference_sequence/hg38_KGGM/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_"
    "MAP2K3_KMT2C_KCNJ18.fasta")

REFERENCE = os.path.abspath(os.path.expanduser(config.get("reference", DEFAULT_REFERENCE)))
REFERENCE_DIR = os.path.dirname(REFERENCE)

print("Reference genome: " + REFERENCE)

# ************************************************************************************************
# Discover existing assembly FASTAs
# ************************************************************************************************

ASSEMBLERS_FOUND, DATASETS_FOUND = glob_wildcards(CWD + r"/assemblies/{assembler}/{dataset}/assembly.fasta")

print("Raw assemblers found:", ASSEMBLERS_FOUND)
print("Raw datasets found:", DATASETS_FOUND)

ALLOWED_OUTPUTS = {
    "flye",
    "goldrush",
    "verkko",
    "ntlink",}

ASSEMBLIES = sorted(
    {
        (assembler, dataset)
        for assembler, dataset
        in zip(ASSEMBLERS_FOUND, DATASETS_FOUND)
        if assembler.lower() in ALLOWED_OUTPUTS})

if not ASSEMBLIES:
	raise ValueError(
		"No completed assembly FASTAs were discovered under "
		"assemblies/{assembler}/{dataset}/assembly.fasta")

print("Discovered assemblies:")

for assembler, dataset in ASSEMBLIES:
    print("  " + assembler + "\t" + dataset)

# ************************************************************************************************
# QUAST targets
# ************************************************************************************************

OUTPUT = [
    f"assembly_quality/quast/{assembler}/{dataset}/report.tsv"
    for assembler, dataset in ASSEMBLIES]

print("QUAST targets:")

for target in OUTPUT:
    print("  " + target)

rule all:
    input:
        OUTPUT


# ************************************************************************************************
# QUAST-LG rule
# ************************************************************************************************

rule quast_assembly:

    input:
        assembly = "assemblies/{assembler}/{dataset}/assembly.fasta",
        reference = REFERENCE

    output:
        report = (
            "assembly_quality/quast/"
        	"{assembler}/{dataset}/report.tsv")

    params:
        outdir = (
            "assembly_quality/quast/"
            "{assembler}/{dataset}")

    log:
        ("assembly_quality/quast/"
        "{assembler}/{dataset}/quast.log")

    threads: 16

    resources:
        mem_gb = 128

    message:
        "Evaluating {wildcards.assembler} {wildcards.dataset} with QUAST-LG"

    shell:
        r"""
        set -eo pipefail

        docker run --rm \
            --cpus {threads} \
            -m {resources.mem_gb}g \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            -v {CWD}:{CWD} \
            -v {REFERENCE_DIR}:{REFERENCE_DIR}:ro \
            --workdir {CWD} \
            {DOCKER_QUAST} \
            /bin/bash -c '
                set -eo pipefail

                mkdir -p "{params.outdir}"

                (
                    printf "Container ID:\t"
                    hostname

                    printf "Start time:\t"
                    date -Is

                    echo "Assembler: {wildcards.assembler}"
                    echo "Dataset: {wildcards.dataset}"
                    echo "Assembly: {input.assembly}"
                    echo "Reference: {input.reference}"
                    echo "Threads: {threads}"
                    echo "Memory: {resources.mem_gb} GB"

                    echo
                    echo "QUAST version:"
                    quast.py --version

                    echo
                    echo "Running QUAST-LG"

                    quast.py \
                        --large \
                        --threads {threads} \
                        --min-contig 500 \
                        --reference "{input.reference}" \
                        --output-dir "{params.outdir}" \
                        "{input.assembly}"

                    echo
                    printf "End time:\t"
                    date -Is

                ) > "{log}" 2>&1

                [[ -s "{output.report}" ]] || {{
                    echo "ERROR: QUAST report.tsv is missing or empty" >> "{log}"
                    exit 101
                }}
            '
        """