# Runs QUAST-LG over existing Flye, GoldRush and Verkko assemblies.
# Doesn't build assemblies, just finds assembly.fasta files other pipelines already produced and scores them.
# Run this from inside the production clone on the shared server, not from a random checkout (see CWD below).

import os


# ************************************************************************************************
# Working repository
# ************************************************************************************************

# CWD is the container --workdir and the base the relative OUTPUT/log paths resolve against.
CWD = os.path.abspath(os.getcwd())

# Flye and GoldRush/Verkko assemblies are located in different collaborators filesystems
YU_ROOT = "/data/genmedbfx/yu_j/lrs_benchmarking"
STOIBER_ROOT = "/home/stoiber_l/smbshare/lrs_benchmarking"

print("Current working directory: " + CWD)
print("Yu repository: " + YU_ROOT)
print("Stoiber assembly repository: " + STOIBER_ROOT)


# ************************************************************************************************
# QUAST configuration
# ************************************************************************************************

#Qiast amd docker image
QUAST_VERSION = "5.3.0"

DOCKER_QUAST = "quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"

print("QUAST version: " + QUAST_VERSION)
print("QUAST Docker image: " + DOCKER_QUAST)


# ************************************************************************************************
# Genome reference
# ************************************************************************************************

# GRCh38 in the Schilling path
DEFAULT_REFERENCE = (
    "/data/genmedbfx/schilling_m/repos/lrs_benchmarking/ref/"
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_"
    "MAP2K3_KMT2C_KCNJ18.fasta"
)

# Override with snakemake --config reference
RAW_REFERENCE = config.get("reference", DEFAULT_REFERENCE)

REFERENCE = os.path.abspath(os.path.expanduser(RAW_REFERENCE))

# Fail until here if the file reference was not found.
if not os.path.isfile(REFERENCE):
    raise ValueError("Reference genome does not exist: " + REFERENCE)

# Bind-mounted read-only into the container so QUAST can read it
REFERENCE_DIR = os.path.dirname(REFERENCE)

print("Reference genome: " + REFERENCE)


# ************************************************************************************************
# Assembly locations
# ************************************************************************************************

# Each assembler root holds one subdir per dataset with an assembly.fasta inside, e.g. <root>/HG002.ont.30x/assembly.fasta.
ASSEMBLY_ROOTS = {
    "flye": os.path.join(YU_ROOT, "assemblies", "flye"),
    "goldrush": os.path.join(STOIBER_ROOT, "assemblies", "goldrush"),
    "verkko": os.path.join(STOIBER_ROOT, "assemblies", "verkko"),
}

for assembler, root in ASSEMBLY_ROOTS.items():
    print(f"{assembler} assembly root: {root}")


# ************************************************************************************************
# Production filtering
# ************************************************************************************************

# Dataset for test the actual coding dev/debug assembly
TEST_MARKERS = (".1k", ".chr21", ".localtest", "smoke")

# The three GIAB samples benchmarked in this study.
PRODUCTION_SAMPLES = {"hg002", "hg003", "hg004"}


# Decides whether a discovered (assembler, dataset) pair should be scored by QUAST.
def is_production_assembly(assembler, dataset):

    assembler = assembler.lower()
    dataset = dataset.lower()

    if any(marker in dataset for marker in TEST_MARKERS):
        return False

    # Verkko is hybrid, one assembly per GIAB sample.
    if assembler == "verkko":
        return dataset in PRODUCTION_SAMPLES

    # Flye and GoldRush are technology-specific, e.g. HG002.ont.30x, HG002.pb.30x.
    if assembler in {"flye", "goldrush"}:
        return ".30x" in dataset

    # Unknown assembler, exclude until there's an explicit rule for it.
    return False


# ************************************************************************************************
# Discover completed assembly FASTAs
# ************************************************************************************************

# Scan disk instead of hard-coding datasets, so the DAG follows whatever assemblies are actually finished.
ASSEMBLIES = []

for assembler, root in ASSEMBLY_ROOTS.items():

    # A root can legitimately be missing (share not mounted, no finished runs yet) — warn and skip, don't fail.
    if not os.path.isdir(root):
        print("WARNING: assembly directory does not exist: " + root)
        continue

    # Finds every dataset subdir that already has an assembly.fasta; in-progress ones are skipped automatically.
    datasets, = glob_wildcards(os.path.join(root, "{dataset}", "assembly.fasta"))

    for dataset in datasets:

        if is_production_assembly(assembler, dataset):
            ASSEMBLIES.append((assembler, dataset))


# Organzie the assemblers ordering in an specific way for across re-runs.
ASSEMBLIES = sorted(set(ASSEMBLIES))


if not ASSEMBLIES:
    raise ValueError(
        "No completed production Flye, GoldRush or Verkko "
        "assembly FASTAs were discovered.")


print("Discovered production assemblies:")

for assembler, dataset in ASSEMBLIES:
    print("  " + assembler + "\t" + dataset)


# ************************************************************************************************
# Resolve assembly FASTA
# ************************************************************************************************

# Input function: given wildcards, returns the actual assembly.fasta path
def assembly_path(wildcards):

    root = ASSEMBLY_ROOTS[wildcards.assembler]

    path = os.path.join(root, wildcards.dataset, "assembly.fasta")

    # Check twice since the earlier filesystem
    if not os.path.isfile(path):
        raise ValueError("Assembly FASTA does not exist: " + path)

    return path


# ************************************************************************************************
# QUAST targets
# ************************************************************************************************

# One report.tsv per discovered production (assembler, dataset) pair; must match output.quast_tsv's pattern
OUTPUT = [
    f"assembly_quality/quast/{assembler}/{dataset}/report.tsv"
    for assembler, dataset in ASSEMBLIES]


print("QUAST targets:")

for target in OUTPUT:
    print("  " + target)


# Default target: schedules one quast_assembly job per production assembly.
rule all:
    input:
        OUTPUT


# ************************************************************************************************
# QUAST-LG
# ************************************************************************************************

# Shell string below goes through str.format(), so any literal "{" or "}" meant for bash (the `|| {{ }}` guards) must be doubled.
rule quast_assembly:
    input:
        assembly=assembly_path,
        reference=REFERENCE

    output:
        quast_tsv="assembly_quality/quast/{assembler}/{dataset}/report.tsv"

    params:
        outdir="assembly_quality/quast/{assembler}/{dataset}"

    log:
        "assembly_quality/quast/{assembler}/{dataset}/quast.log"

    threads: 16

    resources:
        mem_gb=128

    message:
        "Evaluating {wildcards.assembler} {wildcards.dataset} with QUAST-LG"

    shell:
        r"""
        set -eo pipefail

        mkdir -p "{params.outdir}"
        : > "{log}"

        # The assembly may be a symlink to somewhere outside the bind-mounted roots, so dereference it on the host and stage a real copy under YU_ROOT.
        STAGED_ASSEMBLY="{params.outdir}/staged_assembly.fasta"
        trap 'rm -f "$STAGED_ASSEMBLY"' EXIT
        cp -L "{input.assembly}" "$STAGED_ASSEMBLY"

        # YU_ROOT is rw since the staged copy and QUAST's output live there; STOIBER_ROOT and the reference dir are root.
        # --workdir matches CWD so the relative outdir/log paths resolve the same inside and outside the container (CWD must sit under YU_ROOT).
        docker run --rm \
            --cpus {threads} \
            -m {resources.mem_gb}g \
            --tmpfs /tmp:size=50g,exec \
            -u $UID:$(id -g) \
            -v "{YU_ROOT}:{YU_ROOT}" \
            -v "{STOIBER_ROOT}:{STOIBER_ROOT}:ro" \
            -v "{REFERENCE_DIR}:{REFERENCE_DIR}:ro" \
            --workdir "{CWD}" \
            {DOCKER_QUAST} \
            /bin/bash -c '
                set -eo pipefail

                mkdir -p "{params.outdir}"

                printf "Container hostname:\t"
                hostname

                printf "Start time:\t"
                date -Is

                echo "Assembler: {wildcards.assembler}"
                echo "Dataset: {wildcards.dataset}"
                echo "Assembly (original): {input.assembly}"
                echo "Assembly (staged): {params.outdir}/staged_assembly.fasta"
                echo "Reference: {input.reference}"
                echo "Threads: {threads}"
                echo "Memory: {resources.mem_gb} GB"

                echo
                echo "QUAST version:"
                quast.py --version

                echo
                echo "Running QUAST-LG"

                # Runs against the staged copy, never the original symlinked assembly, so the container doesn't need to resolve links outside its mounts.
                quast.py \
                    --large \
                    --threads {threads} \
                    --min-contig 500 \
                    --reference "{input.reference}" \
                    --output-dir "{params.outdir}" \
                    "{params.outdir}/staged_assembly.fasta"

                echo
                printf "End time:\t"
                date -Is

                # QUAST can exit 0 with a broken report (e.g. disk ran out mid-write), so double-check it here. Exit 101 marks it as our own check, not QUAST's.
                [[ -s "{output.quast_tsv}" ]] || {{
                    echo "ERROR: QUAST report.tsv is missing or empty"
                    exit 101
                }}

                grep -q "^N50" "{output.quast_tsv}" || {{
                    echo "ERROR: N50 not found in QUAST report"
                    exit 101
                }}
            ' >> "{log}" 2>&1
        """
