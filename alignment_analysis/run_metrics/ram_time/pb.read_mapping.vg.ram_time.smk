# ************************************************************************************************
# VG Giraffe PacBio HiFi read mapping -- RUNTIME / PEAK-RAM MEASUREMENT ONLY
#
# Re-runs the vg_map_sort rule of pb.read_mapping.vg.smk with the same Docker
# image, indexes, memory ceilings and parameters, only to record wall-clock
# time, peak RAM, CPU time and peak scratch disk. The CRAM is not kept.
#
# Rules (per dataset):
#   vg_run      -> runs the three stages of vg_map_sort, each under
#                  /usr/bin/time -v inside its container:
#                    map   (vg giraffe -> SAM)
#                    sort  (samtools sort SAM -> CRAM)
#                    index (samtools index CRAM)
#                  writing everything to a temp() scratch directory
#   vg_record   -> sums the three stages into {dataset}.ram_time.tsv
#                  (plus {dataset}.raw/ with the raw time and resource-sample files)
#
# The scratch directory is a temp() output, so Snakemake deletes it automatically
# as soon as the record rule succeeds (and also removes it if the run fails).
# Every rule writes its own log ({run,record}.log) per CONSTITUTION VI.3/VI.4.
#
# /usr/bin/time -v runs *inside* the container because Snakemake's
# `benchmark:` directive samples the host process tree and cannot see the
# container's real memory use.
#
# Configurable (--config key=value), no paths hard-coded in the rules;
# both directories must stay inside the repository root (CONSTITUTION I.1):
#   ram_time_threads      (default 64, production used 16; same as the assembler ram_time runs)
#   ram_time_dir          (default aligners_ram_time)          -- ram_time TSVs + logs
#   ram_time_scratch_dir  (default aligners_ram_time/scratch)  -- deleted automatically after each run
# Relative paths are resolved against the repository root (the working directory).
# ************************************************************************************************

import os

CWD = os.getcwd()
print("Current working directory: " + CWD)

try:
    DATASET_FILTER = config["dataset_filter"]
except (KeyError, NameError):
    DATASET_FILTER = None

#################
# RAM / time measurement settings

RAM_TIME_THREADS = int(config.get("ram_time_threads", 64))
RAM_TIME_DIR = os.path.join(CWD, config.get("ram_time_dir", "aligners_ram_time"), "vg_giraffe")
SCRATCH_DIR = os.path.join(CWD, config.get("ram_time_scratch_dir", "aligners_ram_time/scratch"), "vg_giraffe")

for _path in (RAM_TIME_DIR, SCRATCH_DIR):
    if os.path.commonpath([CWD, os.path.realpath(_path)]) != os.path.realpath(CWD):
        raise ValueError("Constitution I.1: ram_time paths must be inside the repository root: " + _path)

print("RAM/time threads: " + str(RAM_TIME_THREADS))
print("RAM/time results: " + RAM_TIME_DIR)
print("RAM/time scratch: " + SCRATCH_DIR)

#################
# Docker image, reference and indexes (identical to pb.read_mapping.vg.smk)

DOCKER_VG = "schimar/lrs-vg:v1.73.0"

LOCAL_REFERENCE = os.path.join(
    CWD,
    "reference",
    "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta",
)

RAW_REFERENCE = config.get("reference", LOCAL_REFERENCE)
REF = os.path.expanduser(RAW_REFERENCE)
if not os.path.isabs(REF):
    REF = os.path.join(CWD, REF)
REF = os.path.abspath(REF)

VG_INDEX_DIR = CWD + "/vg_index"
VG_GBZ       = VG_INDEX_DIR + "/hg38.giraffe.gbz"
VG_DIST      = VG_INDEX_DIR + "/hg38.dist"
VG_MIN       = VG_INDEX_DIR + "/hg38.longread.withzip.min"
VG_ZIPCODES  = VG_INDEX_DIR + "/hg38.longread.zipcodes"

# Background sampler for whole-container memory and scratch disk use, shared
# with the assembler ram_time workflows (see the README).
SAMPLER = os.path.join(CWD, "assemblers", "whole_genome_asm", "ram_time", "sample_container_resources.sh")
if not os.path.isfile(SAMPLER):
    raise FileNotFoundError("Resource sampler not found: " + SAMPLER)

PARSER = os.path.join(CWD, "assembly_analysis", "scripts", "metrics", "parse_assembler_time_v.py")

#####################
# Discover datasets (production filter, plus the test-subset exclusions of the
# assembler ram_time workflows)
DATASETS_FASTQ, = glob_wildcards(CWD + r"/fastq/{dataset,[A-Za-z0-9._-]+}.fastq.gz")
DATASETS = [
    dataset
    for dataset in DATASETS_FASTQ
    if ".pb." in dataset.lower()
    and ".1k" not in dataset.lower()
    and ".chr21." not in dataset.lower()
    and "localtest" not in dataset.lower()
    and "smoke" not in dataset.lower()]
if DATASET_FILTER:
    DATASETS = [d for d in DATASETS if DATASET_FILTER in d]

##############
# Targets

OUTPUT = []

OUTPUT += expand(RAM_TIME_DIR + "/{dataset}.ram_time.tsv", dataset=DATASETS)

rule all:
    input:
        OUTPUT

print("Discover datasets and create wildcards")
print(OUTPUT)

################
# Rules

rule vg_run:
    input:
        fastq    = CWD + "/fastq/{dataset}.fastq.gz",
        ref      = REF,
        gbz      = VG_GBZ,
        dist     = VG_DIST,
        min_idx  = VG_MIN,
        zipcodes = VG_ZIPCODES,

    output:
        scratch = temp(directory(SCRATCH_DIR + "/{dataset}"))

    log:
        RAM_TIME_DIR + "/{dataset}.run.log"

    message:
        "executing {rule} with output {output} and input {input}"

    threads: RAM_TIME_THREADS

    resources:
        mem_mb = 163840,

    shell:
        """
        mkdir -p "$(dirname "{log}")"

        (
            set -eo pipefail

            echo "[$(date -Is)] START vg_run {wildcards.dataset}"
            echo "Container hostname: vg-pb-ram-time"
            echo "Dataset: {wildcards.dataset}"
            echo "Read technology: PacBio HiFi (-b hifi)"
            echo "Threads: {threads}"
            echo "Scratch: {output.scratch}"

            [[ {threads} -eq {RAM_TIME_THREADS} ]] || {{
                echo "ERROR: Snakemake granted {threads} threads, expected {RAM_TIME_THREADS}."
                echo "Run with --cores {RAM_TIME_THREADS} (or more) so the measurement is comparable."
                exit 104;
            }}

            mkdir -p "{output.scratch}"

            echo "Checking /usr/bin/time, python3 and du -sb inside {DOCKER_VG} ..."
            docker run --rm --entrypoint sh "{DOCKER_VG}" \
                -c 'command -v /usr/bin/time && command -v python3 && du -sb /tmp >/dev/null' || {{
                echo "ERROR: /usr/bin/time, python3 or a du supporting -sb (needed by the record rule and the sampler) is not available inside {DOCKER_VG}"
                exit 103;
            }}

            docker stats --no-stream >/dev/null || {{
                echo "ERROR: docker stats does not work on this host; the resource sampler needs it"
                exit 105;
            }}

            # One sampler per stage, each watching that stage's named container.
            SAMPLER_PID=""
            start_sampler() {{
                bash "{SAMPLER}" "ramtime-vg-{wildcards.dataset}-$1" "{output.scratch}" "{output.scratch}/resource_samples.$1.tsv" &
                SAMPLER_PID=$!
                echo "Resource sampler for stage $1 started (pid $SAMPLER_PID)"
            }}
            stop_sampler() {{
                [[ -n "$SAMPLER_PID" ]] || return 0
                kill $SAMPLER_PID 2>/dev/null || true
                wait $SAMPLER_PID 2>/dev/null || true
                SAMPLER_PID=""
            }}
            trap stop_sampler EXIT

            TMP_SAM="{output.scratch}/{wildcards.dataset}.sam"
            CRAM="{output.scratch}/{wildcards.dataset}.cram"
            VG_STDERR="{output.scratch}/vg.stderr"

            # --- Stage 1: vg giraffe -> SAM
            # NB: capture the exit code with set +e -- otherwise set -e kills the
            # subshell at the docker line and the stderr dump below is never reached
            start_sampler map
            set +e
            docker run --rm \
                --name "ramtime-vg-{wildcards.dataset}-map" \
                --hostname "vg-pb-ram-time" \
                --tmpfs /tmp:size=50g,exec \
                -u $UID:$(id -g) \
                --cpus {threads} \
                -m 160g \
                -v {CWD}:{CWD} \
                -v {input.ref}:{input.ref}:ro \
                --entrypoint /usr/bin/time \
                {DOCKER_VG} \
                -v -o "{output.scratch}/map_time_v.txt" \
                vg giraffe \
                -t {threads} \
                -Z {input.gbz} \
                -d {input.dist} \
                -m {input.min_idx} \
                -z {input.zipcodes} \
                -b hifi \
                -f {input.fastq} \
                --output-format SAM \
                -R "ID:{wildcards.dataset}\tSM:{wildcards.dataset}" \
                -N {wildcards.dataset} \
                > "$TMP_SAM" \
                2> "$VG_STDERR"

            VG_EXIT=$?
            set -e
            stop_sampler

            echo "vg giraffe exit code: $VG_EXIT"
            if [[ "$VG_EXIT" -ne 0 ]]; then
                echo "ERROR: vg giraffe failed with exit code $VG_EXIT"
                cat "$VG_STDERR" || true
                exit 101
            fi

            [[ -s "$TMP_SAM" ]] || {{ echo "ERROR: vg giraffe produced empty/missing $TMP_SAM"; exit 101; }}

            # --- Stage 2: samtools sort SAM -> CRAM
            start_sampler sort
            docker run --rm \
                --name "ramtime-vg-{wildcards.dataset}-sort" \
                --hostname "vg-pb-ram-time" \
                --tmpfs /tmp:size=50g,exec \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus 4 \
                -m 16g \
                -v {CWD}:{CWD} \
                -v {input.ref}:{input.ref}:ro \
                --entrypoint /usr/bin/time \
                {DOCKER_VG} \
                -v -o "{output.scratch}/sort_time_v.txt" \
                samtools sort \
                -@ 4 \
                -O CRAM \
                --reference {input.ref} \
                -o "$CRAM" \
                "$TMP_SAM"
            stop_sampler

            # Scratch size with both SAM and CRAM present: the disk peak.
            PEAK_DISK=$(docker run --rm -u $UID:$(id -g) -v {CWD}:{CWD} --entrypoint du {DOCKER_VG} -sb "{output.scratch}" | cut -f1)
            printf 'epoch_seconds\\tcontainer_mem_bytes\\tscratch_disk_bytes\\n%s\\tNA\\t%s\\n' "$(date +%s)" "$PEAK_DISK" \
                > "{output.scratch}/resource_samples.after_sort.tsv"

            rm -f "$TMP_SAM"

            # --- Stage 3: samtools index CRAM
            start_sampler index
            docker run --rm \
                --name "ramtime-vg-{wildcards.dataset}-index" \
                --hostname "vg-pb-ram-time" \
                --tmpfs /tmp:size=50g,exec \
                --workdir /tmp \
                -u $UID:$(id -g) \
                --cpus 4 \
                -m 8g \
                -v {CWD}:{CWD} \
                -v {input.ref}:{input.ref}:ro \
                --entrypoint /usr/bin/time \
                {DOCKER_VG} \
                -v -o "{output.scratch}/index_time_v.txt" \
                samtools index "$CRAM" "$CRAM.crai"
            stop_sampler
            trap - EXIT

            # Merge the per-stage samples into one table (one header).
            {{
                head -n 1 "{output.scratch}/resource_samples.map.tsv"
                for stage in map sort after_sort index; do
                    tail -n +2 "{output.scratch}/resource_samples.$stage.tsv"
                done
            }} > "{output.scratch}/resource_samples.tsv"

            [[ $(awk -F '\\t' 'NR > 1 && $2 != "NA"' "{output.scratch}/resource_samples.tsv" | wc -l) -gt 0 ]] || {{
                echo "ERROR: resource sampler recorded no memory samples ({output.scratch}/resource_samples.tsv)"
                exit 106;
            }}

            [[ -s "$CRAM.crai" ]] || {{ echo "ERROR: CRAI is missing or empty"; exit 101; }}

            if [[ "$(du -b "$CRAM" | cut -f 1)" -le 64 ]]; then
                echo "ERROR: CRAM $CRAM is <=64 bytes -- vg giraffe produced no alignments; failing"
                exit 101
            fi

            for stage in map sort index; do
                [[ -s "{output.scratch}/${{stage}}_time_v.txt" ]] || {{
                    echo "ERROR: /usr/bin/time output for stage $stage is missing or empty"
                    exit 102;
                }}
            done

            echo "[$(date -Is)] END vg_run {wildcards.dataset}"

        ) > "{log}" 2>&1
        """

rule vg_record:
    input:
        scratch = SCRATCH_DIR + "/{dataset}"

    output:
        ram_time = RAM_TIME_DIR + "/{dataset}.ram_time.tsv",
        raw = directory(RAM_TIME_DIR + "/{dataset}.raw")

    log:
        RAM_TIME_DIR + "/{dataset}.record.log"

    message:
        "executing {rule} with output {output} and input {input}"

    shell:
        """
        (
            set -eo pipefail

            echo "[$(date -Is)] START vg_record {wildcards.dataset}"
            echo "Container hostname: vg-record-{wildcards.dataset}"

            # Parser is stdlib-only Python, run inside the aligner's own pinned image.
            docker run --rm \
                --tmpfs /tmp:size=50g,exec \
                --hostname vg-record-{wildcards.dataset} \
                -u $UID:$(id -g) \
                -v {CWD}:{CWD} \
                --entrypoint python3 \
                {DOCKER_VG} \
                "{PARSER}" \
                --aligner vg_giraffe \
                --sample "{wildcards.dataset}" \
                --technology pb \
                --threads {RAM_TIME_THREADS} \
                --time-file "{input.scratch}/map_time_v.txt" "{input.scratch}/sort_time_v.txt" "{input.scratch}/index_time_v.txt" \
                --resource-samples "{input.scratch}/resource_samples.tsv" \
                --output "{output.ram_time}"

            [[ $(wc -l < "{output.ram_time}") -eq 2 ]] || {{
                echo "ERROR: {output.ram_time} must contain exactly one header and one data row"
                exit 101;
            }}

            # Keep the raw time files and resource samples so every number can be re-checked.
            mkdir -p "{output.raw}"
            cp "{input.scratch}"/*_time_v.txt "{input.scratch}"/resource_samples*.tsv "{output.raw}/"

            [[ -s "{output.raw}/resource_samples.tsv" ]] || {{
                echo "ERROR: raw files were not copied to {output.raw}"
                exit 102;
            }}

            echo "[$(date -Is)] END vg_record {wildcards.dataset}"

        ) > "{log}" 2>&1
        """
