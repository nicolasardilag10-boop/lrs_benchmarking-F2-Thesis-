#####################################################################################
# GOLDRUSH ASSEMBLY
# Produces one assembly per sample and sequencing technology.
#####################################################################################

rule goldrush_assembly:
    input:
        validation=VALIDATION_OK,
        reads=fastq_for

    output:
        assembly="results/goldrush/{sample}.{technology}.fasta",
        done=PROGRESS_DIR + "/goldrush/{sample}.{technology}.done"

    threads:
        GOLDRUSH_THREADS

    resources:
        mem_mb=GOLDRUSH_MEM_MB

    params:
        image=DOCKER_IMAGES["goldrush"],
        workdir=lambda wc: str(PROJECT_DIR / "results" / "work" / "goldrush" / f"{wc.sample}.{wc.technology}"),
        prefix=lambda wc: f"{wc.sample}_{wc.technology}_goldrush",
        reads_basename=lambda wc: f"{wc.sample}_{wc.technology}",
        min_read_length=lambda wc: 5000 if wc.technology == "ont" else 10000,
        genome_size=GENOME_SIZE,
        shm_size=GOLDRUSH_SHM_SIZE,
        running=lambda wc: f"{PROGRESS_DIR}/goldrush/{wc.sample}.{wc.technology}.running",
        failed=lambda wc: f"{PROGRESS_DIR}/goldrush/{wc.sample}.{wc.technology}.failed"

    log:
        LOG_DIR + "/goldrush/{sample}.{technology}.log"

    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.workdir}" "$(dirname "{output.assembly}")" "$(dirname "{output.done}")" "$(dirname "{log}")"
        rm -f "{output.done}" "{params.running}" "{params.failed}"

        printf "STATUS=RUNNING\nsample=%s\ntechnology=%s\nstarted=%s\n" \
            "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{params.running}"

        trap 'rc=$?; rm -f "{params.running}"; if [[ $rc -ne 0 ]]; then printf "STATUS=FAILED\nsample=%s\ntechnology=%s\nfailed=%s\n" "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{params.failed}"; fi; exit $rc' EXIT

        exec > "{log}" 2>&1

        echo "GOLDRUSH START"
        echo "Image:       {params.image}"
        echo "Input:       {input.reads}"
        echo "Output:      {output.assembly}"
        echo "Workdir:     {params.workdir}"
        echo "Threads:     {threads}"
        echo "Memory MB:   {resources.mem_mb}"
        echo "Min read bp: {params.min_read_length}"

        READS_FASTQ="{params.workdir}/{params.reads_basename}.fastq"

        if [[ ! -s "$READS_FASTQ" ]]; then
            echo "Preparing uncompressed GoldRush input..."
            gzip -cd "{input.reads}" > "$READS_FASTQ.tmp"
            mv "$READS_FASTQ.tmp" "$READS_FASTQ"
        fi

        docker run --rm \
            --cpus {threads} \
            --memory {resources.mem_mb}m \
            --shm-size "{params.shm_size}" \
            --user "$(id -u):$(id -g)" \
            -e HOME=/tmp \
            -e TMPDIR=/tmp \
            --workdir "{params.workdir}" \
            {DOCKER_MOUNTS} \
            "{params.image}" \
            goldrush run \
            reads="{params.reads_basename}" \
            G={params.genome_size} \
            t={threads} \
            m={params.min_read_length} \
            P=0 \
            p="{params.prefix}" \
            track_time=1

        FINAL_ASSEMBLY=$(find "{params.workdir}/goldrush_intermediate_files" \
            -type f -name "*ntLink-5rounds.polished.fa" -print | sort | tail -n 1)

        if [[ -z "$FINAL_ASSEMBLY" || ! -s "$FINAL_ASSEMBLY" ]]; then
            echo "ERROR: GoldRush final polished assembly was not produced."
            echo "Expected under: {params.workdir}/goldrush_intermediate_files"
            exit 1
        fi

        cp "$FINAL_ASSEMBLY" "{output.assembly}"
        test -s "{output.assembly}"

        printf "STATUS=PASS\nsample=%s\ntechnology=%s\ncompleted=%s\n" \
            "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{output.done}"

        rm -f "{params.running}" "{params.failed}"
        trap - EXIT
        """
