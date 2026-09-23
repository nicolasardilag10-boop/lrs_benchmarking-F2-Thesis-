#####################################################################################
# FLYE ASSEMBLY
#####################################################################################

rule flye_assembly:
    input:
        validation=VALIDATION_OK,
        fastq=fastq_for

    output:
        assembly="results/flye/{sample}.{technology}.fasta",
        done=PROGRESS_DIR + "/flye/{sample}.{technology}.done"

    threads:
        FLYE_THREADS

    resources:
        mem_mb=FLYE_MEM_MB

    params:
        image=DOCKER_IMAGES["flye2"],
        outdir=lambda wc: str(PROJECT_DIR / "results" / "work" / "flye" / f"{wc.sample}.{wc.technology}"),
        read_option=flye_read_option,
        genome_size=GENOME_SIZE,
        running=lambda wc: f"{PROGRESS_DIR}/flye/{wc.sample}.{wc.technology}.running",
        failed=lambda wc: f"{PROGRESS_DIR}/flye/{wc.sample}.{wc.technology}.failed"

    log:
        LOG_DIR + "/flye/{sample}.{technology}.log"

    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}" "$(dirname "{output.assembly}")" "$(dirname "{output.done}")" "$(dirname "{log}")"
        rm -f "{output.done}" "{params.running}" "{params.failed}"

        printf "STATUS=RUNNING\nsample=%s\ntechnology=%s\nstarted=%s\n" \
            "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{params.running}"

        trap 'rc=$?; rm -f "{params.running}"; if [[ $rc -ne 0 ]]; then printf "STATUS=FAILED\nsample=%s\ntechnology=%s\nfailed=%s\n" "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{params.failed}"; fi; exit $rc' EXIT

        exec > "{log}" 2>&1

        echo "FLYE START"
        echo "Image:   {params.image}"
        echo "Input:   {input.fastq}"
        echo "Output:  {output.assembly}"
        echo "Threads: {threads}"

        docker run --rm \
            --cpus {threads} \
            --memory {resources.mem_mb}m \
            --user "$(id -u):$(id -g)" \
            -e HOME=/tmp -e TMPDIR=/tmp \
            --workdir "{PROJECT_DIR_STR}" \
            {DOCKER_MOUNTS} \
            "{params.image}" \
            flye {params.read_option} "{input.fastq}" \
            --genome-size "{params.genome_size}" \
            --threads {threads} \
            --out-dir "{params.outdir}"

        if [[ ! -s "{params.outdir}/assembly.fasta" ]]; then
            echo "ERROR: Flye did not produce {params.outdir}/assembly.fasta"
            exit 1
        fi

        cp "{params.outdir}/assembly.fasta" "{output.assembly}"
        test -s "{output.assembly}"

        printf "STATUS=PASS\nsample=%s\ntechnology=%s\ncompleted=%s\n" \
            "{wildcards.sample}" "{wildcards.technology}" "$(date -Is)" > "{output.done}"

        rm -f "{params.running}" "{params.failed}"
        trap - EXIT
        """
