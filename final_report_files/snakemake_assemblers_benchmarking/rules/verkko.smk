#####################################################################################
# VERKKO HYBRID ASSEMBLY
# Combines PacBio HiFi and ONT reads from the same sample.
#####################################################################################

rule verkko_assembly:
    input:
        validation=VALIDATION_OK,
        hifi=pb_fastq_for_sample,
        ont=ont_fastq_for_sample

    output:
        assembly="results/verkko/{sample}.fasta",
        done=PROGRESS_DIR + "/verkko/{sample}.done"

    threads:
        VERKKO_THREADS

    resources:
        mem_mb=VERKKO_MEM_MB

    params:
        image=DOCKER_IMAGES["verkko2"],
        workdir=lambda wc: str(PROJECT_DIR / "results" / "work" / "verkko" / wc.sample),
        memory_gb=VERKKO_LOCAL_MEMORY_GB,
        running=lambda wc: f"{PROGRESS_DIR}/verkko/{wc.sample}.running",
        failed=lambda wc: f"{PROGRESS_DIR}/verkko/{wc.sample}.failed"

    log:
        LOG_DIR + "/verkko/{sample}.log"

    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.workdir}" "$(dirname "{output.assembly}")" "$(dirname "{output.done}")" "$(dirname "{log}")"
        rm -f "{output.done}" "{params.running}" "{params.failed}"

        printf "STATUS=RUNNING\nsample=%s\nstarted=%s\n" \
            "{wildcards.sample}" "$(date -Is)" > "{params.running}"

        trap 'rc=$?; rm -f "{params.running}"; if [[ $rc -ne 0 ]]; then printf "STATUS=FAILED\nsample=%s\nfailed=%s\n" "{wildcards.sample}" "$(date -Is)" > "{params.failed}"; fi; exit $rc' EXIT

        exec > "{log}" 2>&1

        echo "VERKKO START"
        echo "Image:   {params.image}"
        echo "HiFi:    {input.hifi}"
        echo "ONT:     {input.ont}"
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
            verkko -d "{params.workdir}" \
            --hifi "{input.hifi}" \
            --nano "{input.ont}" \
            --local \
            --local-cpus {threads} \
            --local-memory {params.memory_gb}

        if [[ ! -s "{params.workdir}/assembly.fasta" ]]; then
            echo "ERROR: Verkko did not produce {params.workdir}/assembly.fasta"
            exit 1
        fi

        cp "{params.workdir}/assembly.fasta" "{output.assembly}"
        test -s "{output.assembly}"

        printf "STATUS=PASS\nsample=%s\ncompleted=%s\n" \
            "{wildcards.sample}" "$(date -Is)" > "{output.done}"

        rm -f "{params.running}" "{params.failed}"
        trap - EXIT
        """
