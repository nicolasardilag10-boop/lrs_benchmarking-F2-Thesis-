# ************************************************************************************************
# * Snakefile for truvari anno annotations on fn/fp VCFs
# ************************************************************************************************

import os


#Get current working directory
CWD = os.getcwd()


#Define variables
DB_DIR = "/mnt/storage/db"
FASTA = "references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
ADOTTO_TRF = "assets/tr_annotated.bed.gz"

MAMBA = "/mnt/storage/groups/genetics/VarCAD-dev/external/Miniforge3/condabin/mamba run --live-stream"

#Set FASTA_PATH
FASTA_PATH = DB_DIR + "/" + FASTA
ADOTTO_TRF_PATH = CWD + "/" + ADOTTO_TRF


# Discover existing truvari pbsv result directories with fn/fp VCFs
TRUVARI_DIRS, VCF_TYPES = glob_wildcards(CWD + "/truvari/{truvari_dir,[A-Za-z0-9\-\._]+}/{vcftype,fn|fp}.vcf.gz")

#Logging
print("Database directory: " + DB_DIR)
print("Path to reference genome: " + FASTA_PATH)
print("Current working directory: " + CWD)
print("Truvari dirs: " + ','.join(set(TRUVARI_DIRS)))


# *** Define Output
OUTPUT = []

OUTPUT = OUTPUT + expand(CWD + "/truvari/{truvari_dir}/{vcftype}.anno.vcf.gz",
    zip, truvari_dir=TRUVARI_DIRS, vcftype=VCF_TYPES)

OUTPUT = OUTPUT + expand(CWD + "/truvari/{truvari_dir}/{vcftype}.grm.jl",
    zip, truvari_dir=TRUVARI_DIRS, vcftype=VCF_TYPES)


# ************************************************************************************************

rule all:
    input: OUTPUT


# ************************************************************************************************
# Rules
# ************************************************************************************************

rule truvari_anno_gcpct:
    input:
        vcf_gz="{cwd}/truvari/{truvari_dir}/{vcftype}.vcf.gz",
        fasta=FASTA_PATH
    output:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.gcpct.vcf"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.gcpct.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
        mem_gb=8
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno gcpct \
                -r {input.fasta} \
                -o {output.vcf} \
                {input.vcf_gz}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_trf:
    input:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.gcpct.vcf",
        fasta=FASTA_PATH,
        adotto=ADOTTO_TRF_PATH
    output:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.trf.vcf"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.trf.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 4
    resources:
        mem_gb=16
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno trf \
                -i {input.vcf} \
                -o {output.vcf} \
                -r {input.adotto} \
                -f {input.fasta} \
                -t {threads}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_repmask:
    input:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.trf.vcf"
    output:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.repmask.vcf"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.repmask.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 4
    resources:
        mem_gb=16
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno repmask \
                -i {input.vcf} \
                -o {output.vcf} \
                -T {threads}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_remap:
    input:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.repmask.vcf",
        fasta=FASTA_PATH
    output:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.remap.vcf"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.remap.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
        mem_gb=8
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno remap \
                -r {input.fasta} \
                -o {output.vcf} \
                {input.vcf}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_lcr:
    input:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.remap.vcf"
    output:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.lcr.vcf"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.lcr.vcf.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
        mem_gb=8
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno lcr \
                -o {output.vcf} \
                {input.vcf}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_dpcnt:
    input:
        vcf="{cwd}/truvari/{truvari_dir}/{vcftype}.lcr.vcf"
    output:
        vcf_gz="{cwd}/truvari/{truvari_dir}/{vcftype}.anno.vcf.gz"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.anno.vcf.gz.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 2
    resources:
        mem_gb=8
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n truvari truvari anno dpcnt \
                {input.vcf} \
                | bgzip -c -@ {threads} > {output.vcf_gz}; \
            tabix -p vcf {output.vcf_gz}; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


rule truvari_anno_grm:
    """
    Run truvari anno grm separately — outputs a joblib DataFrame, not a VCF.
    Requires a BWA-indexed reference FASTA.
    """
    input:
        vcf_gz="{cwd}/truvari/{truvari_dir}/{vcftype}.vcf.gz",
        fasta=FASTA_PATH
    output:
        jl="{cwd}/truvari/{truvari_dir}/{vcftype}.grm.jl"
    log:     "{cwd}/truvari/{truvari_dir}/{vcftype}.grm.jl.log"
    message: "executing {rule} with output {output} and input {input}"
    threads: 4
    resources:
        mem_gb=16
    shell:   "umask 0027; \
        srun -p all -c {threads} --mem={resources.mem_gb}GB /bin/bash -c \" \
            printf 'Container ID:\\t'; hostname; \
            printf 'Start time:\\t'; date; \
            umask 0027; \
            {MAMBA} -n bwa bash -c ' \
                {MAMBA} -n truvari truvari anno grm \
                    -i {input.vcf_gz} \
                    -r {input.fasta} \
                    -o {output.jl} \
                    -t {threads}; \
            '; \
            printf 'End time:\\t'; date; \" \
        &> {log};"


# ************************************************************************************************
