#####################################################################################
# INPUT AND RUNTIME VALIDATION
#
# Checks metadata, automatically resolved inputs, the server reference when needed,
# and Docker availability before any assembler is allowed to start.
#####################################################################################


rule validate_inputs:
    output:
        VALIDATION_OK

    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        errors = []

        # Sample sheet
        sample_sheet_path = Path(SAMPLE_SHEET)
        if not sample_sheet_path.is_absolute():
            sample_sheet_path = PROJECT_DIR / sample_sheet_path

        if not sample_sheet_path.exists():
            errors.append(f"Missing sample sheet: {sample_sheet_path}")

        # Input root
        if not INPUT_ROOT.exists():
            errors.append(f"Missing input_root directory: {INPUT_ROOT}")

        # Automatically resolved source FASTQs
        for row in SAMPLE_ROWS:
            sample = row["sample"]
            technology = row["technology"]
            fastq = Path(source_fastq_for_values(sample, technology))

            if not fastq.exists():
                errors.append(f"Missing source FASTQ for {sample} {technology}: {fastq}")
            elif not fastq.is_file():
                errors.append(f"Input is not a regular FASTQ file for {sample} {technology}: {fastq}")

        # Reference is required only for server/WGS mode.
        if INPUT_MODE == "whole_genome":
            if not REFERENCE:
                errors.append("Server/whole_genome mode requires 'reference' in config/server.yaml")
            elif not REFERENCE.exists():
                errors.append(f"Missing reference FASTA: {REFERENCE}")

        # Docker is required whenever an assembler is active.
        if ACTIVE_ASSEMBLERS:
            docker = shutil.which("docker")

            if not docker:
                errors.append("Docker executable was not found in PATH.")
            else:
                docker_info = subprocess.run(
                    [docker, "info"],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True
                )

                if docker_info.returncode != 0:
                    detail = docker_info.stderr.strip().splitlines()
                    detail = detail[-1] if detail else "Docker daemon is not reachable."
                    errors.append(f"Docker daemon is not reachable: {detail}")

        if errors:
            message = "\n".join(f"- {error}" for error in errors)

            raise RuntimeError(
                "\n"
                "ASSEMBLER WORKFLOW VALIDATION FAILED\n\n"
                f"Input mode: {INPUT_MODE}\n"
                f"Input root: {INPUT_ROOT}\n\n"
                "Fix these problems before running the workflow:\n\n"
                f"{message}\n"
            )

        Path(output[0]).write_text(
            "Assembler workflow validation passed.\n"
            f"input_mode={INPUT_MODE}\n"
            f"input_root={INPUT_ROOT}\n"
            f"samples={len(SAMPLES)}\n"
            f"datasets={len(SAMPLE_TECH_PAIRS)}\n"
            f"docker_images={','.join(DOCKER_IMAGES[key] for key in sorted(DOCKER_IMAGES))}\n"
        )
