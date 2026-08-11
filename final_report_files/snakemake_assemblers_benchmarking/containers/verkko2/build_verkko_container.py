#!/usr/bin/env python3


# =============================================================================
# STEP 1 — IMPORT PYTHON MODULES
# =============================================================================

from pathlib import Path

import subprocess
import sys


# =============================================================================
# STEP 2 — DEFINE PROJECT DIRECTORIES
# =============================================================================

PROJECT = (
    Path.home()
    / "lrs_benchmarking"
    / "final_report_files"
    / "snakemake_assemblers_benchmarking"
)

CONTAINERS = PROJECT / "containers"

verkko_context = CONTAINERS / "verkko2"

LOGDIR = PROJECT / "container_build_logs"
LOGDIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# STEP 3 — DEFINE VERKKO VARIABLES
# =============================================================================

verkko_image = "lrs-verkko2:2.3.2"

verkko_log = LOGDIR / "verkko2.build.log"


# =============================================================================
# STEP 4 — CHECK THAT DOCKER IS AVAILABLE
# =============================================================================

print("Checking Docker...")

docker_check = subprocess.run(
    ["docker", "info"],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)

if docker_check.returncode != 0:
    print("ERROR: Docker is not available.")
    print("Start Docker Desktop and enable WSL integration.")
    sys.exit(1)

print("Docker is ready.")


# =============================================================================
# STEP 5 — CREATE THE DOCKER BUILD COMMAND
# =============================================================================

command = [
    "docker",
    "build",
    "--progress=plain",
    "--pull",
    "--tag",
    verkko_image,
    str(verkko_context),
]


# =============================================================================
# STEP 6 — DISPLAY THE BUILD COMMAND
# =============================================================================

print()
print("=" * 60)
print("BUILDING VERKKO2")
print("=" * 60)

print("Command:")
print(" ".join(command))


# =============================================================================
# STEP 7 — RUN THE VERKKO2 DOCKER BUILD
# =============================================================================

with verkko_log.open("w") as log_file:

    log_file.write("Starting Verkko2 Docker build\n")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    # Open log
    #     ↓
    # Start Docker build
    #     ↓
    # Read Docker output line by line
    #     ├── Show it in terminal
    #     └── Save it in verkko2.build.log
    #     ↓
    # Wait for Docker to finish
    #     ↓
    # Store return code

    for line in process.stdout:
        print(line, end="")
        log_file.write(line)

    return_code = process.wait()


# =============================================================================
# STEP 8 — CHECK WHETHER THE BUILD SUCCEEDED
# =============================================================================

if return_code != 0:
    print("ERROR: Verkko2 Docker image build failed.")
    print(f"Check the build log: {verkko_log}")
    sys.exit(return_code)

print("SUCCESS: Verkko2 Docker image created.")


# =============================================================================
# STEP 9 — CHECK THAT VERKKO EXISTS INSIDE THE CONTAINER
# =============================================================================

runtime_command = [
    "docker",
    "run",
    "--rm",
    verkko_image,
    "bash",
    "-lc",
    "command -v verkko",
]

runtime_check = subprocess.run(
    runtime_command,
    capture_output=True,
    text=True,
)

if runtime_check.returncode != 0:
    print("ERROR: Verkko executable was not found inside the container.")
    print(runtime_check.stderr)
    sys.exit(runtime_check.returncode)

print("Verkko executable found:")
print(runtime_check.stdout.strip())


# =============================================================================
# STEP 10 — TEST THE VERKKO COMMAND
# =============================================================================

help_command = [
    "docker",
    "run",
    "--rm",
    verkko_image,
    "verkko",
    "--help",
]

help_check = subprocess.run(
    help_command,
    capture_output=True,
    text=True,
)

if help_check.returncode != 0:
    print("ERROR: Verkko help command failed.")
    print(help_check.stdout)
    print(help_check.stderr)
    sys.exit(help_check.returncode)

print()
print("Verkko started successfully.")
print(help_check.stdout.strip())


# =============================================================================
# STEP 11 — FINAL VALIDATION SUMMARY
# =============================================================================

print()
print("=" * 60)
print("VERKKO2 CONTAINER VALIDATION PASSED")
print("=" * 60)

print()
print("Validation results:")
print("  Docker available            [PASS]")
print("  Docker image built          [PASS]")
print("  Verkko executable found     [PASS]")
print("  Verkko command executed     [PASS]")

print()
print("Docker image:")
print(f"  {verkko_image}")

print()
print("Docker build context:")
print(f"  {verkko_context}")

print()
print("Build log:")
print(f"  {verkko_log}")

print()
print("No genome assembly was executed.")
print()