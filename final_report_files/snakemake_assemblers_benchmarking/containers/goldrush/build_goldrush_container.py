#!/usr/bin/env python3

# =============================================================================
# GOLDRUSH CONTAINER BUILD AND VALIDATION
# =============================================================================
#
# Purpose:
#   Build the local GoldRush Docker image used for validation.
#
# Important:
#   This script DOES NOT run a genome assembly.
#
# Workflow:
#
#   GoldRush Dockerfile
#          ↓
#   docker build
#          ↓
#   lrs-goldrush:1.2.2
#          ↓
#   check executable
#          ↓
#   goldrush help
#
# =============================================================================


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


# Directory containing:
#
# containers/
# ├── flye2/
# ├── goldrush/
# ├── ntlink/
# └── verkko2/

CONTAINERS = PROJECT / "containers"


# GoldRush Docker build context:
#
# containers/goldrush/
# └── Dockerfile

goldrush_context = CONTAINERS / "goldrush"


# Directory where Docker build logs are stored.

LOGDIR = PROJECT / "container_build_logs"

LOGDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# STEP 3 — CHECK THAT DOCKER IS AVAILABLE
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
# STEP 4 — START THE GOLDRUSH BUILD SECTION
# =============================================================================

print()
print("=" * 60)
print("BUILDING GOLDRUSH")
print("=" * 60)


# =============================================================================
# STEP 5 — DEFINE GOLDRUSH VARIABLES
# =============================================================================

# Docker image that will be created locally.

goldrush_image = "lrs-goldrush:1.2.2"


# File where Docker build output will be saved.

goldrush_log = LOGDIR / "goldrush.build.log"


# =============================================================================
# STEP 6 — CREATE THE DOCKER BUILD COMMAND
# =============================================================================

command = [
    "docker",
    "build",
    "--progress=plain",
    "--pull",
    "--tag",
    goldrush_image,
    str(goldrush_context),
]


print("Command:")
print(" ".join(command))




# =============================================================================
# STEP 7 — OPEN THE LOG FILE
# =============================================================================

with goldrush_log.open("w") as log_file:

    log_file.write("Starting GoldRush build\n")


    # =========================================================================
    # STEP 8 — START THE DOCKER BUILD
    # =========================================================================

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


    # =========================================================================
    # STEP 9 — READ DOCKER OUTPUT WHILE IT RUNS
    # =========================================================================

    for line in process.stdout:

        print(
            line,
            end="",
        )

        log_file.write(
            line
        )


    # =========================================================================
    # STEP 10 — WAIT FOR DOCKER TO FINISH
    # =========================================================================

    return_code = process.wait()


# =============================================================================
# STEP 11 — CHECK WHETHER THE BUILD SUCCEEDED
# =============================================================================

if return_code != 0:

    print(
        "ERROR: GoldRush Docker image build failed."
    )

    sys.exit(
        return_code
    )


print(
    "SUCCESS: GoldRush Docker image created."
)

# =============================================================================
# STEP 12 — CHECK GOLDRUSH EXECUTABLE INSIDE THE CONTAINER
# =============================================================================

#Command to verift run time
runtime_command = [
    "docker",
    "run",
    "--rm",
    goldrush_image,
    "bash",
    "-lc",
    "command -v goldrush",
]

#Check the subprocess is true
runtime_check = subprocess.run(
    runtime_command,
    capture_output=True,
    text=True,
)

#Verify its running
if runtime_check.returncode != 0:
    print("ERROR: GoldRush executable was not found inside the container.")
    print(runtime_check.stderr)
    sys.exit(runtime_check.returncode)


print(
    "GoldRush executable found:"
)

print(
    runtime_check.stdout.strip()
)


# =============================================================================
# STEP 13 — CHECK THAT GOLDRUSH CAN START
# =============================================================================

help_command=[
    "docker",
    "run",
    "--rm",
    goldrush_image,
    "goldrush",
    "help",
]

# Execute "goldrush help" inside a temporary container
# and capture both standard output and errors
help_check = subprocess.run(
    help_command,
    capture_output=True,
    text=True
)

# A non-zero return code means GoldRush could not start correctly.
if help_check.returncode != 0:
    print("ERROR: GoldRush help command failed.")
    print(help_check.stderr)
    sys.exit(help_check.returncode)

# Python reaches this point only if GoldRush returned code 0.
print("GoldRush started successfully.")
print(help_check.stdout.strip())

# =============================================================================
# STEP 14 — FINAL VALIDATION SUMMARY
# =============================================================================

# =============================================================================
# STEP 14 — FINAL VALIDATION SUMMARY
# =============================================================================

print()
print("=" * 60)
print("GOLDRUSH CONTAINER VALIDATION PASSED")
print("=" * 60)

print()
print("Validation results:")
print("  Docker available            [PASS]")
print("  Docker image built          [PASS]")
print("  GoldRush executable found   [PASS]")
print("  GoldRush command executed   [PASS]")

print()
print("Docker image:")
print(
    f"  {goldrush_image}" #print directly this string
)

print()
print("Docker build context:")
print(
    f"  {goldrush_context}"
)

print()
print("Build log:")
print(
    f"  {goldrush_log}"
)

print()
print("No genome assembly was executed.")
print()