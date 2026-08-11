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

ntlink_context = CONTAINERS / "ntlink"

LOGDIR = PROJECT / "container_build_logs"
LOGDIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# STEP 3 — DEFINE NTLINK VARIABLES
# =============================================================================

ntlink_image = "lrs-ntlink:1.3.11"

ntlink_log = LOGDIR / "ntlink.build.log"


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
    ntlink_image,
    str(ntlink_context),
]


# =============================================================================
# STEP 6 — DISPLAY THE BUILD COMMAND
# =============================================================================

print()
print("=" * 60)
print("BUILDING NTLINK")
print("=" * 60)

print("Command:")
print(" ".join(command))


# =============================================================================
# STEP 7 — OPEN THE BUILD LOG AND START DOCKER BUILD
# =============================================================================

with ntlink_log.open("w") as log_file:

    log_file.write("Starting ntLink Docker build\n")

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


# =============================================================================
# STEP 8 — READ DOCKER OUTPUT WHILE IT RUNS
# =============================================================================

    for line in process.stdout:
        print(line, end="")
        log_file.write(line)


# =============================================================================
# STEP 9 — WAIT FOR DOCKER TO FINISH
# =============================================================================

    return_code = process.wait()


# =============================================================================
# STEP 10 — CHECK WHETHER THE BUILD SUCCEEDED
# =============================================================================

if return_code != 0:
    print("ERROR: ntLink Docker image build failed.")
    print(f"Check the build log: {ntlink_log}")
    sys.exit(return_code)

print("SUCCESS: ntLink Docker image created.")


# =============================================================================
# STEP 11 — CHECK THAT NTLINK EXISTS INSIDE THE CONTAINER
# =============================================================================

runtime_command = [
    "docker",
    "run",
    "--rm",
    ntlink_image,
    "bash",
    "-lc",
    "command -v ntLink",
]

runtime_check = subprocess.run(
    runtime_command,
    capture_output=True,
    text=True,
)

if runtime_check.returncode != 0:
    print("ERROR: ntLink executable was not found inside the container.")
    print(runtime_check.stderr)
    sys.exit(runtime_check.returncode)

print("ntLink executable found:")
print(runtime_check.stdout.strip())

# =============================================================================
# STEP 12 — TEST THE NTLINK COMMAND
# =============================================================================

help_command=[
    "docker",
    "run",
    "--rm",
    ntlink_image,
    "ntLink",
    "help",
]

help_check= subprocess.run(
    help_command,
    capture_output=True,
    text=True,
)

if help_check.returncode != 0:
    print("ERROR: ntLink help command failed.")
    print(help_check.stderr)
    sys.exit(help_check.returncode)

# =============================================================================
# STEP 13 — FINAL VALIDATION SUMMARY
# =============================================================================

print()
print("=" * 60)
print("NTLINK CONTAINER VALIDATION PASSED")
print("=" * 60)

print()
print("Validation results:")
print("  Docker available            [PASS]")
print("  Docker image built          [PASS]")
print("  ntLink executable found     [PASS]")
print("  ntLink command executed     [PASS]")

print()
print("Docker image:")
print(f"  {ntlink_image}")

print()
print("Docker build context:")
print(f"  {ntlink_context}")

print()
print("Build log:")
print(f"  {ntlink_log}")

print()
print("No scaffolding job was executed.")
print()
