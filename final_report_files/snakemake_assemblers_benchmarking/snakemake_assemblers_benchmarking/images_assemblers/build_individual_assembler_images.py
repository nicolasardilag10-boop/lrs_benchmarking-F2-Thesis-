#!/usr/bin/env python3

# ============================================================
# GENERAL WORKFLOW
# ============================================================
# 1. Import the Python modules we need.
# 2. Define where the project, Dockerfiles, and logs are located.
# 3. Check that Docker is running before trying to build anything.
# 4. Define the information needed to build the Flye2 Docker image.
# 5. Construct the equivalent of the "docker build" Bash command.
# 6. Start Docker from Python using subprocess.Popen().
# 7. Read Docker output line by line while the build is running.
# 8. Show that output in the terminal AND save it in a log file.
# 9. Wait until Docker finishes.
# 10. Check Docker's return code:
#       0     = build succeeded
#       != 0  = build failed
# ============================================================


# ============================================================
# STEP 1 — IMPORT PYTHON MODULES
# ============================================================

# allow to use directories
# Path lets us construct file and directory paths in a safer way
# than manually writing long strings such as "/home/nicolas/..."
from pathlib import Path

# subprocess allows Python to execute external programs,
# such as Docker commands
import subprocess

# sys is used here mainly to stop the Python script with sys.exit()
import sys


# ============================================================
# STEP 2 — DEFINE PROJECT DIRECTORIES
# ============================================================

# PROJECT is the main assembler benchmarking project directory.
#
# Path.home() automatically gives:
# /home/nicolas
#
# Then "/" is used by pathlib to join directories.
PROJECT = (
    Path.home()
    / "lrs_benchmarking"
    / "final_report_files"
    / "snakemake_assemblers_benchmarking"
)

# Directory containing the individual Docker build contexts.
#
# For example:
# containers/flye2/
# containers/goldrush/
# containers/ntlink/
# containers/verkko2/
CONTAINERS = PROJECT / "containers"

# Directory where Docker build logs will be stored.
LOGDIR = PROJECT / "container_build_logs"

# Create the log directory if it does not already exist.
#
# parents=True
#     also creates missing parent directories if necessary.
#
# exist_ok=True
#     does not produce an error if the directory already exists.
LOGDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# STEP 3 — CHECK THAT DOCKER IS AVAILABLE
# ============================================================

print("Checking Docker...")

# Run:
#
# docker info
#
# We only care about whether Docker works or not,
# so we hide stdout and stderr using DEVNULL.
docker_check = subprocess.run(
    ["docker", "info"],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)

# If docker succed return 0; if not 1+ error and stop the run
#
# Docker/Linux programs normally use:
#
# return code 0     -> success
# return code != 0  -> something failed
if docker_check.returncode != 0:
    print("ERROR: Docker is not available")
    print("Start Docker Desktop and enable WSL integration.")
    sys.exit(1)

# Python only reaches this line if Docker returned 0.
print("Docker is ready.")


# ============================================================
# STEP 4 — START THE FLYE2 BUILD SECTION
# ============================================================
#take the string and repeat it 60 times so it will create separations with "==="

print()
print("=" * 60)
print("BUILDING FLYE 2")
print("=" * 60)


#####################
# FLYE2
#####################

# ============================================================
# STEP 5 — DEFINE VARIABLES REGARDING FLYE2
# ============================================================

# variables regarding flye 2

# Docker build context.
#
# Docker will look inside:
# containers/flye2/
#
# for its Dockerfile and any files needed during the build.
flye_context = CONTAINERS / "flye2"

# Name and tag that the new local Docker image will receive.
#
# Result:
# lrs-flye2:2.9.6
flye_image = "lrs-flye2:2.9.6"

# File where the complete Docker build output will be saved.
flye_log = LOGDIR / "flye2.build.log"


# ============================================================
# STEP 6 — CREATE THE DOCKER BUILD COMMAND
# ============================================================

# using the same bash command but for python
#
# This Python list represents approximately this Bash command:
#
# docker build \
#     --progress=plain \
#     --pull \
#     --tag lrs-flye2:2.9.6 \
#     /path/to/containers/flye2
#
# Each element of the list represents one command-line argument.
command = [
    "docker",
    "build",
    "--progress=plain",
    "--pull",
    "--tag",
    flye_image,
    str(flye_context),
]

# Show the command before executing it.
#
# " ".join(command) converts the Python list into a readable
# shell-like command for us to inspect.
print("Command:")
print(" ".join(command))


# ============================================================
# STEP 7 — OPEN THE LOG FILE
# ============================================================

# Save each line
#
# "w" means write mode.
# A new build replaces the previous build log.
#
# "with" automatically closes the file when this block finishes.
with flye_log.open("w") as log_file:

    # Write an initial message into the log before Docker starts.
    log_file.write("Starting Flye build\n")


    # ========================================================
    # STEP 8 — START THE DOCKER BUILD
    # ========================================================

    # Popen starts Docker but allows Python to continue interacting
    # with the process while Docker is still running.
    #
    # This is different from subprocess.run(), which normally waits
    # until the command completely finishes before continuing.
    process = subprocess.Popen(
        command,

        # Allow python to capture the output and can be read
        #
        # PIPE creates a connection between Docker's stdout
        # and our Python program.
        stdout=subprocess.PIPE,

        # Capture the ouput
        #
        # Redirect Docker's stderr into stdout.
        # Therefore warnings, progress information and errors
        # can all be handled through process.stdout.
        stderr=subprocess.STDOUT,

        # Convert Docker output into normal Python strings
        # instead of raw bytes.
        text=True,
    )


    # ========================================================
    # STEP 9 — READ DOCKER OUTPUT WHILE IT RUNS
    # ========================================================

    # process.stdout contains the output arriving from Docker.
    #
    # Every time Docker produces another line:
    #     1. Python stores it temporarily in "line"
    #     2. prints it in the terminal
    #     3. saves it in flye2.build.log
    for line in process.stdout:

        # Show Docker progress live in the terminal.
        #
        # end="" avoids adding an extra newline because Docker's
        # line normally already contains "\n".
        print(line, end="")

        # Save exactly the same Docker line in the build log.
        log_file.write(line)


    # ========================================================
    # STEP 10 — WAIT FOR DOCKER TO FINISH
    # ========================================================

    # Once Docker stops producing output, wait() makes sure the
    # process has completely finished.
    #
    # It returns Docker's exit code:
    #
    # 0     = successful build
    # != 0  = failed build
    return_code = process.wait()


    # ========================================================
    # STEP 11 — CHECK WHETHER THE BUILD SUCCEEDED
    # ========================================================

    if return_code != 0:
        print("ERROR: Flye Docker image build failed")

        # Stop the Python script and preserve Docker's error code.
        sys.exit(return_code)

    # Python reaches here only when return_code == 0.
    print("SUCCESS: Flye Docker image created")