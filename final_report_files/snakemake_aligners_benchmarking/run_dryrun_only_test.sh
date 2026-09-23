#!/usr/bin/env bash
set -euo pipefail

echo "============================================================"
echo " DRY-RUN ONLY TEST: minimap2 + pbmm2 + VACmap + vg Giraffe"
echo "============================================================"

cd "$HOME/lrs_benchmarking/final_report_files/snakemake_aligners_benchmarking"

echo "[1/6] Creating dry-run dummy input files..."
mkdir -p .dryrun_inputs/samples .dryrun_inputs/reference .dryrun_inputs/graph envs

# Tiny valid dummy FASTQ files. They are only for DAG construction.
printf '@dryrun_read\nACGTACGT\n+\nIIIIIIII\n' | gzip -c > .dryrun_inputs/samples/DRYRUN.ont.1k.fastq.gz
printf '@dryrun_read\nACGTACGT\n+\nIIIIIIII\n' | gzip -c > .dryrun_inputs/samples/DRYRUN.pb.1k.fastq.gz

# Tiny dummy reference. It will not be used because this is dry-run only.
cat > .dryrun_inputs/reference/genome.fasta <<'EOF'
>chrDryRun
ACGTACGTACGTACGTACGTACGTACGTACGT
EOF

# Dummy vg graph/index placeholders. Dry-run only.
touch .dryrun_inputs/graph/graph.gbz
touch .dryrun_inputs/graph/graph.min
touch .dryrun_inputs/graph/graph.dist

echo "[2/6] Creating dry-run config with all four aligners active..."
cat > config.dryrun.yaml <<'YAML'
samples_dir: .dryrun_inputs/samples
reference: .dryrun_inputs/reference/genome.fasta

samples:
  - DRYRUN

technologies:
  - ont
  - pb

active_aligners:
  - minimap2
  - pbmm2
  - vacmap
  - vg_giraffe

vacmap_bin: "conda run -n vacmap_env vacmap"

vg_gbz: .dryrun_inputs/graph/graph.gbz
vg_min: .dryrun_inputs/graph/graph.min
vg_dist: .dryrun_inputs/graph/graph.dist
vg_zipcodes: ""
YAML

echo "[3/6] Creating minimal environment YAML files for dry-run..."
cat > envs/alignment.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - python=3.11
  - minimap2
  - pbmm2
  - samtools
  - pandas
  - numpy
  - scipy
  - matplotlib
  - pip
YAML

cat > envs/vg.yaml <<'YAML'
channels:
  - bioconda
  - conda-forge
  - nodefaults

dependencies:
  - vg
  - samtools
YAML

echo "[4/6] Creating dry-run Snakefile copy..."
python - <<'PY'
from pathlib import Path
import re

source = Path("snakemake_aligners.smk")
target = Path("snakemake_aligners.dryrun.smk")

if not source.exists():
    raise SystemExit("Missing snakemake_aligners.smk. Create the main Snakefile first.")

text = source.read_text(encoding="utf-8")

# Replace normal configfile with dry-run configfile.
if re.search(r'configfile:\s*["\']config\.yaml["\']', text):
    text = re.sub(
        r'configfile:\s*["\']config\.yaml["\']',
        'configfile: "config.dryrun.yaml"',
        text,
        count=1,
    )

# Fix older incorrect global conda block if it is still present.
elif 'conda:\n    "envs/alignment.yaml"' in text:
    text = text.replace(
        'conda:\n    "envs/alignment.yaml"',
        'configfile: "config.dryrun.yaml"',
        1,
    )

else:
    # Insert configfile after imports if no configfile is found.
    text = text.replace(
        "from pathlib import Path\n",
        'from pathlib import Path\n\nconfigfile: "config.dryrun.yaml"\n',
        1,
    )

target.write_text(text, encoding="utf-8")
print(f"Created: {target}")
PY

echo "[5/6] Listing rules..."
snakemake -s snakemake_aligners.dryrun.smk --list-rules

echo
echo "[6/6] Running DRY-RUN ONLY. No aligner will execute."
snakemake -s snakemake_aligners.dryrun.smk --sdm conda --cores 4 -np --forceall | tee dryrun_only.log

echo
echo "============================================================"
echo " DRY-RUN FINISHED"
echo "============================================================"
echo "Check this log:"
echo "  dryrun_only.log"
echo
echo "Expected: you should see planned jobs for:"
echo "  align_minimap2"
echo "  align_pbmm2"
echo "  align_vacmap"
echo "  align_vg_giraffe"
echo "  samtools_quality"
echo "  alignment_summary"
echo "  quality_check"
echo "  prepare_plot_data"
echo "  paired_t_tests"
echo "  alignment_figures"
echo
echo "No real samples were aligned."
echo "No BAM files were created."
