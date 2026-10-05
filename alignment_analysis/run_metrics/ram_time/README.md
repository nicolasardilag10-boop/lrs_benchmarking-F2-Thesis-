# Aligner RAM / Time Measurement (VG Giraffe, VACmap)

Companion to [`../../README.md`](../../README.md) and to the assembler version
[`assemblers/whole_genome_asm/ram_time/README.md`](../../../assemblers/whole_genome_asm/ram_time/README.md),
whose structure, sampler and parser this folder reuses. It only measures the
**computational cost** of the two graph/alternative aligners: wall-clock
time, peak RAM, CPU time and peak scratch disk use. It never produces CRAMs
used in the report.

```text
ont.read_mapping.vg.ram_time.smk
pb.read_mapping.vg.ram_time.smk
ont.read_mapping.vacmap.ram_time.smk
pb.read_mapping.vacmap.ram_time.smk
```

Shared with the assembler workflows (not copied):

```text
assemblers/whole_genome_asm/ram_time/sample_container_resources.sh   # background memory/disk sampler
assembly_analysis/scripts/metrics/parse_assembler_time_v.py          # time -v parser (--aligner vg_giraffe|vacmap)
```

## Why this folder exists

The production workflows in the repository root (`ont.read_mapping.vg.smk`,
`pb.read_mapping.vg.smk`, `ont.read_mapping.vacmap.smk`,
`pb.read_mapping.vacmap.smk`) only print `START`/`END` timestamps in their
`map_sort.log`. The tables next to this folder were computed from those
timestamps:

```text
../vg_giraffe_30x_runtime.tsv
../vacmap_30x_runtime.tsv
../vacmap_30x_wallclock_runtime.tsv
../vg_vacmap_30x_runtime.tsv
```

They give wall-clock time only. There is **no peak RAM, CPU time or disk
use** for VG Giraffe and VACmap, and the timestamps cannot provide them.

Each `*.ram_time.smk` re-runs the `vg_map_sort` / `vacmap_map_sort` rule of
its production counterpart with the **same Docker image, version, memory
ceilings and aligner parameters**:

| Workflow | Image | Mapping step | Mapping memory limit | Threads (production -> here) |
|---|---|---|---|---|
| VG ONT | `schimar/lrs-vg:v1.73.0` | `vg giraffe -b r10 -o BAM` | 160 GB | 16 -> 64 |
| VG PacBio HiFi | `schimar/lrs-vg:v1.73.0` | `vg giraffe -b hifi --output-format SAM` | 160 GB | 16 -> 64 |
| VACmap ONT | `schimar/lrs-vacmap:v1.2.0` | `vacmap -mode H` | 48 GB | 16 -> 64 |
| VACmap PacBio HiFi | `schimar/lrs-vacmap:v1.2.0` | `vacmap -mode L` | 48 GB | 16 -> 64 |

Sort (`samtools sort -@ 4`, 4 CPUs, 16 GB) and index (`samtools index`,
4 CPUs, 8 GB) are also identical to production. The only difference is the
mapping thread count: **64 threads** (`ram_time_threads`) instead of the
production 16, the same as the assembler ram_time runs, so all tools are
measured at 64 threads and are directly comparable. Because of this, the new
wall-clock times are expected to be **shorter** than the 16-thread
`*_30x_runtime.tsv` tables next to this folder; `cpu_hours` barely depends on
the thread count and stays comparable with them.

The memory limits are unchanged (160 GB VG, 48 GB VACmap). If a 64-thread
run is killed for memory (docker exit 137, or `container_peak_mem_gb`
close to the limit), that limit, not the aligner, is the bottleneck.

The local root workflows are the reference. The older copies on the upstream
branch `feature/new-mappers-bm` use `-m 64g` for VG and SAM output for VG
ONT. Both were later changed locally (OOM kills, malformed ONT SAM tags).

## How time and RAM are measured

Snakemake's `benchmark:` directive is not used: it samples the host process
tree and cannot see memory used inside a Docker container, so it would
report a near-zero RAM value. Instead, as for the assemblers:

- The production rule has **three stages**, each in its own container. Each
  stage is wrapped in `/usr/bin/time -v` *inside* its container
  (`--entrypoint /usr/bin/time`):

  | Stage | Command | Time file |
  |---|---|---|
  | `map` | `vg giraffe` or `vacmap` -> temporary BAM/SAM | `map_time_v.txt` |
  | `sort` | `samtools sort` -> CRAM | `sort_time_v.txt` |
  | `index` | `samtools index` -> CRAI | `index_time_v.txt` |

  The record step **sums** the elapsed and CPU time of the three stages and
  takes the **maximum** peak RSS across them. This is the same treatment as
  the GoldRush stages, and the sum matches the production `map_sort` rule's
  wall-clock time (`../vacmap_30x_wallclock_runtime.tsv`).
- Before the multi-hour mapping starts, the run step checks that the image
  contains `/usr/bin/time`, `python3` (needed by the record step) and a `du`
  supporting `-sb`, and that `docker stats` works on the host. It fails
  within seconds otherwise.

`time -v` reports the peak RAM of the **largest single process**, so every
stage also runs `sample_container_resources.sh` in the background. It is the
same sampler as for the assemblers, started once per stage on that stage's
named container `ramtime-{vg,vacmap}-<dataset>-{map,sort,index}`:

- Every **30 s** it reads the memory of the **whole container** (all
  processes together, page cache excluded) with `docker stats`.
- Every **5 min** it measures the size of the scratch directory with
  `du -sb` inside the running container (CONSTITUTION II.1).
- One extra disk measurement is taken right after `sort`, while the
  temporary BAM/SAM and the finished CRAM both exist. This is the disk peak,
  since production deletes the BAM/SAM only after sorting.

The per-stage sample files are merged into one `resource_samples.tsv`. The
run fails if it holds no memory sample at all. The `index` stage (and often
`sort`) is shorter than 30 s, so it may have no memory sample of its own;
its peak is still in `peak_rss_gb`.

## Rules: run -> record

| Rule | What it does | Kept? |
|---|---|---|
| `vg_run` / `vacmap_run` | Runs map, sort and index inside a scratch directory declared as `temp(directory(...))` | No |
| `vg_record` / `vacmap_record` | Parses the three time files and the resource samples into `{dataset}.ram_time.tsv` with `parse_assembler_time_v.py --aligner ...`, run with `python3` inside the aligner's own image (CONSTITUTION II.1), and copies the raw files to `{dataset}.raw/` | **Yes** |

- The scratch directory (temporary BAM/SAM, CRAM, CRAI, time files) is
  deleted by Snakemake as soon as the record step succeeds, and also when
  the run fails. Do not pass `--notemp`.
- If the record step fails, the scratch directory is kept so the time files
  can be re-parsed without re-running the mapping.
- The run step refuses to start if Snakemake grants fewer threads than
  `ram_time_threads`, so always pass `--cores 64` (or more).
- Nothing is ever written to `cram/`. `idxstats` and `stats` are not re-run:
  they are QC on the CRAM, not part of the mapping cost.
- VG indexes are **not** rebuilt or timed. The four files in `vg_index/`
  (`hg38.giraffe.gbz`, `hg38.dist`, `hg38.longread.withzip.min`,
  `hg38.longread.zipcodes`) must already exist, otherwise Snakemake stops
  with a missing-input error. Use the production `build_vg_index` rule
  first if needed.

## Configuration

No paths are hard-coded. Everything is resolved relative to the directory
you launch Snakemake from (the repository root). Optional `--config` keys:

| Key | Default | Purpose |
|---|---|---|
| `ram_time_threads` | `64` | Threads given to the aligner and to its `docker --cpus` (sort/index keep 4) |
| `ram_time_dir` | `aligners_ram_time` | Where `.ram_time.tsv` files and rule logs are written |
| `ram_time_scratch_dir` | `aligners_ram_time/scratch` | Temporary mapping workspace, deleted after each run |
| `reference` | `reference/GRCh38_GIABv3_..._KCNJ18.fasta` | Same key and default as production |
| `dataset_filter` | none | Only datasets whose name contains this string, as in production |

Both directories must stay **inside** the repository root
(CONSTITUTION I.1); the workflow refuses to start otherwise.

Example, one sample only:

```bash
snakemake --snakefile alignment_analysis/run_metrics/ram_time/ont.read_mapping.vg.ram_time.smk --cores 64 --config dataset_filter=HG002
```

## Inputs and outputs

Inputs are the production FASTQs, discovered with the production filter
(`.ont.` or `.pb.` in the name), minus the test subsets skipped by the
assembler workflows (`.1k`, `.chr21.`, `localtest`, `smoke`):

```text
fastq/{sample}.{ont,pb}.30x.fastq.gz
reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
vg_index/hg38.*                     # VG only
```

Outputs: one small TSV per dataset, the raw files it was computed from, and
the logs:

```text
aligners_ram_time/{vg_giraffe,vacmap}/{dataset}.ram_time.tsv
aligners_ram_time/{vg_giraffe,vacmap}/{dataset}.raw/   # 3 time -v files + resource_samples*.tsv
aligners_ram_time/{vg_giraffe,vacmap}/{dataset}.{run,record}.log
```

Columns are the same as for the assemblers, except that the first one is
`aligner` (`vg_giraffe` or `vacmap`) instead of `assembler`:
`aligner`, `sample`, `technology`, `threads`, `wall_clock_seconds`,
`wall_clock_hours`, `peak_rss_gb`, `cpu_hours`, `cpu_efficiency`,
`exit_status`, `container_peak_mem_gb`, `peak_scratch_disk_gb`,
`n_mem_samples`, `n_stages_summed` (always 3), `source_time_files`.

- `wall_clock_*`: map + sort + index, at 64 mapping threads. The
  16-thread equivalent is `../vacmap_30x_wallclock_runtime.tsv`.
- `peak_rss_gb`: peak RAM of the largest single process (`time -v`), max
  over the three stages. In practice this is the mapping stage.
- `cpu_hours`: user + system CPU time, summed over the three stages.
- `cpu_efficiency`: CPU time / (wall-clock time x `threads`). Sort and index
  use only 4 threads, so they pull it slightly below the mapping-only value.
- `exit_status`: always 0. If any time file reports a non-zero exit status,
  the record step fails and no row is written.
- `container_peak_mem_gb`: highest sampled memory of the whole container.
- `peak_scratch_disk_gb`: largest scratch size, normally the
  BAM/SAM + CRAM measurement right after sorting. VACmap writes uncompressed
  SAM, so expect it to be much larger than for VG ONT (BAM).
- `n_mem_samples`: number of 30 s memory samples (about 120 per hour).

Mapping-only time (the 64-thread counterpart of `../vg_giraffe_30x_runtime.tsv`
and `../vacmap_30x_runtime.tsv`) is the `Elapsed (wall clock)` line of
`{dataset}.raw/map_time_v.txt`, and is also printed per stage in
`{dataset}.record.log`:

```bash
grep -H 'Elapsed (wall clock)' aligners_ram_time/*/*.raw/map_time_v.txt
```

Which peak RAM to report:

| Tool | Use | Why |
|---|---|---|
| VG Giraffe | `peak_rss_gb` | One multi-threaded `vg` process; `time -v` sees every byte |
| VACmap | `container_peak_mem_gb` | Python tool that runs worker processes in parallel; `peak_rss_gb` sees only the largest one and can undercount |

For VG the two values should be close; a large gap is worth a look in
`{dataset}.raw/resource_samples.tsv`. Sampled values can miss spikes
shorter than 30 s, so they are slightly below the true peak. The `/tmp`
tmpfs (50 GB, as in production) is memory-backed. Files `samtools sort`
writes there can count towards container memory.

## Where to run

Launch from the **repository root** on the host that has `fastq/`,
`reference/`, `vg_index/` and the Docker daemon used for the production
mapping runs.

Check before running:

```bash
cd /path/to/lrs_benchmarking    # replace with the real path on this host
ls -lh fastq/*.30x.fastq.gz vg_index/hg38.*
docker info >/dev/null && echo "Docker daemon reachable"
docker stats --no-stream >/dev/null && echo "docker stats works (needed by the sampler)"
```

Check that both images contain what the run step needs. This is the same
check the workflow does, run by hand:

```bash
for img in schimar/lrs-vg:v1.73.0 schimar/lrs-vacmap:v1.2.0; do docker run --rm --entrypoint sh "$img" -c 'command -v /usr/bin/time && command -v python3 && du -sb /tmp >/dev/null' && echo "$img OK" || echo "$img MISSING a tool"; done
```

### Verify which FASTQs each workflow will use

Every `fastq/*.fastq.gz` whose name contains `.ont.` (or `.pb.`) is used, not
only `30x`, so any extra file there becomes an extra multi-hour job. List
each input with the real file it points to and its size:

```bash
for f in fastq/*.fastq.gz; do case "$f" in *.1k*|*.chr21.*|*localtest*|*smoke*|*SMOKE*) continue;; esac; printf '%-35s %6s  %s\n' "$f" "$(du -hL "$f" | cut -f1)" "$(readlink -f "$f")"; done
```

Then do a dry run of each workflow. At startup, each one prints the
`.ram_time.tsv` targets it found, followed by the job count:

```bash
for smk in alignment_analysis/run_metrics/ram_time/*.ram_time.smk; do echo "=== $smk"; snakemake --snakefile "$smk" --cores 64 --dry-run --quiet rules 2>&1 | grep -E 'ram_time.tsv|_run|_record|Error'; done
```

Expected output, one target per dataset (shown for HG002):

```text
=== .../ont.read_mapping.vacmap.ram_time.smk
['.../aligners_ram_time/vacmap/HG002.ont.30x.ram_time.tsv']
vacmap_record        1
vacmap_run           1
...
```

An empty list (`[]`) means no input was found: you are not in the
repository root, or `fastq/` is missing.

## Commands

Always do a dry run first (add `--dry-run`), then the real run.

### VG Giraffe 1.73.0 (`schimar/lrs-vg:v1.73.0`)

```bash
# ONT
snakemake --snakefile alignment_analysis/run_metrics/ram_time/ont.read_mapping.vg.ram_time.smk --cores 64 --resources mem_mb=163840 --rerun-incomplete --printshellcmds --show-failed-logs

# PacBio HiFi
snakemake --snakefile alignment_analysis/run_metrics/ram_time/pb.read_mapping.vg.ram_time.smk --cores 64 --resources mem_mb=163840 --rerun-incomplete --printshellcmds --show-failed-logs
```

### VACmap 1.2.0 image (`schimar/lrs-vacmap:v1.2.0`)

```bash
# ONT
snakemake --snakefile alignment_analysis/run_metrics/ram_time/ont.read_mapping.vacmap.ram_time.smk --cores 64 --rerun-incomplete --printshellcmds --show-failed-logs

# PacBio HiFi
snakemake --snakefile alignment_analysis/run_metrics/ram_time/pb.read_mapping.vacmap.ram_time.smk --cores 64 --rerun-incomplete --printshellcmds --show-failed-logs
```

`--cores 64` runs one dataset at a time, so no two mappings compete for
CPU or memory and the times stay comparable.

### If a run is interrupted

After an interruption (e.g. Ctrl-C), check that no container is left over.
A leftover one makes the next run fail with "name is already in use":

```bash
docker ps -a --filter name=ramtime- --format '{{.Names}} {{.Status}}'
docker rm -f ramtime-<tool>-<dataset>-<stage>    # only if listed and not wanted
```

## Combining the results

After all runs finish (one `.ram_time.tsv` per dataset across both aligners
and technologies), merge them into one table:

```bash
{
    head -n 1 "$(ls aligners_ram_time/*/*.ram_time.tsv | head -n 1)"
    tail -n +2 -q aligners_ram_time/*/*.ram_time.tsv
} > alignment_analysis/run_metrics/vg_vacmap_30x_ram_time.tsv
```

Thread-hours are `threads x wall_clock_hours`, both already in the table.

## Validation before commit

From the repository root, for every modified workflow:

```bash
git diff --check
snakemake --snakefile alignment_analysis/run_metrics/ram_time/<workflow>.ram_time.smk --cores 64 --dry-run --printshellcmds
```
