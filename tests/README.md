# Tests and validation

| Directory | Purpose |
|---|---|
| `fixtures/` | Small deterministic test inputs |
| `smoke/` | Fast end-to-end checks |
| `workflow_validation/` | Snakemake validation tests |
| `container_validation/` | Container availability and behavior checks |
| `validation/` | General validation records |

Large production inputs and outputs do not belong here. Alignment-specific fixtures may remain next to the relevant analysis pipeline under `alignment_analysis/scripts/30x/tests/fixtures/`.
