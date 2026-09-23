# Structural-variant and aligner validation workspace

This workspace contains validation configurations and workflows supporting structural-variant caller development and aligner checks.

| Location | Purpose |
|---|---|
| [`workflow/aligners/`](workflow/aligners/) | minimap2, pbmm2, VACMap, and VG Giraffe validation rules |
| [`configs/aligner_validation.yaml`](configs/aligner_validation.yaml) | validation configuration |
| [`scripts/check_variant_caller_images.sh`](scripts/check_variant_caller_images.sh) | container-image validation |
| [`docs/`](docs/) | tool inventory and Docker validation notes |
| `results/` | generated validation outputs |
| `logs/` | generated logs |

Production root-level caller workflows remain governed by [`../CONSTITUTION.md`](../CONSTITUTION.md). Analysis-ready variant tables and plots belong in [`../variant_calling_analysis/`](../variant_calling_analysis/README.md).
