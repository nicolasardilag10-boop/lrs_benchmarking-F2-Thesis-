# Repository utilities

This directory contains project-maintenance utilities, not scientific analysis code.

| Script | Purpose |
|---|---|
| [`update_repository_tree.py`](update_repository_tree.py) | Refresh the compact repository explorers embedded in maintained READMEs |
| [`repository_cleanup_audit.py`](repository_cleanup_audit.py) | Report duplicates, empty files, and cleanup candidates without deleting scientific data |
| [`install_git_hooks.sh`](install_git_hooks.sh) | Install the local pre-commit hook that refreshes README navigation |
| [`test_variant_caller_docker_images.sh`](test_variant_caller_docker_images.sh) | Validate availability of variant-caller container images |

Scientific scripts belong inside their domain module, for example `alignment_analysis/scripts/`, `assembly_analysis/scripts/`, or `variant_calling_analysis/scripts/`.
