# Repository tree automation installation

Copy the four repository-relative paths from this kit into the repository root:

- `scripts/update_repository_tree.py`
- `scripts/install_git_hooks.sh`
- `.githooks/pre-commit`
- `.github/workflows/repository-tree.yml`

Then run, from the repository root:

```bash
python3 scripts/update_repository_tree.py
```

Review the README changes:

```bash
git diff -- README.md alignment_analysis/README.md assemblers/README.md assemblers/whole_genome_asm/README.md assemblers/whole_genome_asm/assessment/README.md assembly_analysis/README.md variant_calling_analysis/README.md
```

Install the versioned pre-commit hook:

```bash
bash scripts/install_git_hooks.sh
```

Validate that a second run makes no changes:

```bash
python3 scripts/update_repository_tree.py --check
```

The generator uses `git ls-files --cached`, so ignored local outputs, `.snakemake/`, untracked result folders, temporary files, and other local clutter are not inserted into public README trees. Newly staged files are included automatically before a commit.
