# Reproduction Checklist

## 1. Open the project in WSL

```bash
cd "$HOME/lrs_benchmarking"
conda activate lrs-plots
code .
```

VS Code should show:

```text
WSL: Ubuntu
```

## 2. Confirm the Python interpreter

```bash
which python
```

Expected:

```text
/home/nicolas/miniconda3/envs/lrs-plots/bin/python
```

Do not use:

```text
/usr/bin/python3
```

for the plotting environment.

## 3. Verify required packages

```bash
python -c "
import sys
import pandas
import numpy
import matplotlib
import scipy

print('Python:', sys.executable)
print('Pandas:', pandas.__version__)
print('NumPy:', numpy.__version__)
print('Matplotlib:', matplotlib.__version__)
print('SciPy:', scipy.__version__)
"
```

## 4. Syntax-check a script

```bash
python -m py_compile \
alignment_analysis/scripts/plot_general_technology_comparison.py
```

No output means the syntax check passed.

## 5. Run a plotting script

```bash
python \
alignment_analysis/scripts/plot_general_technology_comparison.py
```

## 6. Confirm outputs

```bash
ls -lh alignment_analysis/figures/
```

## 7. Validate before interpreting

Confirm:

- expected row count,
- no duplicate sample–technology–configuration rows,
- no missing paired values,
- correct category order,
- correct metric units,
- correct output paths.

## 8. Commit the work

```bash
git status --short
git add -A
git diff --cached --stat

git commit -m "Update long-read benchmarking documentation and figures"

git log -1 --oneline
git push origin "$(git branch --show-current)"
git status
```

Expected final state:

```text
nothing to commit, working tree clean
```

## 9. Preview the stylish documentation

```bash
python -m pip install -r requirements-docs.txt
mkdocs serve
```

Open the local address shown in the terminal.