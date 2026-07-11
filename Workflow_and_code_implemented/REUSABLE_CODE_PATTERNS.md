# Reusable Code Patterns

## 1. Define paths

```python
from pathlib import Path

PROJECT = Path.home() / "lrs_benchmarking"

INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "alignment_summary.tsv"
)
```

## 2. Read a TSV table

```python
import pandas as pd

data = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)
```

## 3. Convert a metric to numeric

```python
METRIC = "error_rate"

data[METRIC] = pd.to_numeric(
    data[METRIC],
    errors="raise",
)
```

## 4. Filter with Boolean conditions

```python
selected_conditions = (
    (
        (
            data["read_technology"] == "ONT"
        )
        & data["configuration"].isin([
            "mm2-ont",
            "pbmm2-ont",
        ])
    )
    |
    (
        (data["read_technology"] == "PacBio")
        & data["configuration"].isin([
            "mm2-pb",
            "pbmm2-pb",
        ])
    )
)
```

## 5. Select rows and columns

```python
plot_data = data.loc[
    selected_conditions,
    [
        "sample",
        "read_technology",
        "configuration",
        METRIC,
    ],
].copy()
```

## 6. Map configurations to aligners

```python
configuration_to_aligner = {
    "mm2-ont": "mm2",
    "mm2-pb": "mm2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
}

plot_data["aligner_group"] = (
    plot_data["configuration"]
    .map(configuration_to_aligner)
)
```

## 7. Detect duplicates

```python
duplicates = plot_data.duplicated(
    subset=[
        "sample",
        "read_technology",
        "configuration",
    ],
    keep=False,
)

if duplicates.any():
    raise ValueError(
        "Duplicated observations were found."
    )
```

## 8. Reshape paired data

```python
paired_data = plot_data.pivot(
    index="sample",
    columns="configuration",
    values=METRIC,
)
```

## 9. Check missing values

```python
if paired_data.isna().any().any():
    raise ValueError(
        "The paired table contains missing values."
    )
```

## 10. Convert columns to NumPy arrays

```python
first_values = (
    paired_data["mm2-ont"]
    .to_numpy(dtype=float)
)

second_values = (
    paired_data["pbmm2-ont"]
    .to_numpy(dtype=float)
)
```

## 11. Run a paired t-test

```python
from scipy.stats import ttest_rel

t_statistic, p_value = ttest_rel(
    first_values,
    second_values,
)
```

## 12. Create a two-panel figure

```python
import matplotlib.pyplot as plt

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(12.5, 6.5),
    sharey=True,
)
```

## 13. Draw paired points and lines

```python
axis.plot(
    x_positions,
    sample_values,
    linestyle=":",
    linewidth=1.3,
    color=sample_color,
    alpha=0.70,
)

axis.scatter(
    x_positions,
    sample_values,
    s=46,
    color=sample_color,
    edgecolor="black",
    linewidth=0.35,
)
```

## 14. Save high-quality outputs

```python
OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,
)

figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
)

plt.close(figure)
```

!!! note "Core pattern"
    Most final scripts follow the same sequence: define paths → read table → validate columns → filter → reshape → test → plot → save.