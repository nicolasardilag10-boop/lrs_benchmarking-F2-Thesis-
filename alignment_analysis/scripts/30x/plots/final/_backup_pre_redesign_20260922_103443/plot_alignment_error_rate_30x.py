from pathlib import Path
import sys


try:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import pandas as pd

except ModuleNotFoundError as error:
    raise SystemExit(
        "A required plotting package is missing: "
        f"{error.name}\n"
        f"Python interpreter: {sys.executable}\n"
    ) from error


# ============================================================
# File paths and plotting configuration
# ============================================================

PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "final/alignment_benchmark_30x.tsv")

OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/01_alignment_error_rate_30x.png")

OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "final/01_alignment_error_rate_30x.pdf")

METRIC = "error_percent"

ALIGNER_ORDER = [
    "minimap2",
    "pbmm2",
    "VACmap",
    "VG Giraffe",]

TECHNOLOGY_ORDER = [
    "ONT",
    "PacBio",]

SAMPLE_ORDER = [
    "HG002",
    "HG003",
    "HG004",]

SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",}

SAMPLE_OFFSETS = {
    "HG002": -0.22,
    "HG003": 0.0,
    "HG004": 0.22,}


# ============================================================
# Read, validate and select the benchmark data
# ============================================================

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"The input table was not found:\n{INPUT_FILE}")

data = pd.read_csv(
    INPUT_FILE,
    sep="\t",)

required_columns = {
    "sample",
    "read_technology",
    "configuration",
    METRIC,}

missing_columns = required_columns.difference(data.columns)

if missing_columns:
    raise ValueError(
        "The input table is missing these columns: "
        f"{sorted(missing_columns)}")

data[METRIC] = pd.to_numeric(
    data[METRIC],
    errors="raise",)

selected_conditions = (
    (
        (data["read_technology"] == "ONT")
        & data["configuration"].isin([
            "mm2-ont",
            "pbmm2-ont",
            "vacmap-ont",
            "vg-ont",
        ])
    )
    |
    (
        (data["read_technology"] == "PacBio")
        & data["configuration"].isin([
            "mm2-pb",
            "pbmm2-pb",
            "vacmap-pb",
            "vg-pb",
        ])
    )
)

plot_data = data.loc[
    selected_conditions,
    [
        "sample",
        "read_technology",
        "configuration",
        METRIC,
    ],
].copy()

configuration_to_aligner = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
    "vacmap-ont": "VACmap",
    "vacmap-pb": "VACmap",
    "vg-ont": "VG Giraffe",
    "vg-pb": "VG Giraffe",}

plot_data["aligner"] = (
    plot_data["configuration"]
    .map(configuration_to_aligner))

if plot_data["aligner"].isna().any():
    raise ValueError(
        "At least one selected configuration could not be mapped "
        "to an aligner.")

duplicated_rows = plot_data.duplicated(
    subset=[
        "sample",
        "read_technology",
        "aligner",
    ],
    keep=False,)

if duplicated_rows.any():
    raise ValueError(
        "Duplicated sample/technology/aligner observations found:\n"
        + plot_data.loc[duplicated_rows].to_string(index=False))


# Print the exact observations used in the figure.
plotting_table = (
    plot_data[
        [
            "read_technology",
            "aligner",
            "sample",
            METRIC,
        ]]
    .rename(
        columns={
            "read_technology": "technology",
            METRIC: "error_percent",})
    .sort_values(
        [
            "technology",
            "aligner",
            "sample",
        ]))

print("Python interpreter:")
print(sys.executable)
print()
print("General technology comparison plotting table:")
print(plotting_table.to_string(index=False))

expected_combinations = pd.MultiIndex.from_product(
    [
        TECHNOLOGY_ORDER,
        ALIGNER_ORDER,
        SAMPLE_ORDER,
    ],
    names=[
        "technology",
        "aligner",
        "sample",
    ])

observed_combinations = pd.MultiIndex.from_frame(
    plotting_table[
        [
            "technology",
            "aligner",
            "sample",
        ]])

missing_combinations = expected_combinations.difference(
    observed_combinations)

print()
print("Missing sample/aligner combinations:")
if len(missing_combinations) == 0:
    print("None")
else:
    for combination in missing_combinations:
        print("  - " + " | ".join(combination))


# ============================================================
# Create the boxplot general technology comparison
# ============================================================

maximum_error = plot_data[METRIC].max()
y_axis_maximum = maximum_error * 1.18

figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(13.0, 6.8),)

aligner_positions = list(
    range(len(ALIGNER_ORDER)))

for panel_index, technology in enumerate(TECHNOLOGY_ORDER):
    axis = axes[panel_index]

    boxplot_values = []
    boxplot_positions = []

    for aligner_index, aligner in enumerate(ALIGNER_ORDER):
        aligner_data = plot_data[
            (
                plot_data["read_technology"] == technology)
                & (plot_data["aligner"] == aligner)
        ]

        values = (
            aligner_data
            .set_index("sample")[METRIC]
            .reindex(SAMPLE_ORDER))

        available_values = values.dropna().to_numpy()

        if len(available_values) > 0:
            boxplot_values.append(available_values)
            boxplot_positions.append(aligner_index)

        for sample in SAMPLE_ORDER:
            value = values.loc[sample]
            x_position = (
                aligner_index
                + SAMPLE_OFFSETS[sample])

            if pd.isna(value):
                axis.annotate(
                    "NA",
                    xy=(
                        aligner_index,
                        y_axis_maximum * 0.025),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="#666666",)
                continue

            axis.scatter(
                x_position,
                value,
                color=SAMPLE_COLORS[sample],
                edgecolor="white",
                linewidth=0.7,
                s=58,
                label=(
                    sample
                    if aligner_index == 0
                    else "_nolegend_"),
                zorder=3,)

            axis.annotate(
                f"{value:.2f}%",
                xy=(
                    x_position,
                    value),
                xytext=(
                    0,
                    6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                zorder=4,)

    if boxplot_values:
        axis.boxplot(
            boxplot_values,
            positions=boxplot_positions,
            widths=0.58,
            patch_artist=False,
            showfliers=False,
            boxprops={
                "color": "black",
                "linewidth": 1.0,},
            whiskerprops={
                "color": "black",
                "linewidth": 1.0,},
            capprops={
                "color": "black",
                "linewidth": 1.0,},
            medianprops={
                "color": "black",
                "linewidth": 1.2,},
            zorder=2,)

    axis.set_title(
        "ONT" if technology == "ONT" else "PacBio HiFi",
        fontsize=16,
        pad=16,)
    axis.text(
        -0.08,
        1.03,
        "a" if technology == "ONT" else "b",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        fontweight="bold",)
    axis.set_ylim(
        0,
        y_axis_maximum,)
    axis.set_xticks(
        aligner_positions,
        ALIGNER_ORDER,)
    axis.grid(
        axis="y",
        color="#D9D9D9",
        linewidth=0.7,
        linestyle="-",
        zorder=0,)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("black")
    axis.spines["bottom"].set_color("black")
    axis.tick_params(
        axis="both",
        colors="black",
        labelsize=11,)
    axis.set_facecolor("white")

axes[0].set_ylabel(
    "Alignment error rate (%)",
    fontsize=13,)

legend_handles, legend_labels = axes[0].get_legend_handles_labels()
figure.legend(
    handles=legend_handles,
    labels=legend_labels,
    frameon=False,
    title="GIAB sample",
    title_fontsize=13,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.99),
    columnspacing=1.5,
    fontsize=12,)

figure.patch.set_facecolor("white")
figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.14,
    top=0.84,
    wspace=0.12,)

OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,)

figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",)

figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",)

plt.close(figure)

print()
print("Figure created successfully:")
print(OUTPUT_PNG)
print(OUTPUT_PDF)