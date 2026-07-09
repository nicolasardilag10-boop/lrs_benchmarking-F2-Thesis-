from pathlib import Path
import sys


# ------------------------------------------------------------
# STEP 1: Import the required Python packages
# ------------------------------------------------------------

# The try/except block provides a clearer message when
# VS Code runs the script with the wrong Python interpreter.
try:
    # Matplotlib creates and saves the graph.
    import matplotlib

    # "Agg" is a non-interactive backend.
    # It allows the script to save figures without needing
    # to open a separate graphical window.
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    # NumPy handles numerical arrays and calculations.
    import numpy as np

    # Pandas loads, filters, and reorganizes the TSV table.
    import pandas as pd

    # ttest_rel performs a paired t-test.
    #
    # It is appropriate here because every original input
    # has one mm2 result and one corresponding pbmm2 result.
    from scipy.stats import ttest_rel

except ModuleNotFoundError as error:
    raise SystemExit(
        "\nA required Python package is missing.\n"
        f"Missing package: {error.name}\n"
        f"Python interpreter used: {sys.executable}\n\n"
        "VS Code must use this interpreter:\n"
        "/home/nicolas/miniconda3/envs/lrs-plots/bin/python\n"
    ) from error


# Print the interpreter used by the VS Code Run button.
#
# The path should contain:
# miniconda3/envs/lrs-plots
print("Python interpreter:")
print(sys.executable)


# ------------------------------------------------------------
# STEP 2: Define the project and file paths
# ------------------------------------------------------------

# Path.home() returns:
# /home/nicolas
#
# The / operator joins directory names safely.
PROJECT = Path.home() / "lrs_benchmarking"


# Input table containing the alignment statistics.
INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "alignment_summary.tsv"
)


# Output PNG figure.
OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "general_technology_comparison.png"
)


# Output PDF figure.
OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "general_technology_comparison.pdf"
)


# Stop if the input table does not exist.
if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"The input table was not found:\n{INPUT_FILE}"
    )


print()
print("Input file:")
print(INPUT_FILE)


# ------------------------------------------------------------
# STEP 3: Read the TSV table
# ------------------------------------------------------------

# sep="\t" means that the columns are separated by tabs.
data = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)


print()
print("Original table dimensions:")
print(data.shape)


# ------------------------------------------------------------
# STEP 4: Select the metric
# ------------------------------------------------------------

# error_rate uses decimal values.
#
# Example:
# 0.03 = 3% error
METRIC = "error_rate"


# These columns are required for the analysis.
required_columns = {
    "sample",
    "read_technology",
    "configuration",
    METRIC,
}


# Determine whether any required columns are missing.
missing_columns = required_columns.difference(
    data.columns
)


if missing_columns:
    raise ValueError(
        "The input table is missing these columns: "
        f"{sorted(missing_columns)}"
    )


# Convert the selected metric to numeric values.
#
# errors="raise" stops the script if a value cannot
# be converted into a number.
data[METRIC] = pd.to_numeric(
    data[METRIC],
    errors="raise",
)


# ------------------------------------------------------------
# STEP 5: Select technology-appropriate configurations
# ------------------------------------------------------------

# ONT comparison:
# minimap2 map-ont versus pbmm2 SUBREAD
#
# PacBio comparison:
# minimap2 map-hifi versus pbmm2 CCS/HIFI
#
# & means AND
# | means OR
selected_conditions = (
    (
        (data["read_technology"] == "ONT")
        & (
            data["configuration"].isin([
                "mm2-ont",
                "pbmm2-ont",
            ])
        )
    )
    |
    (
        (data["read_technology"] == "PacBio")
        & (
            data["configuration"].isin([
                "mm2-pb",
                "pbmm2-pb",
            ])
        )
    )
)


# data.loc selects:
#
# 1. Rows where selected_conditions is True.
# 2. Only the listed columns.
#
# .copy() creates an independent DataFrame.
plot_data = data.loc[
    selected_conditions,
    [
        "sample",
        "read_technology",
        "configuration",
        METRIC,
    ],
].copy()


# Expected:
#
# 3 samples
# × 2 technologies
# × 2 aligners
# = 12 rows
if len(plot_data) != 12:
    raise ValueError(
        f"Expected 12 selected rows, "
        f"but found {len(plot_data)}."
    )


print()
print("Selected plotting data:")

print(
    plot_data.sort_values(
        [
            "sample",
            "read_technology",
            "configuration",
        ]
    ).to_string(
        index=False,
    )
)


# ------------------------------------------------------------
# STEP 6: Convert four configurations into two aligner groups
# ------------------------------------------------------------

# The table contains four configuration names.
#
# This dictionary reduces them to two general aligner groups:
# mm2 and pbmm2.
configuration_to_aligner = {
    "mm2-ont": "mm2",
    "mm2-pb": "mm2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
}


# .map() replaces each configuration according to
# the dictionary above.
plot_data["aligner_group"] = (
    plot_data["configuration"]
    .map(configuration_to_aligner)
)


# Check that every configuration was mapped.
if plot_data["aligner_group"].isna().any():
    raise ValueError(
        "At least one configuration could not be mapped "
        "to an aligner group."
    )


# ------------------------------------------------------------
# STEP 7: Create original input-file identifiers
# ------------------------------------------------------------

# Short suffixes used in the legend.
technology_suffix = {
    "ONT": "ont",
    "PacBio": "pb",
}


# Create labels such as:
#
# HG002_ont
# HG002_pb
# HG003_ont
# HG003_pb
plot_data["orig_file"] = (
    plot_data["sample"]
    + "_"
    + plot_data["read_technology"].map(
        technology_suffix
    )
)


if plot_data["orig_file"].isna().any():
    raise ValueError(
        "At least one original input label "
        "could not be created."
    )


# ------------------------------------------------------------
# STEP 8: Check for duplicate measurements
# ------------------------------------------------------------

# Each combination of:
#
# original input + aligner
#
# must occur exactly once.
duplicated_pairs = plot_data.duplicated(
    subset=[
        "orig_file",
        "aligner_group",
    ],
    keep=False,
)


if duplicated_pairs.any():
    duplicated_data = plot_data.loc[
        duplicated_pairs
    ]

    raise ValueError(
        "Duplicated input-aligner observations were found:\n"
        + duplicated_data.to_string(
            index=False
        )
    )


# ------------------------------------------------------------
# STEP 9: Validate the paired structure
# ------------------------------------------------------------

# Count how many different aligners exist for every input.
#
# Every input should have:
# one mm2 result
# one pbmm2 result
#
# Therefore, each count should equal 2.
pair_counts = (
    plot_data.groupby(
        "orig_file"
    )["aligner_group"]
    .nunique()
)


print()
print("Number of aligners per input:")
print(pair_counts)


# .all() returns True only when every count equals 2.
if not (pair_counts == 2).all():
    raise ValueError(
        "Every input file must contain exactly "
        "one mm2 result and one pbmm2 result."
    )


# ------------------------------------------------------------
# STEP 10: Reshape the table into paired columns
# ------------------------------------------------------------

# Before pivot:
#
# orig_file    aligner_group    error_rate
# HG002_ont    mm2              value
# HG002_ont    pbmm2            value
#
# After pivot:
#
# orig_file        mm2        pbmm2
# HG002_ont        value      value
#
# This structure guarantees that corresponding observations
# are placed on the same row.
paired_data = (
    plot_data.pivot(
        index="orig_file",
        columns="aligner_group",
        values=METRIC,
    )
    .reindex(
        columns=[
            "mm2",
            "pbmm2",
        ]
    )
    .sort_index()
)


print()
print("Paired mm2 and pbmm2 values:")
print(paired_data)


# Check for missing values.
if paired_data.isna().any().any():
    raise ValueError(
        "The paired table contains a missing "
        "mm2 or pbmm2 value."
    )


# Six paired inputs are expected:
#
# HG002_ont
# HG002_pb
# HG003_ont
# HG003_pb
# HG004_ont
# HG004_pb
if len(paired_data) != 6:
    raise ValueError(
        f"Expected 6 paired input files, "
        f"but found {len(paired_data)}."
    )


# ------------------------------------------------------------
# STEP 11: Convert the paired columns into NumPy arrays
# ------------------------------------------------------------

# Convert the mm2 Pandas column into a NumPy array.
mm2_values = (
    paired_data["mm2"]
    .to_numpy(
        dtype=float,
    )
)


# Convert the pbmm2 Pandas column into a NumPy array.
pbmm2_values = (
    paired_data["pbmm2"]
    .to_numpy(
        dtype=float,
    )
)


print()
print("mm2 values:")
print(mm2_values)

print()
print("pbmm2 values:")
print(pbmm2_values)


# ------------------------------------------------------------
# STEP 12: Calculate the means
# ------------------------------------------------------------

# Calculate one arithmetic mean for each aligner.
mm2_mean = np.mean(
    mm2_values
)

pbmm2_mean = np.mean(
    pbmm2_values
)


mean_values = np.array([
    mm2_mean,
    pbmm2_mean,
])


print()
print("Mean mm2 error rate:")
print(f"{mm2_mean:.6f}")

print()
print("Mean pbmm2 error rate:")
print(f"{pbmm2_mean:.6f}")


# ------------------------------------------------------------
# STEP 13: Perform the paired t-test
# ------------------------------------------------------------

# ttest_rel compares the paired values.
#
# Each position in mm2_values corresponds to the same
# input-file position in pbmm2_values.
t_statistic, p_value = ttest_rel(
    mm2_values,
    pbmm2_values,
)


# Convert the p-value into significance stars.
#
# p < 0.001  -> ***
# p < 0.01   -> **
# p < 0.05   -> *
# otherwise  -> ns
if p_value < 0.001:
    significance_label = "***"
elif p_value < 0.01:
    significance_label = "**"
elif p_value < 0.05:
    significance_label = "*"
else:
    significance_label = "ns"


# Create the text displayed above the boxes.
if p_value < 0.0001:
    p_value_text = (
        f"Paired t-test, p = {p_value:.2e} "
        f"({significance_label})"
    )
else:
    p_value_text = (
        f"Paired t-test, p = {p_value:.4f} "
        f"({significance_label})"
    )


print()
print("Paired t-test:")
print(f"t-statistic = {t_statistic:.6f}")
print(f"p-value = {p_value:.6f}")
print(f"Significance = {significance_label}")


# ------------------------------------------------------------
# STEP 14: Prepare plotting values and positions
# ------------------------------------------------------------

# Matplotlib expects one array for every box.
boxplot_values = [
    mm2_values,
    pbmm2_values,
]


# Horizontal locations of the two boxes.
x_positions = np.array([
    1.0,
    2.0,
])


# Combine all observations for automatic axis positioning.
all_values = np.concatenate([
    mm2_values,
    pbmm2_values,
])


minimum_value = all_values.min()
maximum_value = all_values.max()
data_range = maximum_value - minimum_value


# Prevent positioning errors if all values are identical.
if data_range == 0:
    data_range = 1.0


# ------------------------------------------------------------
# STEP 15: Create the figure
# ------------------------------------------------------------

figure, axis = plt.subplots(
    figsize=(
        8.5,
        5.8,
    )
)


# ------------------------------------------------------------
# STEP 16: Draw white box plots
# ------------------------------------------------------------

boxplot = axis.boxplot(
    boxplot_values,

    # Place mm2 at x=1 and pbmm2 at x=2.
    positions=x_positions,

    # Width of each box.
    widths=0.52,

    # Allow control of the box face color.
    patch_artist=True,

    # Do not add automatic outlier markers because every
    # observation will be plotted manually.
    showfliers=False,

    # White box interiors with black borders.
    boxprops={
        "facecolor": "white",
        "edgecolor": "black",
        "linewidth": 1.2,
    },

    # Black median lines.
    medianprops={
        "color": "black",
        "linewidth": 1.6,
    },

    # Black whiskers.
    whiskerprops={
        "color": "black",
        "linewidth": 1.2,
    },

    # Black whisker caps.
    capprops={
        "color": "black",
        "linewidth": 1.2,
    },
)


# ------------------------------------------------------------
# STEP 17: Add paired points and dotted lines
# ------------------------------------------------------------

# tab10 provides distinct colors for the six inputs.
color_map = plt.get_cmap(
    "tab10"
)


# Extract input labels in the same order as paired_data.
input_names = paired_data.index.tolist()


# Draw one paired line for each original input.
for index, input_name in enumerate(
    input_names
):
    # Extract the two corresponding values:
    # mm2 and pbmm2.
    paired_values = (
        paired_data.loc[
            input_name,
            [
                "mm2",
                "pbmm2",
            ],
        ]
        .to_numpy(
            dtype=float,
        )
    )

    # Choose one color for the complete pair.
    input_color = color_map(
        index
    )

    # Draw the dotted connecting line.
    axis.plot(
        x_positions,
        paired_values,
        linestyle=":",
        linewidth=1.3,
        color=input_color,
        alpha=0.75,
        zorder=2,
    )

    # Draw the paired points.
    axis.scatter(
        x_positions,
        paired_values,
        s=42,
        color=input_color,
        edgecolor="black",
        linewidth=0.35,
        label=input_name,
        zorder=3,
    )


# ------------------------------------------------------------
# STEP 18: Display the mean inside each box
# ------------------------------------------------------------

# Place the numerical mean at its real vertical position.
#
# The text has a small white background so that it remains
# readable if a line passes behind it.
for x_position, mean_value in zip(
    x_positions,
    mean_values,
):
    label_y = mean_value + data_range * 0.01

    axis.text(
        x_position,
        label_y,
        f"{mean_value:.3f}",
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=9,
        color="black",
        bbox={
            "facecolor": "none",
            "edgecolor": "none",
            "alpha": 0.90,
            "pad": 1.5,
        },
        zorder=5,
    )


# ------------------------------------------------------------
# STEP 19: Add the paired t-test bracket
# ------------------------------------------------------------

# Place the bracket above the highest observation.
bracket_y = (
    maximum_value
    + data_range * 0.16
)


# Height of the small vertical ends of the bracket.
bracket_height = (
    data_range * 0.035
)


# Draw the bracket:
#
# └────────────────┘
axis.plot(
    [
        1.0,
        1.0,
        2.0,
        2.0,
    ],
    [
        bracket_y,
        bracket_y + bracket_height,
        bracket_y + bracket_height,
        bracket_y,
    ],
    color="black",
    linewidth=1.2,
    clip_on=False,
)


# Place the p-value above the bracket.
axis.text(
    1.5,
    bracket_y + bracket_height + data_range * 0.025,
    p_value_text,
    horizontalalignment="center",
    verticalalignment="bottom",
    fontsize=10,
    color="black",
)


# ------------------------------------------------------------
# STEP 20: Format the axes and title
# ------------------------------------------------------------

# Define x-axis positions and labels.
axis.set_xticks(
    x_positions
)

axis.set_xticklabels([
    "mm2",
    "pbmm2",
])


# Axis labels.
axis.set_xlabel(
    "Aligner"
)

axis.set_ylabel(
    "Alignment error rate"
)


# Figure title.
axis.set_title(
    "Alignment error rate: minimap2 versus pbmm2"
)


# Remove all gridlines.
axis.grid(
    False
)


# Keep only the left and bottom axis borders.
axis.spines["top"].set_visible(
    False
)

axis.spines["right"].set_visible(
    False
)


# Define horizontal limits around the boxes.
axis.set_xlim(
    0.55,
    2.45,
)


# Increase the upper limit so the p-value annotation
# is fully visible.
axis.set_ylim(
    minimum_value - data_range * 0.12,
    maximum_value + data_range * 0.34,
)


# ------------------------------------------------------------
# STEP 21: Add the legend
# ------------------------------------------------------------

axis.legend(
    title="Input file",

    # Place the legend to the right of the plotting area.
    bbox_to_anchor=(
        1.06,
        0.5,
    ),

    loc="center left",

    # Remove the legend border.
    frameon=False,
)


# ------------------------------------------------------------
# STEP 22: Save the figure
# ------------------------------------------------------------

# Create the output directory if it does not exist.
OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# Improve spacing between figure elements.
figure.tight_layout()


# Save a high-resolution PNG.
figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)


# Save a vector-quality PDF.
figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
)


# Close the figure and release memory.
plt.close(
    figure
)


# ------------------------------------------------------------
# STEP 23: Confirm successful completion
# ------------------------------------------------------------

print()
print("Graph created successfully.")

print()
print("PNG:")
print(OUTPUT_PNG)

print()
print("PDF:")
print(OUTPUT_PDF)