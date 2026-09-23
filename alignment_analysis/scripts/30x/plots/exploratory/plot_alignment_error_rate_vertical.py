from pathlib import Path
import sys


# ============================================================
# STEP 1: IMPORT THE REQUIRED PACKAGES
# ============================================================

# This try/except block gives a clear message when VS Code
# runs the script with the wrong Python environment.
try:
    # Matplotlib creates and saves the graph.
    import matplotlib

    # "Agg" allows the script to save figures without opening
    # a separate graphical window.
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    # Pandas reads, filters, and reorganizes the TSV table.
    import pandas as pd

except ModuleNotFoundError as error:
    raise SystemExit(
        "\nA required package is missing.\n"
        f"Missing package: {error.name}\n"
        f"Python interpreter used: {sys.executable}\n\n"
        "Select the lrs-plots Python environment in VS Code.\n"
    ) from error


# Print the interpreter used by the VS Code Run button.
#
# It should contain:
# miniconda3/envs/lrs-plots
print("Python interpreter:")
print(sys.executable)


# ============================================================
# STEP 2: DEFINE INPUT AND OUTPUT PATHS
# ============================================================

# Discover the project root from this script location.
PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())


# Input table containing the alignment statistics.
INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv")


# Output PNG figure.
OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "alignment_error_rate_30x_samtools_vertical.png")


# Output PDF figure.
OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "alignment_error_rate_30x_samtools_vertical.pdf")


# Stop if the input table does not exist.
if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"The input table was not found:\n{INPUT_FILE}")


print()
print("Input file:")
print(INPUT_FILE)


# ============================================================
# STEP 3: READ THE TSV TABLE
# ============================================================

# sep="\t" means the table columns are separated by tabs.
data = pd.read_csv(
    INPUT_FILE,
    sep="\t",)


print()
print("Original table dimensions:")
print(data.shape)


# The complete table should contain:
#
# 3 samples
# × 2 technologies
# × 4 configurations
# = 24 rows
if len(data) != 24:
    raise ValueError(
        f"Expected 24 rows, but found {len(data)}.")


# ============================================================
# STEP 4: DEFINE THE METRIC
# ============================================================

# error_percent already represents the error rate as a
# percentage.
#
# Example:
# error_rate    = 0.036
# error_percent = 3.6
METRIC = "error_percent"


# These columns are required for the graph.
required_columns = {
    "sample",
    "read_technology",
    "configuration",
    METRIC,}


# Check whether any required column is missing.
missing_columns = required_columns.difference(
    data.columns)


if missing_columns:
    raise ValueError(
        "The input table is missing these columns: "
        f"{sorted(missing_columns)}")


# Convert the selected error percentage into numeric values.
#
# errors="raise" stops the script if a value cannot
# be converted into a number.
data[METRIC] = pd.to_numeric(
    data[METRIC],
    errors="raise",)


# ============================================================
# STEP 5: SELECT TECHNOLOGY-MATCHED CONFIGURATIONS
# ============================================================

# For a fair comparison, use the preset matched to each
# sequencing technology.
#
# Minimap2:
# ONT     -> mm2-ont
# PacBio  -> mm2-pb
#
# pbmm2:
# ONT     -> pbmm2-ont
# PacBio  -> pbmm2-pb
selected_conditions = (
    (
        (data["read_technology"] == "ONT")
        & (
            data["configuration"].isin([
                "mm2-ont",
                "pbmm2-ont",])))|
    (
        (data["read_technology"] == "PacBio")
        & (
            data["configuration"].isin([
                "mm2-pb",
                "pbmm2-pb",
            ]))))


# Select only the rows and columns required for the graph.
plot_data = data.loc[
    selected_conditions,
    ["sample", "read_technology", "configuration",
        METRIC,],].copy()


# Expected:
#
# 3 samples
# × 2 technologies
# × 2 aligners
# = 12 rows
if len(plot_data) != 12:
    raise ValueError(
        f"Expected 12 selected rows, "
        f"but found {len(plot_data)}.")


# ============================================================
# STEP 6: CREATE THE TWO ALIGNER GROUPS
# ============================================================

# Convert the four configuration labels into two aligner
# categories used as the panel titles.
configuration_to_aligner = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",}


plot_data["aligner_group"] = (
    plot_data["configuration"]
    .map(configuration_to_aligner))


# Stop if a configuration could not be mapped.
if plot_data["aligner_group"].isna().any():
    raise ValueError(
        "At least one configuration could not be mapped "
        "to an aligner group.")


# ============================================================
# STEP 7: CREATE orig_file LABELS
# ============================================================

# Use shorter labels for the sequencing technologies.
technology_suffix = {
    "ONT": "ont",
    "PacBio": "pb",}


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
        technology_suffix))


if plot_data["orig_file"].isna().any():
    raise ValueError(
        "At least one orig_file label could not be created.")


# ============================================================
# STEP 8: CHECK FOR DUPLICATED RESULTS
# ============================================================

# Every input file should have one result for mm2
# and one result for pbmm2.
duplicated_rows = plot_data.duplicated(
    subset=[
        "orig_file",
        "aligner_group",
    ],
    keep=False,)


if duplicated_rows.any():
    duplicated_data = plot_data.loc[
        duplicated_rows]

    raise ValueError(
        "Duplicated input-aligner observations were found:\n"
        + duplicated_data.to_string(
            index=False))


print()
print("Data selected for plotting:")

print(
    plot_data.sort_values(
        [
            "aligner_group",
            "orig_file",
        ]
    ).to_string(
        index=False,))


# ============================================================
# STEP 9: DEFINE THE Y-AXIS ORDER
# ============================================================

# This order recreates the arrangement shown in your
# reference image.
#
# The first item appears at the top after invert_yaxis().
orig_file_order = [
    "HG004_pb",
    "HG004_ont",
    "HG003_pb",
    "HG003_ont",
    "HG002_pb",
    "HG002_ont",]


# Check that all six input labels exist.
missing_files = set(
    orig_file_order
).difference(
    plot_data["orig_file"])


if missing_files:
    raise ValueError(
        "These expected input files are missing: "
        f"{sorted(missing_files)}")


# ============================================================
# STEP 10: CREATE ONE TABLE FOR EACH ALIGNER
# ============================================================

# Create the minimap2 panel table.
mm2_data = (
    plot_data.loc[
        plot_data["aligner_group"] == "minimap2"
    ]
    .set_index("orig_file")
    .reindex(orig_file_order))


# Create the pbmm2 panel table.
pbmm2_data = (
    plot_data.loc[
        plot_data["aligner_group"] == "pbmm2"
    ]
    .set_index("orig_file")
    .reindex(orig_file_order))


# Check for missing values after reordering.
if mm2_data[METRIC].isna().any():
    raise ValueError(
        "The mm2 plotting table contains missing values.")


if pbmm2_data[METRIC].isna().any():
    raise ValueError(
        "The pbmm2 plotting table contains missing values.")


print()
print("mm2 panel:")
print(mm2_data[[METRIC]])

print()
print("pbmm2 panel:")
print(pbmm2_data[[METRIC]])


# ============================================================
# STEP 11: DEFINE THE BAR COLORS
# ============================================================

# The colors resemble your reference figure.
#
# ONT is shown in coral/red.
# PacBio is shown in teal.
technology_colors = {
    "ONT": "#F45B5B",
    "PacBio": "#20B2B8",}


# Create one color for every y-axis entry.
bar_colors = [
    technology_colors["PacBio"],
    technology_colors["ONT"],
    technology_colors["PacBio"],
    technology_colors["ONT"],
    technology_colors["PacBio"],
    technology_colors["ONT"],]


# ============================================================
# STEP 12: DETERMINE A COMMON X-AXIS RANGE
# ============================================================

# Use the same y-axis range in both panels so the bar
# heights can be compared directly.
maximum_error = plot_data[METRIC].max()


# Add 5% empty space above the tallest bar.
y_axis_maximum = maximum_error * 1.05


# ============================================================
# STEP 13: CREATE THE TWO-PANEL FIGURE
# ============================================================

figure, axes = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(
        9.5,
        8.5,),
    sharex=True,)


# Give readable names to the two plotting areas.
axis_mm2 = axes[0]
axis_pbmm2 = axes[1]


# ============================================================
# STEP 14: DRAW THE mm2 BARS
# ============================================================

axis_mm2.bar(
    x=orig_file_order,
    height=mm2_data[METRIC],
    color=bar_colors,
    edgecolor="none",)


# ============================================================
# STEP 15: DRAW THE pbmm2 BARS
# ============================================================

axis_pbmm2.bar(
    x=orig_file_order,
    height=pbmm2_data[METRIC],
    color=bar_colors,
    edgecolor="none",)


# Add error percentage above each bar.
for axis, panel_data in zip(
    axes,
    [
        mm2_data,
        pbmm2_data,
    ],):
    for x_position, value in zip(
        orig_file_order,
        panel_data[METRIC],):
        axis.annotate(
            f"{value:.2f}%",
            xy=(x_position, value),
            xytext=(0, 4),
            textcoords="offset points",
            horizontalalignment="center",
            verticalalignment="bottom",
            fontsize=8,)


# ============================================================
# STEP 16: CREATE THE FACET-STYLE PANEL HEADERS
# ============================================================

# Instead of using an ordinary title, place the aligner name
# inside a white rectangular strip with a black border.
for axis, panel_title in zip(
    axes,
    [
        "minimap2",
        "pbmm2",
    ],):
    axis.text(
        0.5,
        1.015,
        panel_title,

        # Position relative to the axis.
        transform=axis.transAxes,

        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=11,

        # Draw the facet-style rectangular header.
        bbox={
            "facecolor": "white",
            "edgecolor": "black",
            "linewidth": 1.1,

            # A wider padding creates a strip-like appearance.
            "boxstyle": "square,pad=0.45",},)


# ============================================================
# STEP 17: FORMAT BOTH PANELS
# ============================================================

for axis in axes:

    # Use the same numerical range in both panels.
    axis.set_ylim(
        0,
        y_axis_maximum,)

    # Remove gridlines to match the clean reference style.
    axis.grid(
        False)

    # Keep black left and bottom borders.
    axis.spines["left"].set_color(
        "black")

    axis.spines["bottom"].set_color(
        "black")

    # Remove unnecessary outer borders.
    axis.spines["top"].set_visible(
        False)

    axis.spines["right"].set_visible(
        False)

    # Make the tick labels black.
    axis.tick_params(
        axis="both",
        colors="black",)

    # White plotting background.
    axis.set_facecolor(
        "white")


# The y-axis title is shared by both panels.
axis_mm2.set_ylabel(
    "Alignment error rate (%)")


# ============================================================
# STEP 18: ADD SHARED TITLE AND X-AXIS LABEL
# ============================================================

# Place the title near the upper-left side of the figure,
# similar to the reference image.
#figure.suptitle(
#    "Alignment error rate per file",
#    x=0.5,
#    y=0.98,                                           <----subtittle was taking out
#    horizontalalignment="center",
#    fontsize=15,)


# Add one shared x-axis label below both panels.
figure.supxlabel(
    "Sample",
    fontsize=12,)


# White figure background.
figure.patch.set_facecolor(
    "white")


# ============================================================
# STEP 19: ADJUST SPACING
# ============================================================

# wspace controls the space between the two panels.
#
# A small value produces the narrow separation seen
# in the reference figure.
figure.subplots_adjust(
    left=0.11,
    right=0.98,
    bottom=0.11,
    top=0.88,
    hspace=0.22,)


# ============================================================
# STEP 20: SAVE THE FIGURE
# ============================================================

# Create the output directory if necessary.
OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,)


# Save a high-resolution PNG.
figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",)


# Save a vector-quality PDF.
figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",)


# Close the figure and release memory.
plt.close(
    figure)


# ============================================================
# STEP 21: CONFIRM SUCCESSFUL COMPLETION
# ============================================================

print()
print("Figure created successfully.")

print()
print("PNG:")
print(OUTPUT_PNG)

print()
print("PDF:")
print(OUTPUT_PDF)