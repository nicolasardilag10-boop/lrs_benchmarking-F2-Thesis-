from pathlib import Path
import sys


# ============================================================
# STEP 1: IMPORT THE REQUIRED PACKAGES
# ============================================================

# This block gives a clearer message when VS Code uses
# the wrong Python interpreter.
try:
    # Matplotlib creates and saves the figure.
    import matplotlib

    # Agg saves figures without opening a separate window.
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    # NumPy handles numerical arrays and calculations.
    import numpy as np

    # Pandas reads, filters, and reshapes the TSV table.
    import pandas as pd

    # ttest_rel performs paired t-tests.
    from scipy.stats import ttest_rel

except ModuleNotFoundError as error:
    raise SystemExit(
        "\nA required Python package is missing.\n"
        f"Missing package: {error.name}\n"
        f"Python interpreter used: {sys.executable}\n\n"
        "VS Code should use:\n"
        "/home/nicolas/miniconda3/envs/lrs-plots/bin/python\n"
    ) from error


# Print the interpreter used by the VS Code Run button.
#
# The path should contain:
# miniconda3/envs/lrs-plots
print("Python interpreter:")
print(sys.executable)


# ============================================================
# STEP 2: DEFINE PROJECT AND FILE PATHS
# ============================================================

# Main project directory:
# /home/nicolas/lrs_benchmarking
PROJECT = Path.home() / "lrs_benchmarking"


# Input alignment summary table.
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
    / "alignment_error_rate_by_technology.png"
)


# Output PDF figure.
OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "alignment_error_rate_by_technology.pdf"
)


# Stop when the input table cannot be found.
if not INPUT_FILE.exists():
    raise FileNotFoundError(
        f"The input table was not found:\n{INPUT_FILE}"
    )


print()
print("Input file:")
print(INPUT_FILE)


# ============================================================
# STEP 3: READ THE TSV TABLE
# ============================================================

# sep="\t" means the table columns are separated by tabs.
data = pd.read_csv(
    INPUT_FILE,
    sep="\t",
)


print()
print("Original table dimensions:")
print(data.shape)


# The complete table should contain:
#
# 3 samples
# × 2 sequencing technologies
# × 4 configurations
# = 24 rows
if len(data) != 24:
    raise ValueError(
        f"Expected 24 rows, but found {len(data)}."
    )


# ============================================================
# STEP 4: DEFINE THE METRIC AND CATEGORY ORDER
# ============================================================

# error_rate contains decimal values.
#
# Example:
# 0.03 = 3% error
METRIC = "error_rate"


# Order of the four boxes inside each panel.
configuration_order = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]


# Order of biological samples.
sample_order = [
    "HG002",
    "HG003",
    "HG004",
]


# Order of the two sequencing-technology panels.
technology_order = [
    "ONT",
    "PacBio",
]


# Human-readable panel titles.
technology_titles = {
    "ONT": "ONT reads",
    "PacBio": "PacBio reads",
}


# Check that all required columns exist.
required_columns = {
    "sample",
    "read_technology",
    "configuration",
    METRIC,
}


missing_columns = required_columns.difference(
    data.columns
)


if missing_columns:
    raise ValueError(
        "The table is missing these required columns: "
        f"{sorted(missing_columns)}"
    )


# Convert error_rate into numeric values.
#
# errors="raise" stops the script when an invalid
# numerical value is found.
data[METRIC] = pd.to_numeric(
    data[METRIC],
    errors="raise",
)


# ============================================================
# STEP 5: SELECT THE DATA REQUIRED FOR THE FIGURE
# ============================================================

# Keep only the four required configurations.
plot_data = data.loc[
    data["configuration"].isin(
        configuration_order
    ),
    [
        "sample",
        "read_technology",
        "configuration",
        METRIC,
    ],
].copy()


# All 24 rows should remain.
if len(plot_data) != 24:
    raise ValueError(
        f"Expected 24 plotting rows, "
        f"but found {len(plot_data)}."
    )


# Check for duplicated measurements.
#
# Every sample-technology-configuration combination
# must occur exactly once.
duplicated_rows = plot_data.duplicated(
    subset=[
        "sample",
        "read_technology",
        "configuration",
    ],
    keep=False,
)


if duplicated_rows.any():
    duplicated_data = plot_data.loc[
        duplicated_rows
    ]

    raise ValueError(
        "Duplicated measurements were found:\n"
        + duplicated_data.to_string(
            index=False
        )
    )


print()
print("Data selected for plotting:")

print(
    plot_data.sort_values(
        [
            "read_technology",
            "sample",
            "configuration",
        ]
    ).to_string(
        index=False,
    )
)


# ============================================================
# STEP 6: CREATE ONE PAIRED TABLE FOR EACH TECHNOLOGY
# ============================================================

# This dictionary will contain:
#
# paired_tables["ONT"]
# paired_tables["PacBio"]
paired_tables = {}


for technology in technology_order:

    # Select only one sequencing technology.
    technology_data = plot_data.loc[
        plot_data["read_technology"] == technology
    ].copy()

    # Each technology should have:
    #
    # 3 samples × 4 configurations = 12 rows
    if len(technology_data) != 12:
        raise ValueError(
            f"Expected 12 rows for {technology}, "
            f"but found {len(technology_data)}."
        )

    # Reshape the table.
    #
    # Before:
    #
    # sample    configuration    error_rate
    # HG002     mm2-ont          value
    # HG002     mm2-pb           value
    #
    # After:
    #
    # sample    mm2-ont   mm2-pb   pbmm2-ont   pbmm2-pb
    # HG002     value     value     value         value
    paired_table = (
        technology_data.pivot(
            index="sample",
            columns="configuration",
            values=METRIC,
        )
        .reindex(
            index=sample_order,
            columns=configuration_order,
        )
    )

    # Stop if a sample or configuration is missing.
    if paired_table.isna().any().any():
        raise ValueError(
            f"The paired {technology} table "
            "contains missing values."
        )

    paired_tables[technology] = paired_table


print()
print("ONT paired table:")
print(paired_tables["ONT"])

print()
print("PacBio paired table:")
print(paired_tables["PacBio"])


# ============================================================
# STEP 7: DEFINE THE STATISTICAL COMPARISONS
# ============================================================

# Each comparison contains:
#
# first configuration
# second configuration
# x position of first box
# x position of second box
#
# These tests will be performed separately for ONT
# and PacBio.
comparisons = [
    (
        "mm2-ont",
        "pbmm2-ont",
        1.0,
        3.0,
    ),
    (
        "mm2-pb",
        "pbmm2-pb",
        2.0,
        4.0,
    ),
    (
        "mm2-ont",
        "pbmm2-pb",
        1.0,
        4.0,
    ),
]


# ============================================================
# STEP 8: DEFINE HELPER FUNCTIONS
# ============================================================

def format_p_value(p_value):
    """
    Convert a numerical p-value into a readable label.
    """

    if p_value < 0.001:
        significance = "***"
    elif p_value < 0.01:
        significance = "**"
    elif p_value < 0.05:
        significance = "*"
    else:
        significance = ""

    if p_value < 0.0001:
        base_text = "p < 0.0001"
    elif p_value < 0.01:
        base_text = f"p = {p_value:.4f}"
    else:
        base_text = f"p = {p_value:.2f}"

    if significance:
        return f"{base_text} {significance}"

    return base_text


def add_bracket(
    axis,
    x1,
    x2,
    y,
    height,
    label,
):
    """
    Draw a comparison bracket and its p-value label.

    x1 and x2:
        Positions of the compared boxes.

    y:
        Vertical starting position.

    height:
        Height of the short bracket ends.

    label:
        Statistical result written above the bracket.
    """

    # Draw this shape:
    #
    # └────────────────┘
    axis.plot(
        [
            x1,
            x1,
            x2,
            x2,
        ],
        [
            y,
            y + height,
            y + height,
            y,
        ],
        color="black",
        linewidth=1.1,
        clip_on=False,
        zorder=7,
    )

    # Place the p-value above the bracket.
    axis.text(
        (x1 + x2) / 2,
        y + height,
        label,
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=9,
        color="black",
        zorder=8,
    )


# ============================================================
# STEP 9: DETERMINE THE COMMON Y-AXIS RANGE
# ============================================================

# Combine every error-rate value from both technologies.
#
# A common y-axis allows direct visual comparison between
# ONT and PacBio.
all_values = plot_data[METRIC].to_numpy(
    dtype=float
)


global_minimum = all_values.min()
global_maximum = all_values.max()
global_range = global_maximum - global_minimum


# Protect against the unlikely situation where every
# value is identical.
if global_range == 0:
    global_range = 1.0


# First bracket position.
first_bracket_y = (
    global_maximum
    + global_range * 0.08
)


# Height of the bracket ends.
bracket_height = (
    global_range * 0.025
)


# Space between the three brackets.
bracket_spacing = (
    global_range * 0.11
)


# ============================================================
# STEP 10: CREATE THE TWO-PANEL FIGURE
# ============================================================

# sharey=True gives both panels the same y-axis scale.
#
# This is important when comparing ONT with PacBio.
figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(
        14.0,
        6.5,
    ),
    sharey=True,
)


# Horizontal box positions.
x_positions = np.array([
    1.0,
    2.0,
    3.0,
    4.0,
])


# Use one consistent color per biological sample.
color_map = plt.get_cmap(
    "tab10"
)


sample_colors = {
    sample: color_map(index)
    for index, sample in enumerate(sample_order)
}


# ============================================================
# STEP 11: DRAW EACH TECHNOLOGY PANEL
# ============================================================

for axis, technology in zip(
    axes,
    technology_order,
):

    # Retrieve the paired table for this technology.
    paired_table = paired_tables[
        technology
    ]

    # Convert the four configuration columns into
    # four NumPy arrays.
    boxplot_values = [
        paired_table[configuration]
        .to_numpy(dtype=float)
        for configuration in configuration_order
    ]

    # --------------------------------------------------------
    # Draw four white box plots
    # --------------------------------------------------------

    axis.boxplot(
        boxplot_values,
        positions=x_positions,
        widths=0.66,
        patch_artist=True,
        showfliers=False,

        # White box interiors.
        boxprops={
            "facecolor": "white",
            "edgecolor": "black",
            "linewidth": 1.2,
        },

        # Black median line.
        medianprops={
            "color": "black",
            "linewidth": 1.8,
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

    # --------------------------------------------------------
    # Draw paired sample points and dotted lines
    # --------------------------------------------------------

    for sample in sample_order:

        # Extract the four values from one biological sample.
        sample_values = (
            paired_table.loc[
                sample,
                configuration_order,
            ]
            .to_numpy(
                dtype=float
            )
        )

        sample_color = sample_colors[
            sample
        ]

        # Dotted line connecting the same sample
        # across the four configurations.
        axis.plot(
            x_positions,
            sample_values,
            linestyle=":",
            linewidth=1.3,
            color=sample_color,
            alpha=0.70,
            zorder=2,
        )

        # Individual sample points.
        axis.scatter(
            x_positions,
            sample_values,
            s=46,
            color=sample_color,
            edgecolor="black",
            linewidth=0.35,

            # Add legend labels only on the first panel
            # to avoid duplicated legend entries.
            label=(
                sample
                if technology == "ONT"
                else None
            ),

            zorder=3,
        )

    # --------------------------------------------------------
    # Run and draw the paired t-tests
    # --------------------------------------------------------

    print()
    print(f"{technology} paired t-tests:")

    for comparison_index, (
        first_configuration,
        second_configuration,
        first_x,
        second_x,
    ) in enumerate(comparisons):

        # Extract paired values in the same sample order.
        first_values = (
            paired_table[first_configuration]
            .to_numpy(dtype=float)
        )

        second_values = (
            paired_table[second_configuration]
            .to_numpy(dtype=float)
        )

        # Perform paired t-test.
        t_statistic, p_value = ttest_rel(
            first_values,
            second_values,
        )

        print(
            f"{first_configuration} vs "
            f"{second_configuration}: "
            f"t = {t_statistic:.4f}, "
            f"p = {p_value:.6f}"
        )

        # Give each bracket a separate vertical level.
        bracket_y = (
            first_bracket_y
            + comparison_index * bracket_spacing
        )

        add_bracket(
            axis=axis,
            x1=first_x,
            x2=second_x,
            y=bracket_y,
            height=bracket_height,
            label=format_p_value(
                p_value
            ),
        )

    # --------------------------------------------------------
    # Format the individual panel
    # --------------------------------------------------------

    axis.set_title(
        technology_titles[
            technology
        ],
        fontsize=12,
    )

    axis.set_xticks(
        x_positions
    )

    axis.set_xticklabels(
        configuration_order,
        rotation=0,
    )

    axis.set_xlabel(
        "Alignment configuration"
    )

    # Remove gridlines.
    axis.grid(
        False
    )

    # Remove top and right borders.
    axis.spines["top"].set_visible(
        False
    )

    axis.spines["right"].set_visible(
        False
    )

    # Add horizontal space around the boxes.
    axis.set_xlim(
        0.45,
        4.55,
    )


# ============================================================
# STEP 12: FORMAT THE SHARED FIGURE
# ============================================================

# Shared y-axis label.
axes[0].set_ylabel(
    "Alignment error rate"
)


# Use one common y-axis range.
axes[0].set_ylim(
    global_minimum - global_range * 0.12,
    global_maximum + global_range * 0.50,
)


# Shared figure title centered at the bottom.
figure.text(
    0.5,
    0.025,
    "Alignment",
    ha="center",
    va="bottom",
    fontsize=14,
)


handles, labels = axes[0].get_legend_handles_labels()


# Attach the shared legend to the PacBio panel.
#
# axes[1] means the second panel:
# axes[0] = ONT
# axes[1] = PacBio
axes[1].legend(
    handles,
    labels,
    title="Sample",

    # Position the legend centered vertically on the right
    # side of the plotting area (just outside the panel).
    bbox_to_anchor=(
        1.01,
        0.50,
    ),

    loc="center left",

    # Remove additional space between the axis and legend.
    borderaxespad=0.0,

    # Remove the legend border.
    frameon=False,
)


# Allow the two panels to use more of the figure width.
#
# right=0.91 moves the right edge of the PacBio panel
# farther right, reducing unused blank space.
figure.subplots_adjust(
    left=0.08,
    right=0.91,
    bottom=0.15,
    top=0.84,
    wspace=0.18,
)


# ============================================================
# STEP 13: SAVE THE FIGURE
# ============================================================

# Create the figures directory when necessary.
OUTPUT_PNG.parent.mkdir(
    parents=True,
    exist_ok=True,
)


# Save high-resolution PNG.
figure.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)


# Save vector PDF.
figure.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
)


# Close the figure and release memory.
plt.close(
    figure
)


# ============================================================
# STEP 14: CONFIRM SUCCESSFUL COMPLETION
# ============================================================

print()
print("Two-panel graph created successfully.")

print()
print("PNG:")
print(OUTPUT_PNG)

print()
print("PDF:")
print(OUTPUT_PDF)