from pathlib import Path
import sys
import numpy as np


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

    # SciPy performs the paired t-tests.
    from scipy.stats import ttest_rel

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


# Input table containing the combined 30x alignment statistics.
INPUT_FILE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "source/alignment_summary_30x_Samtools_Christian.tsv")


# Output PNG figure.
#
# A new name is used so the previous minimap2/pbmm2 figure
# is not overwritten.
OUTPUT_PNG = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "unmapped_bases_cigar_30x_samtools_grouped_all_mappers.png")


# Output PDF figure.
OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "unmapped_bases_cigar_30x_samtools_grouped_all_mappers.pdf")


# Output table containing the paired statistical tests.
OUTPUT_TESTS = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/statistical_tests/unmapped_bases_cigar_30x_pairwise_tests.tsv")


# Technology-specific figure outputs.
OUTPUT_TECHNOLOGY_PNG = {
    "ONT": OUTPUT_PNG.with_name(
        "unmapped_bases_cigar_30x_ont_grouped.png"),
    "PacBio": OUTPUT_PNG.with_name(
        "unmapped_bases_cigar_30x_pacbio_grouped.png"),}

OUTPUT_TECHNOLOGY_PDF = {
    "ONT": OUTPUT_PDF.with_name(
        "unmapped_bases_cigar_30x_ont_grouped.pdf"),
    "PacBio": OUTPUT_PDF.with_name(
        "unmapped_bases_cigar_30x_pacbio_grouped.pdf"),}


OUTPUT_PACBIO_PAIRED_PNG = OUTPUT_PNG.with_name(
    "unmapped_bases_cigar_30x_pacbio_aligner_paired.png")

OUTPUT_PACBIO_PAIRED_PDF = OUTPUT_PDF.with_name(
    "unmapped_bases_cigar_30x_pacbio_aligner_paired.pdf")


OUTPUT_PAIRED_PANELS_PNG = OUTPUT_PNG.with_name(
    "unmapped_bases_cigar_30x_aligner_paired_panels.png")

OUTPUT_PAIRED_PANELS_PDF = OUTPUT_PDF.with_name(
    "unmapped_bases_cigar_30x_aligner_paired_panels.pdf")


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


# A complete technology-matched benchmark contains:
#
# 3 GIAB samples
# × 2 sequencing technologies
# × 4 aligners
# = 24 rows
#
# At the current project stage, VG ONT results may still be
# unavailable. Missing observations are therefore reported as
# warnings rather than being converted to zero or stopping the
# plotting workflow.
COMPLETE_BENCHMARK_ROWS = 24


# ============================================================
# STEP 4: DEFINE THE METRIC
# ============================================================

# Plot the percentage of bases not represented by aligned CIGAR
# operations.
METRIC = "cigar_unmapped_percent"


# These columns are required for the graph.
required_columns = {
    "sample",
    "read_technology",
    "configuration",
    "total_length",
    "bases_mapped_cigar",}


# Check whether any required column is missing.
missing_columns = required_columns.difference(
    data.columns)


if missing_columns:
    raise ValueError(
        "The input table is missing these columns: "
        f"{sorted(missing_columns)}")


# Convert the values used to calculate the unmapped-base percentage.
#
# errors="raise" stops the script if a value cannot
# be converted into a number.
data["total_length"] = pd.to_numeric(
    data["total_length"],
    errors="raise",)

data["bases_mapped_cigar"] = pd.to_numeric(
    data["bases_mapped_cigar"],
    errors="raise",)


# Unaligned bases (%) =
# 100 x [1 - bases_mapped_cigar / total_length].
# This is the existing calculation expressed in equivalent form.
data[METRIC] = (
    data["total_length"]
    - data["bases_mapped_cigar"]
) / data["total_length"] * 100


# ============================================================
# STEP 5: SELECT TECHNOLOGY-MATCHED CONFIGURATIONS
# ============================================================

# For the main benchmark, use only the configuration matched
# to the sequencing technology.
#
# Minimap2:
# ONT          -> mm2-ont
# PacBio HiFi  -> mm2-pb
#
# pbmm2:
# ONT          -> pbmm2-ont
# PacBio HiFi  -> pbmm2-pb
#
# VACmap:
# ONT          -> vacmap-ont
# PacBio HiFi  -> vacmap-pb
#
# VG Giraffe:
# ONT          -> vg-ont
# PacBio HiFi  -> vg-pb
#
# This removes the cross-preset rows such as ONT with mm2-pb
# from the main technology-matched comparison.
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


# Select only the rows and columns required for the graph.
plot_data = data.loc[
    selected_conditions,
    [
        "sample",
        "read_technology",
        "configuration",
        METRIC,
    ],
].copy()


print()
print(
    "Technology-matched rows available: "
    f"{len(plot_data)} / {COMPLETE_BENCHMARK_ROWS}")


if len(plot_data) != COMPLETE_BENCHMARK_ROWS:
    print(
        "WARNING: the complete four-aligner benchmark would "
        f"contain {COMPLETE_BENCHMARK_ROWS} technology-matched "
        f"rows, but {len(plot_data)} are currently available."
    )


# ============================================================
# STEP 6: CREATE THE FOUR ALIGNER GROUPS
# ============================================================

# Convert the eight technology-specific configuration labels
# into four aligner categories used in the comparison.
configuration_to_aligner = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",

    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",

    "vacmap-ont": "VACmap",
    "vacmap-pb": "VACmap",

    "vg-ont": "VG Giraffe",
    "vg-pb": "VG Giraffe",}


plot_data["aligner_group"] = (
    plot_data["configuration"]
    .map(configuration_to_aligner))


# Stop if a selected configuration could not be mapped.
if plot_data["aligner_group"].isna().any():
    raise ValueError(
        "At least one selected configuration could not be mapped "
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

# Every sample/technology pair should have at most one result
# for each aligner in the technology-matched comparison.
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
# STEP 9: DEFINE THE SAMPLE ORDER
# ============================================================

# The three GIAB benchmark samples are kept in the same order
# used in the previous grouped figure.
sample_order = [
    "HG002",
    "HG003",
    "HG004",]


# Stop if one of the GIAB samples is completely absent.
missing_samples = set(
    sample_order
).difference(
    plot_data["sample"])


if missing_samples:
    raise ValueError(
        "These expected GIAB samples are missing completely: "
        f"{sorted(missing_samples)}")


# ============================================================
# STEP 10: CREATE ONE LOOKUP TABLE FOR ALL ALIGNERS
# ============================================================

# Reorganize the selected data so every row corresponds to one
# GIAB sample and every column corresponds to one
# aligner/technology condition.
#
# Missing combinations, such as VG ONT while those runs are
# unavailable, remain NaN. They are NEVER interpreted as 0%.
condition_table = (
    plot_data
    .pivot(
        index="sample",
        columns=[
            "aligner_group",
            "read_technology",
        ],
        values=METRIC)
    .reindex(sample_order)
)


print()
print("Technology-matched CIGAR-unmapped-base table:")
print(condition_table)


# ============================================================
# STEP 11: DEFINE THE BAR COLORS
# ============================================================

# Keep the existing minimap2/pbmm2 colors and add two related
# color pairs for VACmap and VG Giraffe.
#
# A darker/lighter pair distinguishes ONT and PacBio HiFi
# within each aligner.
condition_colors = {
    ("minimap2", "ONT"): "#0072B2",
    ("minimap2", "PacBio"): "#56B4E9",

    ("pbmm2", "ONT"): "#D55E00",
    ("pbmm2", "PacBio"): "#E69F00",

    ("VACmap", "ONT"): "#009E73",
    ("VACmap", "PacBio"): "#66C2A5",

    ("VG Giraffe", "ONT"): "#CC79A7",
    ("VG Giraffe", "PacBio"): "#AA4499",}


# Define the presentation order in the grouped graph.
condition_order = [
    ("minimap2", "ONT"),
    ("minimap2", "PacBio"),
    ("pbmm2", "ONT"),
    ("pbmm2", "PacBio"),
    ("VACmap", "ONT"),
    ("VACmap", "PacBio"),
    ("VG Giraffe", "ONT"),
    ("VG Giraffe", "PacBio"),]


# ============================================================
# STEP 12: BUILD THE AVAILABLE BAR SERIES
# ============================================================

# A complete final benchmark will produce eight series:
#
# minimap2 ONT
# minimap2 PacBio HiFi
# pbmm2 ONT
# pbmm2 PacBio HiFi
# VACmap ONT
# VACmap PacBio HiFi
# VG Giraffe ONT
# VG Giraffe PacBio HiFi
#
# If one complete condition is currently absent (for example
# VG Giraffe ONT), the series is skipped and a warning is
# printed. Individual missing sample values remain NaN.
bar_series = []


for aligner, technology in condition_order:

    condition_key = (aligner, technology)

    if condition_key not in condition_table.columns:
        print(
            "WARNING: no data available for "
            f"{aligner} {technology}; "
            "this bar series will be omitted."
        )
        continue

    values = (
        condition_table[condition_key]
        .reindex(sample_order)
    )

    if values.isna().all():
        print(
            "WARNING: all values are missing for "
            f"{aligner} {technology}; "
            "this bar series will be omitted."
        )
        continue

    technology_label = (
        "PacBio HiFi"
        if technology == "PacBio"
        else "ONT"
    )

    bar_series.append(
        (
            f"{aligner} {technology_label}",
            values,
            condition_colors[condition_key],
        )
    )


if not bar_series:
    raise ValueError(
        "No technology-matched alignment series are available "
        "for plotting.")


# ============================================================
# STEP 13: DETERMINE A COMMON Y-AXIS RANGE
# ============================================================

# Use the largest available CIGAR-unmapped percentage from the selected
# technology-matched results.
maximum_unmapped_percent = plot_data[METRIC].max()


# Add room for the percentage annotations above the tallest bar.
y_axis_maximum = maximum_unmapped_percent * 1.18


# ============================================================
# STEP 14: RUN PAIRED T-TESTS
# ============================================================

# Each test compares matched GIAB samples.
#
# Two types of comparisons are generated:
#
# 1. Pairwise aligner comparisons within the same technology.
# 2. ONT vs PacBio HiFi within the same aligner.
#
# A comparison is performed only when all three GIAB samples
# have values for both conditions. Therefore, comparisons that
# require VG ONT are automatically skipped until those results
# become available.
aligner_order = [
    "minimap2",
    "pbmm2",
    "VACmap",
    "VG Giraffe",]

technology_order = [
    "ONT",
    "PacBio",]


paired_comparisons = []


# ------------------------------------------------------------
# 14A. Compare aligners within each sequencing technology
# ------------------------------------------------------------

for technology in technology_order:
    for first_index in range(len(aligner_order)):
        for second_index in range(
            first_index + 1,
            len(aligner_order)):

            first_aligner = aligner_order[first_index]
            second_aligner = aligner_order[second_index]

            paired_comparisons.append(
                (
                    (first_aligner, technology),
                    (second_aligner, technology),
                )
            )


# ------------------------------------------------------------
# 14B. Compare ONT vs PacBio HiFi within each aligner
# ------------------------------------------------------------

for aligner in aligner_order:
    paired_comparisons.append(
        (
            (aligner, "ONT"),
            (aligner, "PacBio"),
        )
    )


test_results = []


for first_condition, second_condition in paired_comparisons:

    if first_condition not in condition_table.columns:
        continue

    if second_condition not in condition_table.columns:
        continue

    first_values = (
        condition_table[first_condition]
        .reindex(sample_order)
    )

    second_values = (
        condition_table[second_condition]
        .reindex(sample_order)
    )

    complete_pairs = (
        first_values.notna()
        & second_values.notna()
    )

    number_of_pairs = int(
        complete_pairs.sum())


    # The current statistical design uses HG002, HG003, and
    # HG004 as three matched observations. Do not run a paired
    # t-test when fewer than all three samples are available.
    if number_of_pairs != len(sample_order):
        print(
            "WARNING: paired t-test skipped for "
            f"{first_condition} vs {second_condition}: "
            f"{number_of_pairs}/3 complete GIAB pairs."
        )
        continue


    first_complete = first_values.loc[
        complete_pairs]

    second_complete = second_values.loc[
        complete_pairs]


    test = ttest_rel(
        first_complete,
        second_complete,)


    first_label = (
        f"{first_condition[0]} "
        + (
            "PacBio HiFi"
            if first_condition[1] == "PacBio"
            else "ONT"
        )
    )

    second_label = (
        f"{second_condition[0]} "
        + (
            "PacBio HiFi"
            if second_condition[1] == "PacBio"
            else "ONT"
        )
    )


    test_results.append({
        "comparison": (
            f"{first_label} vs {second_label}"),
        "n_pairs": number_of_pairs,
        "mean_first_cigar_unmapped_percent": first_complete.mean(),
        "mean_second_cigar_unmapped_percent": second_complete.mean(),
        "mean_difference_cigar_unmapped_percent": (
            first_complete
            - second_complete
        ).mean(),
        "t_statistic": test.statistic,
        "p_value": test.pvalue,})


test_results = pd.DataFrame(
    test_results)


# ------------------------------------------------------------
# 14C. Apply Holm multiple-testing correction
# ------------------------------------------------------------

if not test_results.empty:

    number_of_tests = len(
        test_results)

    ordered_indices = (
        test_results["p_value"]
        .sort_values()
        .index
        .tolist()
    )

    adjusted_values = {}
    previous_adjusted = 0.0


    # Holm's step-down procedure:
    #
    # sorted p(1) <= p(2) <= ... <= p(m)
    #
    # adjusted p(i) =
    # max(previous adjusted p,
    #     (m - i + 1) * p(i))
    #
    # and values are capped at 1.
    for rank, index in enumerate(
        ordered_indices,
        start=1,):

        raw_p_value = float(
            test_results.loc[
                index,
                "p_value"])

        multiplier = (
            number_of_tests
            - rank
            + 1)

        adjusted_p_value = min(
            1.0,
            raw_p_value * multiplier)

        adjusted_p_value = max(
            previous_adjusted,
            adjusted_p_value)

        adjusted_values[index] = (
            adjusted_p_value)

        previous_adjusted = (
            adjusted_p_value)


    test_results[
        "holm_adjusted_p_value"
    ] = pd.Series(
        adjusted_values)


    print()
    print(
        "Paired t-tests across HG002, HG003, and HG004:")
    print(
        test_results.to_string(
            index=False,
            float_format=lambda value: f"{value:.4g}"))


    # Save the statistical results so the values used in the
    # report can be traced back to the plotting workflow.
    OUTPUT_TESTS.parent.mkdir(
        parents=True,
        exist_ok=True,)

    test_results.to_csv(
        OUTPUT_TESTS,
        sep="\t",
        index=False,)

else:
    print()
    print(
        "WARNING: no complete paired comparisons were "
        "available for statistical testing.")


# ============================================================
# STEP 15: CREATE THE GROUPED FIGURE
# ============================================================

# A wider figure is used because the expanded benchmark can
# contain as many as eight mapper/technology bar series.
figure, axis = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=(
        12.5,
        6.8,),)


# ============================================================
# STEP 16: DRAW THE GROUPED BARS
# ============================================================

# Draw adjacent bars for each sample group.
x_positions = list(
    range(
        len(sample_order)))


# Scale the bar width automatically according to the number of
# currently available mapper/technology series.
number_of_series = len(
    bar_series)

group_width = 0.84

bar_width = (
    group_width
    / number_of_series)


# Center all available bars around each GIAB sample.
bar_offsets = [
    (
        series_index
        - (
            number_of_series
            - 1
        ) / 2
    )
    * bar_width

    for series_index
    in range(
        number_of_series)
]


for (
    series_name,
    values,
    color,
), offset in zip(
    bar_series,
    bar_offsets,):


    x_for_series = [
        position + offset
        for position
        in x_positions]


    bars = axis.bar(
        x_for_series,
        values.to_numpy(),
        width=bar_width * 0.94,
        color=color,
        edgecolor="white",
        linewidth=0.6,
        label=series_name,)


    # Keep the percentage annotation above every available bar.
    #
    # Missing values are intentionally not labelled because they
    # represent unavailable data, not a zero percentage.
    for bar, value in zip(
        bars,
        values):

        if pd.isna(
            value):
            continue

        axis.annotate(
            f"{value:.2f}%",
            xy=(
                bar.get_x()
                + bar.get_width()
                / 2,
                value),
            xytext=(
                0,
                4),
            textcoords=(
                "offset points"),
            horizontalalignment=(
                "center"),
            verticalalignment=(
                "bottom"),
            fontsize=10,
            rotation=0,)


# ============================================================
# STEP 17: FORMAT THE UNMAPPED-BASE GROUPED PLOT
# ============================================================

axis.set_ylim(
    0,
    y_axis_maximum,)

axis.set_xticks(
    x_positions,
    sample_order,)


# Keep all mapper/technology labels visible in a dedicated strip
# above the plotting area.
figure.legend(
    frameon=False,
    ncols=4,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.99),
    columnspacing=1.5,
    fontsize=10,)


axis.grid(False)
axis.spines["left"].set_color("black")
axis.spines["right"].set_visible(False)
axis.spines["top"].set_visible(False)
axis.tick_params(
    axis="both",
    colors="black",
    labelsize=10,)
axis.set_facecolor("white")


# ============================================================
# STEP 18: ADD SHARED TITLE AND X-AXIS LABEL
# ============================================================

# Place the title near the upper-left side of the figure,
# similar to the reference image.
#
# The title remains commented out, as in the previous version.
#figure.suptitle(
#    "CIGAR-unmapped bases per file (%)",
#    x=0.5,
#    y=0.98,
#    horizontalalignment="center",
#    fontsize=15,)


# Add the x-axis label below the grouped bars.
figure.supxlabel(
    "GIAB sample",
    fontsize=12,
    y=0.05,)


axis.set_ylabel(
    "CIGAR-unmapped bases (%)",
    fontsize=11,)


# White figure background.
figure.patch.set_facecolor(
    "white")


# ============================================================
# STEP 19: ADJUST SPACING
# ============================================================

figure.subplots_adjust(
    left=0.10,
    right=0.98,
    bottom=0.14,
    top=0.80,)


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
# STEP 21: CREATE TECHNOLOGY-SPECIFIC FIGURES
# ============================================================

def create_technology_plot(technology):
    """Create one unmapped-base percentage plot for one technology."""

    technology_bar_series = []

    for aligner in [
        "minimap2",
        "pbmm2",
        "VACmap",
        "VG Giraffe",
    ]:

        condition_key = (
            aligner,
            technology,)

        if condition_key not in condition_table.columns:
            continue

        values = (
            condition_table[condition_key]
            .reindex(sample_order))

        if values.isna().all():
            continue

        technology_label = (
            "PacBio HiFi"
            if technology == "PacBio"
            else "ONT")

        technology_bar_series.append(
            (
                f"{aligner} {technology_label}",
                values,
                condition_colors[condition_key],))

    if not technology_bar_series:
        print(
            f"WARNING: no data available for {technology}; "
            "technology-specific plot will be skipped.")
        return

    technology_figure, technology_axis = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(
            12.5,
            6.8,),)

    number_of_series = len(
        technology_bar_series)

    bar_width = (
        0.84
        / number_of_series)

    bar_offsets = [
        (
            series_index
            - (
                number_of_series
                - 1
            ) / 2
        )
        * bar_width
        for series_index
        in range(number_of_series)]

    x_positions = list(
        range(len(sample_order)))

    for (
        series_name,
        values,
        color,
    ), offset in zip(
        technology_bar_series,
        bar_offsets,):

        x_for_series = [
            position + offset
            for position in x_positions]

        bars = technology_axis.bar(
            x_for_series,
            values.to_numpy(),
            width=bar_width * 0.94,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            label=series_name,)

        for bar, value in zip(
            bars,
            values):

            if pd.isna(value):
                continue

            technology_axis.annotate(
                f"{value:.2f}%",
                xy=(
                    bar.get_x()
                    + bar.get_width()
                    / 2,
                    value),
                xytext=(
                    0,
                    4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,)

    maximum_unmapped_percent = plot_data.loc[
        plot_data["read_technology"] == technology,
        METRIC,
    ].max()

    technology_axis.set_ylim(
        0,
        maximum_unmapped_percent * 1.18,)
    technology_axis.set_xticks(
        x_positions,
        sample_order,)
    technology_axis.set_ylabel(
        "CIGAR-unmapped bases (%)",
        fontsize=11,)
    technology_axis.set_xlabel(
        "GIAB sample",
        fontsize=12,)

    technology_figure.legend(
        frameon=False,
        ncols=4,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        columnspacing=1.5,
        fontsize=10,)

    technology_axis.grid(False)
    technology_axis.spines["left"].set_color("black")
    technology_axis.spines["top"].set_visible(False)
    technology_axis.spines["right"].set_visible(False)
    technology_axis.tick_params(
        axis="both",
        colors="black",
        labelsize=10,)
    technology_axis.set_facecolor("white")

    technology_figure.patch.set_facecolor("white")
    technology_figure.subplots_adjust(
        left=0.10,
        right=0.98,
        bottom=0.14,
        top=0.80,)

    technology_figure.savefig(
        OUTPUT_TECHNOLOGY_PNG[technology],
        dpi=300,
        bbox_inches="tight",
        facecolor="white",)
    technology_figure.savefig(
        OUTPUT_TECHNOLOGY_PDF[technology],
        bbox_inches="tight",
        facecolor="white",)

    plt.close(technology_figure)

    print()
    print(f"{technology} PNG:")
    print(OUTPUT_TECHNOLOGY_PNG[technology])
    print(f"{technology} PDF:")
    print(OUTPUT_TECHNOLOGY_PDF[technology])


for technology in [
    "ONT",
    "PacBio",
]:
    create_technology_plot(technology)


# ============================================================
# STEP 22: CREATE THE PACBIO PAIRED ALIGNER FIGURE
# ============================================================

def create_pacbio_paired_plot():
    """Plot each PacBio sample across the four aligners."""

    aligner_order = [
        "minimap2",
        "pbmm2",
        "VACmap",
        "VG Giraffe",]

    sample_colors = {
        "HG002": "#0072B2",
        "HG003": "#D55E00",
        "HG004": "#009E73",}

    paired_figure, paired_axis = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(
            10.5,
            6.8,),)

    x_positions = list(
        range(len(aligner_order)))

    pacbio_values = []

    for sample in sample_order:
        values = []

        for aligner in aligner_order:
            condition_key = (
                aligner,
                "PacBio",)

            if condition_key not in condition_table.columns:
                values.append(float("nan"))
                continue

            values.append(
                condition_table.loc[
                    sample,
                    condition_key])

        values = pd.Series(
            values,
            index=aligner_order,
            dtype=float,)

        pacbio_values.append(values)

        complete_values = values.notna()
        paired_axis.plot(
            [
                x_position
                for x_position, present
                in zip(x_positions, complete_values)
                if present
            ],
            values.loc[complete_values].to_numpy(),
            color=sample_colors[sample],
            linewidth=1.1,
            alpha=0.75,
            zorder=1,)

        paired_axis.scatter(
            [
                x_position
                for x_position, present
                in zip(x_positions, complete_values)
                if present
            ],
            values.loc[complete_values].to_numpy(),
            color=sample_colors[sample],
            edgecolor="white",
            linewidth=0.8,
            s=58,
            label=sample,
            zorder=2,)

    maximum_unmapped_percent = max(
        values.max()
        for values in pacbio_values)

    paired_axis.set_ylim(
        0,
        maximum_unmapped_percent * 1.18,)
    paired_axis.set_xticks(
        x_positions,
        aligner_order,)
    paired_axis.set_ylabel(
        "CIGAR-unmapped bases (%)",
        fontsize=11,)
    paired_axis.set_xlabel(
        "Aligner",
        fontsize=12,)
    paired_axis.set_title(
        "PacBio HiFi aligner comparison",
        fontsize=14,
        pad=18,)

    paired_axis.text(
        0.99,
        0.98,
        "Lower is better",
        transform=paired_axis.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        color="#444444",)

    paired_axis.legend(
        title="GIAB sample",
        title_fontsize=13,
        frameon=False,
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        columnspacing=1.5,
        fontsize=10,)

    paired_axis.grid(
        axis="y",
        color="#D9D9D9",
        linewidth=0.7,
        linestyle="-",)
    paired_axis.spines["top"].set_visible(False)
    paired_axis.spines["right"].set_visible(False)
    paired_axis.spines["left"].set_color("black")
    paired_axis.tick_params(
        axis="both",
        colors="black",
        labelsize=10,)
    paired_axis.set_facecolor("white")
    paired_figure.patch.set_facecolor("white")
    paired_figure.subplots_adjust(
        left=0.10,
        right=0.98,
        bottom=0.14,
        top=0.82,)

    paired_figure.savefig(
        OUTPUT_PACBIO_PAIRED_PNG,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",)
    paired_figure.savefig(
        OUTPUT_PACBIO_PAIRED_PDF,
        bbox_inches="tight",
        facecolor="white",)

    plt.close(paired_figure)

    print()
    print("PacBio paired PNG:")
    print(OUTPUT_PACBIO_PAIRED_PNG)
    print("PacBio paired PDF:")
    print(OUTPUT_PACBIO_PAIRED_PDF)


create_pacbio_paired_plot()


# ============================================================
# STEP 23: CREATE TWO-PANEL MEAN ± SD BAR FIGURE
# ============================================================

def create_mean_sd_panel_plot():
    """
    Create ONT and PacBio HiFi panels.

    Each bar represents the mean CIGAR-unaligned-base percentage
    across HG002, HG003, and HG004.

    Error bars represent ±1 sample standard deviation (SD).
    """

    aligner_order = [
        "minimap2",
        "pbmm2",
        "VACmap",
        "VG Giraffe",
    ]

    technology_order = [
        "ONT",
        "PacBio",
    ]

    # --------------------------------------------------------
    # Use the same mapper/technology colors already defined
    # earlier in the script.
    #
    # minimap2: blue
    # pbmm2: orange
    # VACmap: green
    # VG Giraffe: purple
    # --------------------------------------------------------

    # Collect all values first so a common y-axis can be used.
    all_values = plot_data.loc[
        plot_data["read_technology"].isin(technology_order),
        METRIC,
    ].dropna()

    if all_values.empty:
        raise ValueError(
            "No values are available for the mean ± SD plot."
        )

    # --------------------------------------------------------
    # Create the two panels
    # --------------------------------------------------------
    figure, axes = plt.subplots(
        nrows=1,
        ncols=2,
        sharey=True,
        figsize=(13.0, 6.8),
    )

    x_positions = np.arange(len(aligner_order))

    # Keep track of the largest mean + SD for the y-axis.
    upper_values = []

    # --------------------------------------------------------
    # Draw ONT and PacBio panels
    # --------------------------------------------------------
    for panel_index, technology in enumerate(technology_order):

        axis = axes[panel_index]

        means = []
        standard_deviations = []
        colors = []

        for aligner in aligner_order:

            condition_key = (
                aligner,
                technology,
            )

            # ----------------------------------------------
            # Handle completely missing aligner/technology
            # combinations.
            # ----------------------------------------------
            if condition_key not in condition_table.columns:

                means.append(np.nan)
                standard_deviations.append(np.nan)

                # Still assign the intended color.
                colors.append(
                    condition_colors[condition_key]
                )

                continue

            values = (
                condition_table[condition_key]
                .reindex(sample_order)
                .dropna()
                .astype(float)
            )

            n_samples = len(values)

            if n_samples == 0:
                means.append(np.nan)
                standard_deviations.append(np.nan)

            else:
                mean_value = values.mean()

                # Sample SD across HG002/HG003/HG004.
                #
                # ddof=1 gives the conventional sample
                # standard deviation.
                if n_samples >= 2:
                    sd_value = values.std(ddof=1)
                else:
                    sd_value = 0.0

                means.append(mean_value)
                standard_deviations.append(sd_value)

                upper_values.append(
                    mean_value + sd_value)

            colors.append(
                condition_colors[condition_key])

        means = np.array(
            means,
            dtype=float,)

        standard_deviations = np.array(
            standard_deviations,
            dtype=float,)

        # ----------------------------------------------------
        # Draw each bar separately so missing conditions can
        # be represented explicitly as NA.
        # ----------------------------------------------------
        for (
            x_position,
            aligner,
            mean_value,
            sd_value,
            color,
        ) in zip(
            x_positions,
            aligner_order,
            means,
            standard_deviations,
            colors,
        ):

            if np.isnan(mean_value):

                axis.text(
                    x_position,
                    all_values.max() * 0.04,
                    "NA",
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    color="#666666",
                )

                continue

            # ------------------------------------------------
            # Bar = mean
            # Error bar = ±1 SD
            # ------------------------------------------------
            bar = axis.bar(
                x_position,
                mean_value,
                width=0.62,
                color=color,
                edgecolor="white",
                linewidth=0.8,
                yerr=sd_value,
                capsize=6,
                error_kw={
                    "elinewidth": 1.4,
                    "capthick": 1.4,
                    "ecolor": "black",
                },
                zorder=2,
            )

            # ------------------------------------------------
            # Mean value above each error bar
            # ------------------------------------------------
            label_height = (
                mean_value
                + sd_value
            )

            axis.annotate(
                f"{mean_value:.2f}%",
                xy=(
                    x_position,
                    label_height,
                ),
                xytext=(
                    0,
                    8,
                ),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=10,)

        # ----------------------------------------------------
        # Panel formatting
        # ----------------------------------------------------
        axis.set_title(
            "ONT"
            if technology == "ONT"
            else "PacBio HiFi",
            fontsize=16,
            pad=18,)

        axis.text(
            -0.08,
            1.03,
            "a"
            if technology == "ONT"
            else "b",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=17,
            fontweight="bold",)

        axis.set_xticks(
            x_positions,
            aligner_order,)

        axis.set_xlabel(
            "Aligner",
            fontsize=14,)

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
            labelsize=12,)

        axis.set_facecolor("white")

    # --------------------------------------------------------
    # Common y-axis
    # --------------------------------------------------------
    if upper_values:
        maximum_with_sd = max(upper_values)
    else:
        maximum_with_sd = all_values.max()

    axes[0].set_ylim(
        0,
        maximum_with_sd * 1.22,
    )

    axes[0].set_ylabel(
    "Mapped and unmapped reads (%)",
    fontsize=13,)

    # --------------------------------------------------------
    # Add explanatory text instead of a sample legend.
    # --------------------------------------------------------
    figure.text(
        0.5,
        0.94,
        "Mean ± SD across HG002, HG003, and HG004",
        ha="center",
        va="center",
        fontsize=12,
    )

    figure.patch.set_facecolor(
        "white"
    )

    figure.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.14,
        top=0.84,
        wspace=0.12,
    )

    # --------------------------------------------------------
    # New output names
    # --------------------------------------------------------
    output_mean_sd_png = OUTPUT_PNG.with_name(
        "bar_mapped_unmapped_percentage.png"
    )

    output_mean_sd_pdf = OUTPUT_PDF.with_name(
        "bar_mapped_unmapped_percentage.png.pdf"
    )

    figure.savefig(
        output_mean_sd_png,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    figure.savefig(
        output_mean_sd_pdf,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(
        figure
    )

    # --------------------------------------------------------
    # Print numerical summary for verification
    # --------------------------------------------------------
    print()
    print("Mean ± SD summary:")

    summary_rows = []

    for technology in technology_order:

        for aligner in aligner_order:

            condition_key = (
                aligner,
                technology,
            )

            if condition_key not in condition_table.columns:
                continue

            values = (
                condition_table[condition_key]
                .reindex(sample_order)
                .dropna()
                .astype(float)
            )

            if values.empty:
                continue

            summary_rows.append({
                "technology": (
                    "PacBio HiFi"
                    if technology == "PacBio"
                    else "ONT"
                ),
                "aligner": aligner,
                "n": len(values),
                "mean_percent": values.mean(),
                "sd_percent": (
                    values.std(ddof=1)
                    if len(values) >= 2
                    else 0.0
                ),
            })

    summary_table = pd.DataFrame(
        summary_rows
    )

    print(
        summary_table.to_string(
            index=False,
            float_format=lambda value: f"{value:.4f}",
        )
    )

    print()
    print("Mean ± SD PNG:")
    print(output_mean_sd_png)

    print()
    print("Mean ± SD PDF:")
    print(output_mean_sd_pdf)


create_mean_sd_panel_plot()


# ============================================================
# STEP 24: CONFIRM SUCCESSFUL COMPLETION
# ============================================================

print()
print(
    "Figure created successfully.")


print()
print(
    "Technology-matched rows plotted:")
print(
    len(
        plot_data))


print()
print(
    "Bar series plotted:")
print(
    len(
        bar_series))


print()
print(
    "PNG:")
print(
    OUTPUT_PNG)


print()
print(
    "PDF:")
print(
    OUTPUT_PDF)


if not test_results.empty:
    print()
    print(
        "Paired-test table:")
    print(
        OUTPUT_TESTS)
