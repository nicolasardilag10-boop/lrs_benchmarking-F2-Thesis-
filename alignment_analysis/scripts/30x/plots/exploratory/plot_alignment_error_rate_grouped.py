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
    / "alignment_error_rate_30x_samtools_grouped_all_mappers.png")


# Output PDF figure.
OUTPUT_PDF = (
    PROJECT
    / "alignment_analysis"
    / "figures"
    / "30x"
    / "exploratory"
    / "alignment_error_rate_30x_samtools_grouped_all_mappers.pdf")


# Output table containing the paired statistical tests.
OUTPUT_TESTS = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/statistical_tests/alignment_error_rate_30x_pairwise_tests.tsv")


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
print("Technology-matched alignment-error table:")
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

# Use the largest available error rate from the selected
# technology-matched results.
maximum_error = plot_data[METRIC].max()


# Add room for the percentage annotations above the tallest bar.
y_axis_maximum = maximum_error * 1.18


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
        "mean_first_percent": first_complete.mean(),
        "mean_second_percent": second_complete.mean(),
        "mean_difference_percent": (
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

# Separate the technology-matched bars into two panels while
# grouping the samples within each aligner.
figure, axes = plt.subplots(
    nrows=1,
    ncols=2,
    sharey=True,
    figsize=(9.1, 4.76),)

aligner_order = [
    "minimap2",
    "pbmm2",
    "VACmap",
    "VG Giraffe",]
aligner_positions = list(range(len(aligner_order)))
technology_order = [
    "ONT",
    "PacBio",]

sample_colors = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",}

sample_offsets = {
    "HG002": -0.22,
    "HG003": 0.0,
    "HG004": 0.22,}

sample_label_offsets = {
    "HG002": (-4, 4),
    "HG003": (0, 16),
    "HG004": (4, 4),}

plotting_table = (
    plot_data[
        [
            "read_technology",
            "aligner_group",
            "sample",
            METRIC,
        ]]
    .rename(
        columns={
            "read_technology": "technology",
            "aligner_group": "aligner",
            METRIC: "error_percent",})
    .sort_values(
        [
            "technology",
            "aligner",
            "sample",
        ]))

print()
print("Final plotting table:")
print(plotting_table.to_string(index=False))

expected_combinations = pd.MultiIndex.from_product(
    [
        technology_order,
        aligner_order,
        sample_order,
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
print("Missing combinations:")
if len(missing_combinations) == 0:
    print("None")
else:
    for combination in missing_combinations:
        print("  - " + " | ".join(combination))

for panel_index, technology in enumerate(technology_order):
    axis = axes[panel_index]
    bar_width = 0.19

    for aligner_index, aligner in enumerate(aligner_order):
        condition_key = (aligner, technology)
        values = (
            condition_table[condition_key]
            .reindex(sample_order)
            if condition_key in condition_table.columns
            else pd.Series(index=sample_order, dtype=float))

        for sample in sample_order:
            value = values.loc[sample]
            x_position = (
                aligner_index
                + sample_offsets[sample] * bar_width / 0.19)
            bar_value = 0 if pd.isna(value) else value
            bars = axis.bar(
                x_position,
                bar_value,
                width=bar_width,
                color=sample_colors[sample],
                edgecolor="white",
                linewidth=0.7,
                label=(
                    sample
                    if aligner_index == 0
                    else "_nolegend_"),)

            if pd.isna(value):
                axis.annotate(
                    "NA",
                    xy=(x_position, y_axis_maximum * 0.025),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="#666666",)
                continue

            axis.annotate(
                f"{value:.2f}",
                xy=(
                    bars[0].get_x() + bars[0].get_width() / 2,
                    value,),
                xytext=sample_label_offsets[sample],
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=8,
                rotation=45,
                rotation_mode="anchor",)

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
    axis.set_ylim(0, y_axis_maximum)
    axis.set_xticks(aligner_positions, aligner_order)
    axis.grid(False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("black")
    axis.spines["bottom"].set_color("black")
    axis.tick_params(axis="both", colors="black", labelsize=11)
    axis.set_facecolor("white")

axes[0].set_ylabel("Alignment error rate (%)", fontsize=13)

legend_handles, legend_labels = axes[0].get_legend_handles_labels()
figure.legend(
    handles=legend_handles,
    labels=legend_labels,
    title="GIAB sample",
    frameon=False,
    ncols=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.965),
    columnspacing=1.0,
    handletextpad=0.4,
    fontsize=9,
    title_fontsize=9,)

figure.patch.set_facecolor("white")
figure.subplots_adjust(
    left=0.08,
    right=0.98,
    bottom=0.20,
    top=0.80,
    wspace=0.12,)


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
