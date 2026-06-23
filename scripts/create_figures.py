from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel


INPUT_FILE = Path("results/mapping_comparison.csv")
OUTPUT_DIRECTORY = Path("figures/alignment_benchmark")
SUPPLEMENTARY_TABLE = Path(
    "results/supplementary_alignment_table.csv"
)


ALIGNMENT_ORDER = [
    "mm2-ont",
    "mm2-pb",
    "pbmm2-ont",
    "pbmm2-pb",
]


TECHNOLOGY_LABELS = {
    "ont": "ONT",
    "pacbio": "PacBio HiFi",
}


def load_and_prepare_data() -> pd.DataFrame:
    """Load the combined CSV and create plotting columns."""

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}\n"
            "Run compare_mappers.py first."
        )

    data = pd.read_csv(INPUT_FILE)

    required_columns = {
        "sample",
        "read_technology",
        "mapper",
        "preset",
        "total_reads",
        "mapped_reads",
        "unmapped_reads",
        "mapping_rate_percent",
        "error_rate_percent",
        "bases_mapped_cigar",
        "mismatches",
    }

    missing_columns = required_columns.difference(data.columns)

    if missing_columns:
        raise ValueError(
            "The CSV is missing these columns: "
            + ", ".join(sorted(missing_columns))
        )

    numeric_columns = [
        "total_reads",
        "mapped_reads",
        "unmapped_reads",
        "mapping_rate_percent",
        "error_rate_percent",
        "bases_mapped_cigar",
        "mismatches",
    ]

    for column in numeric_columns:
        data[column] = pd.to_numeric(
            data[column],
            errors="coerce",
        )

    mapper_abbreviation = {
        "minimap2": "mm2",
        "pbmm2": "pbmm2",
    }

    data["alignment"] = (
        data["mapper"].map(mapper_abbreviation)
        + "-"
        + data["preset"]
    )

    technology_abbreviation = data[
        "read_technology"
    ].replace(
        {"pacbio": "pb"}
    )

    data["orig_file"] = (
        data["sample"]
        + "_"
        + technology_abbreviation
    )

    data["mapped_bases_minus_mismatches"] = (
        data["bases_mapped_cigar"]
        - data["mismatches"]
    )

    global_minimum = data[
        "mapped_bases_minus_mismatches"
    ].min()

    data["mapped_bases_difference"] = (
        data["mapped_bases_minus_mismatches"]
        - global_minimum
    )

    return data


def ordered_values(
    data: pd.DataFrame,
    alignment: str,
    metric: str,
) -> pd.Series:
    """Return non-missing values for one alignment."""

    return data.loc[
        data["alignment"] == alignment,
        metric,
    ].dropna()


def create_error_rate_figures(
    data: pd.DataFrame,
) -> None:
    """
    Create one general comparison of minimap2 and pbmm2
    using all six sample-technology combinations.

    Comparisons:

    HG002 ONT:
        mm2-ont versus pbmm2-ont

    HG002 PacBio:
        mm2-pb versus pbmm2-pb

    The same comparisons are performed for HG003 and HG004.

    Every dotted line connects the same sample and sequencing
    technology between minimap2 and pbmm2.
    """

    # Keep only technology-matched presets.
    matched = data[
        (
            (data["read_technology"] == "ont")
            & (data["preset"] == "ont")
        )
        |
        (
            (data["read_technology"] == "pacbio")
            & (data["preset"] == "pb")
        )
    ].copy()

    # Create identifiers such as HG002_ont and HG002_pb.
    matched["pair_id"] = (
        matched["sample"]
        + "_"
        + matched["read_technology"].replace(
            {"pacbio": "pb"}
        )
    )

    # Convert percentages back to decimal error rates.
    # Example: 2.8% becomes 0.028.
    matched["error_rate_plot"] = (
        matched["error_rate_percent"] / 100
    )

    # Reshape so each row contains the minimap2 and pbmm2
    # values from the same original dataset.
    paired = matched.pivot(
        index="pair_id",
        columns="mapper",
        values="error_rate_plot",
    )

    paired = paired.dropna(
        subset=["minimap2", "pbmm2"]
    )

    expected_pairs = {
        "HG002_ont",
        "HG002_pb",
        "HG003_ont",
        "HG003_pb",
        "HG004_ont",
        "HG004_pb",
    }

    observed_pairs = set(paired.index)

    missing_pairs = expected_pairs - observed_pairs
    unexpected_pairs = observed_pairs - expected_pairs

    if missing_pairs:
        raise ValueError(
            "Missing paired datasets: "
            f"{sorted(missing_pairs)}"
        )

    if unexpected_pairs:
        raise ValueError(
            "Unexpected paired datasets: "
            f"{sorted(unexpected_pairs)}"
        )

    if len(paired) != 6:
        raise ValueError(
            "Expected 6 complete paired datasets, "
            f"but found {len(paired)}."
        )

    print(
        "\nDatasets included in the general "
        "error-rate comparison:"
    )

    print(paired.to_string())

    print(
        "\nNumber of complete paired datasets: "
        f"{len(paired)}"
    )

    # Paired test because each minimap2 value is matched
    # with the pbmm2 value from the same dataset.
    test_result = ttest_rel(
        paired["minimap2"],
        paired["pbmm2"],
    )

    p_value = test_result.pvalue

    minimap2_median = paired[
        "minimap2"
    ].median()

    pbmm2_median = paired[
        "pbmm2"
    ].median()

    figure, axis = plt.subplots(
        figsize=(9, 6)
    )

    positions = [1, 2]

    axis.boxplot(
        [
            paired["minimap2"],
            paired["pbmm2"],
        ],
        positions=positions,
        widths=0.55,
        tick_labels=["mm2", "pbmm2"],
        showfliers=False,
    )

    color_map = plt.get_cmap("tab10")

    # Connect the same dataset between minimap2 and pbmm2.
    for index, (pair_id, row) in enumerate(
        paired.iterrows()
    ):
        axis.plot(
            positions,
            [
                row["minimap2"],
                row["pbmm2"],
            ],
            marker="o",
            linestyle=":",
            linewidth=1.4,
            alpha=0.75,
            color=color_map(index),
            label=pair_id,
        )

    # Display median values.
    axis.text(
        1,
        minimap2_median + 0.00035,
        f"{minimap2_median:.3f}",
        horizontalalignment="center",
    )

    axis.text(
        2,
        pbmm2_median + 0.00035,
        f"{pbmm2_median:.3f}",
        horizontalalignment="center",
    )

    maximum_value = paired[
        ["minimap2", "pbmm2"]
    ].to_numpy().max()

    bracket_height = maximum_value + 0.0012

    # Statistical-comparison bracket.
    axis.plot(
        [1, 1, 2, 2],
        [
            bracket_height - 0.0003,
            bracket_height,
            bracket_height,
            bracket_height - 0.0003,
        ],
        linewidth=1,
    )

    axis.text(
        1.5,
        bracket_height + 0.00015,
        f"Paired t-test, p = {p_value:.3f}",
        horizontalalignment="center",
    )

    axis.set_title(
        "Alignment mismatch rate: "
        "minimap2 versus pbmm2"
    )

    axis.set_xlabel("Aligner")

    axis.set_ylabel(
        "Mismatch rate among mapped bases"
    )

    axis.legend(
        title="Original dataset",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    axis.grid(
        axis="y",
        alpha=0.20,
    )

    figure.tight_layout()

    output_file = (
        OUTPUT_DIRECTORY
        / "error_rate_all_samples_matched.png"
    )

    figure.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(figure)

    print(
        "\nGeneral error-rate figure saved at: "
        f"{output_file}"
    )


def create_mapping_rate_figures(
    data: pd.DataFrame,
) -> None:
    """
    Compare minimap2 and pbmm2 using the technology-matched
    presets only:

    ONT reads -> ONT preset
    PacBio reads -> PB preset
    """

    matched = data[
        (
            (data["read_technology"] == "ont")
            & (data["preset"] == "ont")
        )
        |
        (
            (data["read_technology"] == "pacbio")
            & (data["preset"] == "pb")
        )
    ].copy()

    for technology in ["ont", "pacbio"]:
        subset = matched[
            matched["read_technology"] == technology
        ].copy()

        samples = sorted(
            subset["sample"].unique()
        )

        x_positions = np.arange(len(samples))
        bar_width = 0.36

        minimap2_values = []
        pbmm2_values = []

        for sample in samples:
            sample_data = subset[
                subset["sample"] == sample
            ]

            minimap2_value = sample_data.loc[
                sample_data["mapper"] == "minimap2",
                "mapping_rate_percent",
            ]

            pbmm2_value = sample_data.loc[
                sample_data["mapper"] == "pbmm2",
                "mapping_rate_percent",
            ]

            minimap2_values.append(
                minimap2_value.iloc[0]
                if not minimap2_value.empty
                else np.nan
            )

            pbmm2_values.append(
                pbmm2_value.iloc[0]
                if not pbmm2_value.empty
                else np.nan
            )

        figure, axis = plt.subplots(
            figsize=(8, 5)
        )

        axis.bar(
            x_positions - bar_width / 2,
            minimap2_values,
            width=bar_width,
            label="minimap2",
        )

        axis.bar(
            x_positions + bar_width / 2,
            pbmm2_values,
            width=bar_width,
            label="pbmm2",
        )

        technology_name = (
            TECHNOLOGY_LABELS[technology]
        )

        axis.set_title(
            f"Mapping rate — {technology_name} reads"
        )

        axis.set_xlabel("Sample")
        axis.set_ylabel("Mapped reads (%)")

        axis.set_xticks(x_positions)
        axis.set_xticklabels(samples)

        axis.legend(title="Mapper")

        axis.grid(
            axis="y",
            alpha=0.25,
        )

        minimum_value = np.nanmin(
            minimap2_values + pbmm2_values
        )

        axis.set_ylim(
            max(0, minimum_value - 1),
            100,
        )

        figure.tight_layout()

        output_file = (
            OUTPUT_DIRECTORY
            / f"mapping_rate_{technology}_matched.png"
        )

        figure.savefig(
            output_file,
            dpi=300,
        )

        plt.close(figure)


def create_mapped_bases_by_sample(
    data: pd.DataFrame,
) -> None:
    """
    Create one mapped-bases figure for each sample.

    The plotted value is:

        bases mapped (CIGAR)
        - mismatches
        - global minimum
    """

    for sample in sorted(
        data["sample"].unique()
    ):
        subset = data[
            data["sample"] == sample
        ].copy()

        subset["alignment"] = pd.Categorical(
            subset["alignment"],
            categories=ALIGNMENT_ORDER,
            ordered=True,
        )

        subset = subset.sort_values(
            "alignment"
        )

        figure, axis = plt.subplots(
            figsize=(8.5, 5)
        )

        for technology in ["ont", "pacbio"]:
            technology_data = subset[
                subset["read_technology"]
                == technology
            ].copy()

            technology_data = (
                technology_data.set_index(
                    "alignment"
                )
            )

            x_values = []
            y_values = []

            for position, alignment in enumerate(
                ALIGNMENT_ORDER,
                start=1,
            ):
                if alignment in technology_data.index:
                    value = technology_data.loc[
                        alignment,
                        "mapped_bases_difference",
                    ]

                    if isinstance(value, pd.Series):
                        value = value.iloc[0]

                    x_values.append(position)
                    y_values.append(value)

            axis.plot(
                x_values,
                y_values,
                marker="o",
                linestyle=":",
                label=TECHNOLOGY_LABELS[
                    technology
                ],
            )

        axis.set_title(
            f"{sample}: mapped CIGAR bases "
            "minus mismatches"
        )

        axis.set_xlabel(
            "Alignment configuration"
        )

        axis.set_ylabel(
            "Difference from global minimum "
            "mapped bases"
        )

        axis.set_xticks(
            range(
                1,
                len(ALIGNMENT_ORDER) + 1,
            )
        )

        axis.set_xticklabels(
            ALIGNMENT_ORDER
        )

        axis.legend(
            title="Read technology"
        )

        axis.grid(
            axis="y",
            alpha=0.25,
        )

        figure.tight_layout()

        output_file = (
            OUTPUT_DIRECTORY
            / f"mapped_bases_{sample}.png"
        )

        figure.savefig(
            output_file,
            dpi=300,
        )

        plt.close(figure)


def create_mapped_bases_by_mapper(
    data: pd.DataFrame,
) -> None:
    """
    Create one box-and-point figure for minimap2
    and another for pbmm2.
    """

    for mapper in ["minimap2", "pbmm2"]:
        subset = data[
            data["mapper"] == mapper
        ].copy()

        mapper_prefix = (
            "mm2"
            if mapper == "minimap2"
            else "pbmm2"
        )

        mapper_order = [
            f"{mapper_prefix}-ont",
            f"{mapper_prefix}-pb",
        ]

        boxplot_values = [
            ordered_values(
                subset,
                alignment,
                "mapped_bases_difference",
            )
            for alignment in mapper_order
        ]

        figure, axis = plt.subplots(
            figsize=(7, 5)
        )

        axis.boxplot(
            boxplot_values,
            tick_labels=mapper_order,
            showfliers=False,
        )

        for position, alignment in enumerate(
            mapper_order,
            start=1,
        ):
            alignment_data = subset[
                subset["alignment"] == alignment
            ].sort_values("orig_file")

            if alignment_data.empty:
                continue

            offsets = np.linspace(
                -0.08,
                0.08,
                len(alignment_data),
            )

            axis.scatter(
                position + offsets,
                alignment_data[
                    "mapped_bases_difference"
                ],
            )

        axis.set_title(
            "Mapped CIGAR bases by preset — "
            f"{mapper}"
        )

        axis.set_xlabel(
            "Alignment configuration"
        )

        axis.set_ylabel(
            "Difference from global minimum "
            "mapped bases"
        )

        axis.grid(
            axis="y",
            alpha=0.25,
        )

        figure.tight_layout()

        output_file = (
            OUTPUT_DIRECTORY
            / f"mapped_bases_{mapper}.png"
        )

        figure.savefig(
            output_file,
            dpi=300,
        )

        plt.close(figure)


def create_supplementary_table(
    data: pd.DataFrame,
) -> None:
    """Create the 30x alignment supplementary table."""

    columns = [
        "filename",
        "sample",
        "read_technology",
        "coverage",
        "reference",
        "mapper",
        "preset",
        "total_reads",
        "mapped_reads",
        "unmapped_reads",
        "mapping_rate_percent",
        "error_rate_percent",
        "average_read_length",
        "average_quality",
        "bases_mapped_cigar",
        "mismatches",
    ]

    available_columns = [
        column
        for column in columns
        if column in data.columns
    ]

    supplementary = data[
        available_columns
    ].sort_values(
        [
            "sample",
            "read_technology",
            "mapper",
            "preset",
        ]
    )

    supplementary.to_csv(
        SUPPLEMENTARY_TABLE,
        index=False,
    )


def main() -> None:
    """Create all figures and the supplementary table."""

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = load_and_prepare_data()

    create_error_rate_figures(data)
    create_mapping_rate_figures(data)
    create_mapped_bases_by_sample(data)
    create_mapped_bases_by_mapper(data)
    create_supplementary_table(data)

    print("\nFigures created successfully.")

    print(
        f"Figures saved in: "
        f"{OUTPUT_DIRECTORY}"
    )

    print(
        "Supplementary table saved at: "
        f"{SUPPLEMENTARY_TABLE}"
    )


if __name__ == "__main__":
    main()
