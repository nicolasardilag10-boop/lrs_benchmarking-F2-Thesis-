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


PROJECT = next(parent for parent in Path(__file__).resolve().parents if (parent / "alignment_analysis").is_dir())

# Prefer the requested Excel file when it exists. The current project
# contains the equivalent TSV, so the script uses that file as a fallback.
INPUT_XLSX = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "final/alignment_benchmark_30x.xlsx")
INPUT_TSV = (
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
    / "final/03_read_mapping_yield_30x.png")
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")
OUTPUT_SVG = OUTPUT_PNG.with_suffix(".svg")
OUTPUT_TABLE = (
    PROJECT
    / "alignment_analysis"
    / "tables"
    / "30x"
    / "derived/plot_data/read_retention_mapping_yield_30x.tsv")

ALIGNER_ORDER = [
    "minimap2",
    "pbmm2",
    "VACMap",
    "VG Giraffe",]
TECHNOLOGY_ORDER = [
    "ONT",
    "PacBio",]
SAMPLE_ORDER = [
    "HG002",
    "HG003",
    "HG004",]

CONFIGURATION_TO_ALIGNER = {
    "mm2-ont": "minimap2",
    "mm2-pb": "minimap2",
    "pbmm2-ont": "pbmm2",
    "pbmm2-pb": "pbmm2",
    "vacmap-ont": "VACMap",
    "vacmap-pb": "VACMap",
    "vg-ont": "VG Giraffe",
    "vg-pb": "VG Giraffe",}

SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",}
SAMPLE_MARKERS = {
    "HG002": "o",
    "HG003": "^",
    "HG004": "s",}


def read_input_table() -> pd.DataFrame:
    """Read the Excel input or the equivalent tab-separated table."""
    if INPUT_XLSX.exists():
        print(f"Input file: {INPUT_XLSX}")
        return pd.read_excel(INPUT_XLSX)

    if INPUT_TSV.exists():
        print(f"Input file: {INPUT_TSV}")
        return pd.read_csv(INPUT_TSV, sep="\t")

    raise FileNotFoundError(
        "Neither the requested Excel file nor the equivalent TSV was found:\n"
        f"{INPUT_XLSX}\n{INPUT_TSV}")


def prepare_data(data: pd.DataFrame) -> pd.DataFrame:
    required_columns = {
        "sample",
        "read_technology",
        "configuration",
        "raw_total_sequences",
        "reads_mapped",}
    missing_columns = required_columns.difference(data.columns)

    if missing_columns:
        raise ValueError(
            "The input table is missing these required columns: "
            f"{sorted(missing_columns)}")

    selected_conditions = (
        (
            (data["read_technology"] == "ONT")
            & data["configuration"].isin([
                "mm2-ont",
                "pbmm2-ont",
                "vacmap-ont",
                "vg-ont",])
        )
        |
        (
            (data["read_technology"] == "PacBio")
            & data["configuration"].isin([
                "mm2-pb",
                "pbmm2-pb",
                "vacmap-pb",
                "vg-pb",])
        ))

    plot_data = data.loc[
        selected_conditions,
        [
            "sample",
            "read_technology",
            "configuration",
            "raw_total_sequences",
            "reads_mapped",
        ],
    ].copy()

    plot_data["aligner"] = (
        plot_data["configuration"]
        .map(CONFIGURATION_TO_ALIGNER))

    for column in ["raw_total_sequences", "reads_mapped"]:
        plot_data[column] = pd.to_numeric(
            plot_data[column],
            errors="raise")

    duplicated_rows = plot_data.duplicated(
        subset=["sample", "read_technology", "aligner"],
        keep=False)
    if duplicated_rows.any():
        raise ValueError(
            "Duplicated sample/technology/aligner observations found:\n"
            + plot_data.loc[duplicated_rows].to_string(index=False))

    expected_rows = len(SAMPLE_ORDER) * len(TECHNOLOGY_ORDER) * len(ALIGNER_ORDER)
    if len(plot_data) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} technology-matched rows, "
            f"but found {len(plot_data)}.")

    # Diagnostic proxy only: until the true FASTQ read count is available,
    # use the largest aligner-specific raw count for each sample/technology.
    # Replace this with the actual original FASTQ read count for final results.
    input_proxy = (
        plot_data
        .groupby(["sample", "read_technology"])["raw_total_sequences"]
        .transform("max"))
    plot_data["input_reads_proxy"] = input_proxy
    plot_data["read_retention_percent"] = (
        100 * plot_data["raw_total_sequences"]
        / plot_data["input_reads_proxy"])
    plot_data["input_normalized_mapped_reads_percent"] = (
        100 * plot_data["reads_mapped"]
        / plot_data["input_reads_proxy"])

    return plot_data[
        [
            "sample",
            "read_technology",
            "aligner",
            "raw_total_sequences",
            "reads_mapped",
            "input_reads_proxy",
            "read_retention_percent",
            "input_normalized_mapped_reads_percent",
        ]
    ]


def print_summary(plot_data: pd.DataFrame) -> None:
    summary = (
        plot_data
        .groupby(["read_technology", "aligner"], sort=False)
        .agg(
            retention_mean=("read_retention_percent", "mean"),
            retention_min=("read_retention_percent", "min"),
            retention_max=("read_retention_percent", "max"),
            mapping_mean=(
                "input_normalized_mapped_reads_percent",
                "mean"),
            mapping_min=(
                "input_normalized_mapped_reads_percent",
                "min"),
            mapping_max=(
                "input_normalized_mapped_reads_percent",
                "max"),
        )
        .reset_index())
    summary["aligner"] = pd.Categorical(
        summary["aligner"],
        categories=ALIGNER_ORDER,
        ordered=True)
    summary = summary.sort_values(
        ["read_technology", "aligner"])

    print()
    print("Mean and range across HG002, HG003, and HG004:")
    print(summary.to_string(index=False, float_format="%.4f"))


def plot_data_figure(plot_data: pd.DataFrame) -> None:
    metrics = [
        (
            "read_retention_percent",
            "Reads retained from input (%)"),
        (
            "input_normalized_mapped_reads_percent",
            "Input-normalized mapped reads (%)"),]

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14.0, 10.0),
        sharex=True)

    x_positions = list(range(len(ALIGNER_ORDER)))

    for row_index, (metric, y_label) in enumerate(metrics):
        metric_values = plot_data[metric].to_numpy(dtype=float)
        metric_minimum = metric_values.min()
        metric_maximum = metric_values.max()
        metric_range = metric_maximum - metric_minimum
        padding = max(metric_range * 0.12, 0.15)
        y_minimum = min(100 - padding, metric_minimum - padding)
        y_maximum = max(100 + padding, metric_maximum + padding)

        for column_index, technology in enumerate(TECHNOLOGY_ORDER):
            axis = axes[row_index, column_index]
            technology_data = plot_data.loc[
                plot_data["read_technology"] == technology]

            for sample in SAMPLE_ORDER:
                sample_data = (
                    technology_data.loc[
                        technology_data["sample"] == sample]
                    .set_index("aligner")
                    .reindex(ALIGNER_ORDER))
                values = sample_data[metric].to_numpy(dtype=float)

                axis.scatter(
                    x_positions,
                    values,
                    s=58,
                    color=SAMPLE_COLORS[sample],
                    marker=SAMPLE_MARKERS[sample],
                    edgecolor="black",
                    linewidth=0.35,
                    label=sample if row_index == 0 and column_index == 0 else None,
                    zorder=3)

            mean_values = (
                technology_data
                .groupby("aligner")[metric]
                .mean()
                .reindex(ALIGNER_ORDER))
            for x_position, mean_value in zip(
                    x_positions,
                    mean_values.to_numpy(dtype=float)):
                axis.plot(
                    x_position,
                    mean_value,
                    color="#333333",
                    marker="_",
                    markersize=18,
                    markeredgewidth=2.4,
                    linestyle="none",
                    label="_nolegend_",
                    zorder=4)

            axis.axhline(
                100,
                color="#555555",
                linewidth=1.0,
                linestyle="--",
                zorder=1)
            axis.set_ylim(y_minimum, y_maximum)
            axis.set_xticks(x_positions, ALIGNER_ORDER, rotation=45)
            axis.grid(
                axis="y",
                color="#D9D9D9",
                linewidth=0.7,
                linestyle="-",
                zorder=0)
            axis.set_axisbelow(True)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.spines["left"].set_color("black")
            axis.spines["bottom"].set_color("black")
            axis.tick_params(
                axis="both",
                colors="black",
                labelsize=12)
            axis.set_facecolor("white")

            if row_index == 0:
                axis.set_title(
                    "PacBio HiFi" if technology == "PacBio" else technology,
                    fontsize=15,
                    pad=12)
            if column_index == 0:
                axis.set_ylabel(y_label, fontsize=13)

        axes[row_index, 0].text(
            -0.08,
            1.04,
            "A" if row_index == 0 else "B",
            transform=axes[row_index, 0].transAxes,
            fontsize=17,
            fontweight="bold",
            ha="left",
            va="bottom")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        title="GIAB sample",
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        frameon=False,
        fontsize=12,
        title_fontsize=13)
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.08,
        top=0.90,
        hspace=0.28,
        wspace=0.12)

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    figure.savefig(OUTPUT_PDF, bbox_inches="tight")
    figure.savefig(OUTPUT_SVG, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    print("Python interpreter:")
    print(sys.executable)
    data = prepare_data(read_input_table())
    data.to_csv(OUTPUT_TABLE, sep="\t", index=False)
    print_summary(data)
    plot_data_figure(data)
    print()
    print("Created figure files:")
    print(OUTPUT_PNG)
    print(OUTPUT_PDF)
    print(OUTPUT_SVG)
    print("Created plotting table:")
    print(OUTPUT_TABLE)


if __name__ == "__main__":
    main()