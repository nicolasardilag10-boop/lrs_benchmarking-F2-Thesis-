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
    / "final/02_aligned_base_yield_30x.png")
OUTPUT_PDF = OUTPUT_PNG.with_suffix(".pdf")

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
    if INPUT_XLSX.exists():
        print(f"Input file: {INPUT_XLSX}")
        return pd.read_excel(INPUT_XLSX)

    if INPUT_TSV.exists():
        print(f"Input file: {INPUT_TSV}")
        return pd.read_csv(INPUT_TSV, sep="\t")

    raise FileNotFoundError(
        "Neither the requested Excel file nor the equivalent TSV was found:\n"
        f"{INPUT_XLSX}\n{INPUT_TSV}")


def prepare_data(data: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    required_columns = {
        "sample",
        "read_technology",
        "configuration",
        "total_length",
        "bases_mapped_cigar",}
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
            "total_length",
            "bases_mapped_cigar",
        ] + (["input_fastq_total_bases"]
             if "input_fastq_total_bases" in data.columns else []),
    ].copy()
    plot_data["aligner"] = plot_data["configuration"].map(
        CONFIGURATION_TO_ALIGNER)

    for column in ["total_length", "bases_mapped_cigar"]:
        plot_data[column] = pd.to_numeric(
            plot_data[column],
            errors="raise")

    uses_fastq_bases = "input_fastq_total_bases" in plot_data.columns
    if uses_fastq_bases:
        plot_data["input_fastq_total_bases"] = pd.to_numeric(
            plot_data["input_fastq_total_bases"],
            errors="raise")
        if plot_data["input_fastq_total_bases"].isna().any():
            raise ValueError(
                "input_fastq_total_bases contains missing values.")
        input_bases = plot_data["input_fastq_total_bases"]
        print("Using actual input_fastq_total_bases as denominator.")
    else:
        # Temporary diagnostic proxy only. Replace this with the actual
        # original FASTQ base count for the final publication calculation.
        input_bases = (
            plot_data
            .groupby(["sample", "read_technology"])["total_length"]
            .transform("max"))
        print(
            "WARNING: input_fastq_total_bases is unavailable. Using the "
            "maximum total_length across aligners as a temporary proxy.")

    plot_data["input_bases_denominator"] = input_bases
    plot_data["aligned_bases_gb"] = (
        plot_data["bases_mapped_cigar"] / 1e9)
    plot_data["unaligned_bases_gb"] = (
        (input_bases - plot_data["bases_mapped_cigar"]) / 1e9)

    duplicated_rows = plot_data.duplicated(
        subset=["sample", "read_technology", "aligner"],
        keep=False)
    if duplicated_rows.any():
        raise ValueError(
            "Duplicated sample/technology/aligner observations found:\n"
            + plot_data.loc[duplicated_rows].to_string(index=False))

    expected_rows = (
        len(SAMPLE_ORDER)
        * len(TECHNOLOGY_ORDER)
        * len(ALIGNER_ORDER))
    if len(plot_data) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} technology-matched rows, "
            f"but found {len(plot_data)}.")

    return plot_data, uses_fastq_bases


def create_figure(plot_data: pd.DataFrame, uses_fastq_bases: bool) -> None:
    metrics = [
        (
            "aligned_bases_gb",
            "CIGAR-aligned bases (Gb)"),
        (
            "unaligned_bases_gb",
            "Unaligned bases (Gb)"),]

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
        padding = max(metric_range * 0.12, metric_maximum * 0.02)
        y_minimum = metric_minimum - padding
        y_maximum = metric_maximum + padding

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
                    label=(
                        sample
                        if row_index == 0 and column_index == 0
                        else None),
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

            axis.set_ylim(y_minimum, y_maximum)
            axis.set_xticks(x_positions, ALIGNER_ORDER)
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
                labelsize=10)
            axis.set_facecolor("white")

            if row_index == 0:
                axis.set_title(
                    "PacBio HiFi" if technology == "PacBio" else technology,
                    fontsize=15,
                    pad=12)
            if column_index == 0:
                axis.set_ylabel(y_label, fontsize=12)

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
        fontsize=11,
        title_fontsize=12)
    figure.text(
        0.5,
        0.015,
        "Denominator: actual FASTQ bases"
        if uses_fastq_bases
        else "TEMPORARY DIAGNOSTIC PROXY: denominator = max(total_length) per sample/technology",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#555555")
    figure.patch.set_facecolor("white")
    figure.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.10,
        top=0.90,
        hspace=0.28,
        wspace=0.12)

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    figure.savefig(OUTPUT_PDF, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    print("Python interpreter:")
    print(sys.executable)
    plot_data, uses_fastq_bases = prepare_data(read_input_table())
    create_figure(plot_data, uses_fastq_bases)
    print()
    print(f"Created PNG: {OUTPUT_PNG}")
    print(f"Created PDF: {OUTPUT_PDF}")


if __name__ == "__main__":
    main()