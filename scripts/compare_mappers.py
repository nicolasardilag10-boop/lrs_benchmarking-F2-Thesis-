from pathlib import Path
import csv
import re


# Folder containing the samtools statistics files
INPUT_DIRECTORY = Path("samtools_stats_30x_Christian")

# Folder and file where the combined results will be saved
OUTPUT_DIRECTORY = Path("results")
OUTPUT_FILE = OUTPUT_DIRECTORY / "mapping_comparison.csv"


# Metrics that we want to extract from every samtools stats file
METRICS_TO_EXTRACT = {
    "raw total sequences": "total_reads",
    "reads mapped": "mapped_reads",
    "reads unmapped": "unmapped_reads",
    "reads MQ0": "mq0_reads",
    "non-primary alignments": "non_primary_alignments",
    "supplementary alignments": "supplementary_alignments",
    "bases mapped (cigar)": "bases_mapped_cigar",
    "mismatches": "mismatches",
    "error rate": "error_rate",
    "average length": "average_read_length",
    "average quality": "average_quality",
}


def convert_value(value: str):
    """
    Convert a text value to an integer or float when possible.
    Otherwise, return it as text.
    """
    value = value.strip()

    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def read_samtools_stats(file_path: Path) -> dict:
    """
    Read the SN lines from one samtools stats file.
    """
    extracted_metrics = {}

    with file_path.open("r", encoding="utf-8") as stats_file:
        for line in stats_file:
            line = line.strip()

            if not line.startswith("SN"):
                continue

            # Remove the initial SN text
            content = line.removeprefix("SN").strip()

            if ":" not in content:
                continue

            metric_name, value = content.split(":", maxsplit=1)

            metric_name = metric_name.strip()

            # Remove comments beginning with #
            value = value.split("#", maxsplit=1)[0].strip()

            if metric_name in METRICS_TO_EXTRACT:
                output_column = METRICS_TO_EXTRACT[metric_name]
                extracted_metrics[output_column] = convert_value(value)

    return extracted_metrics


def parse_filename(file_path: Path) -> dict:
    """
    Extract sample, sequencing technology, reference,
    mapper and preset from the filename.

    Example:
    HG002_ont_30x.hg38.mm2-ont.cram.stats.SN.txt
    """
    filename = file_path.name

    pattern = re.compile(
        r"^(?P<sample>HG\d+)_"
        r"(?P<read_technology>ont|pb)_"
        r"(?P<coverage>\d+x)\."
        r"(?P<reference>[^.]+)\."
        r"(?P<mapper_preset>mm2|pbmm2)-"
        r"(?P<preset>ont|pb)\."
    )

    match = pattern.search(filename)

    if not match:
        raise ValueError(f"Could not interpret filename: {filename}")

    metadata = match.groupdict()

    if metadata["mapper_preset"] == "mm2":
        metadata["mapper"] = "minimap2"
    else:
        metadata["mapper"] = "pbmm2"

    # Replace the abbreviated PacBio label
    if metadata["read_technology"] == "pb":
        metadata["read_technology"] = "pacbio"

    return metadata


def calculate_derived_metrics(row: dict) -> dict:
    """
    Calculate percentages from the raw statistics.
    """
    total_reads = row.get("total_reads", 0)
    mapped_reads = row.get("mapped_reads", 0)
    unmapped_reads = row.get("unmapped_reads", 0)
    mq0_reads = row.get("mq0_reads", 0)
    error_rate = row.get("error_rate", 0)

    if total_reads:
        row["mapping_rate_percent"] = mapped_reads / total_reads * 100
        row["unmapped_rate_percent"] = unmapped_reads / total_reads * 100
    else:
        row["mapping_rate_percent"] = None
        row["unmapped_rate_percent"] = None

    if mapped_reads:
        row["mq0_rate_percent"] = mq0_reads / mapped_reads * 100
    else:
        row["mq0_rate_percent"] = None

    row["error_rate_percent"] = error_rate * 100

    return row


def collect_all_results(input_directory: Path) -> list[dict]:
    """
    Process every samtools SN statistics file.
    """
    results = []

    files = sorted(input_directory.glob("*.cram.stats.SN.txt"))

    if not files:
        raise FileNotFoundError(
            f"No statistics files were found in {input_directory}"
        )

    for file_path in files:
        metadata = parse_filename(file_path)
        metrics = read_samtools_stats(file_path)

        row = {
            "filename": file_path.name,
            **metadata,
            **metrics,
        }

        row = calculate_derived_metrics(row)
        results.append(row)

    return results


def save_results_to_csv(results: list[dict], output_file: Path):
    """
    Save all rows in a CSV file.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    column_order = [
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
        "unmapped_rate_percent",
        "mq0_reads",
        "mq0_rate_percent",
        "non_primary_alignments",
        "supplementary_alignments",
        "bases_mapped_cigar",
        "mismatches",
        "error_rate",
        "error_rate_percent",
        "average_read_length",
        "average_quality",
    ]

    with output_file.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=column_order,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(results)


def print_summary(results: list[dict]):
    """
    Print a simple summary in the terminal.
    """
    print(f"\nProcessed {len(results)} statistics files.\n")

    for row in results:
        print(
            f"{row['sample']:5} | "
            f"{row['read_technology']:7} | "
            f"{row['mapper']:8} | "
            f"preset={row['preset']:3} | "
            f"mapped={row['mapping_rate_percent']:.3f}% | "
            f"error={row['error_rate_percent']:.3f}%"
        )


def main():
    results = collect_all_results(INPUT_DIRECTORY)
    save_results_to_csv(results, OUTPUT_FILE)
    print_summary(results)

    print(f"\nCSV saved at: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()