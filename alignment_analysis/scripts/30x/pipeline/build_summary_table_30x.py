#!/usr/bin/env python3
"""Build the alignment performance summary table (thesis Table 2).

One row per technology x aligner, values are mean +/- SD across HG002,
HG003 and HG004. Metric definitions match the final 30x figures:

- Mapped reads / bases and CIGAR-aligned yield use the verified raw input
  (the total that aligners retaining unmapped reads agree on) as the
  denominator, so VACmap's dropped unmapped reads count as unmapped.
- MQ0 reads = reads_mq0 / reads_mapped (as plot_mq0_reads_30x.py).
- Runtime = wall-clock hours; peak RAM = measured peak RSS (blank where
  not measured).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT = next(
    parent for parent in Path(__file__).resolve().parents
    if (parent / "alignment_analysis").is_dir()
)
TABLE_DIR = PROJECT / "alignment_analysis" / "tables" / "30x"
INPUT_TSV = TABLE_DIR / "final" / "alignment_benchmark_30x.tsv"
OUTPUT_TSV = TABLE_DIR / "derived" / "alignment_summary_table_30x.tsv"

ALIGNER_ORDER = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]
TECHNOLOGY_ORDER = ["ONT", "PacBio"]
TECHNOLOGY_TITLES = {"ONT": "ONT", "PacBio": "PacBio HiFi"}

# metric column -> decimals shown
METRICS = {
    "error_percent": 2,
    "mapped_reads_percent": 2,
    "mapped_bases_percent": 2,
    "cigar_yield_percent": 2,
    "mq0_reads_percent": 2,
    "runtime_hours": 2,
    "peak_ram_gb": 1,
}


def verified_raw_input(data: pd.DataFrame, column: str) -> pd.Series:
    retaining = data.loc[data["reads_unmapped"] > 0]
    agreed = retaining.groupby(["sample", "read_technology"])[column].agg(["nunique", "first"])
    if (agreed["nunique"] != 1).any():
        raise ValueError(f"Aligners retaining unmapped reads disagree on {column}:\n{agreed}")
    keys = pd.MultiIndex.from_frame(data[["sample", "read_technology"]])
    values = agreed["first"].reindex(keys).to_numpy()
    if pd.isna(values).any():
        raise ValueError(f"No retaining aligner for some sample/technology ({column}).")
    return pd.Series(values, index=data.index, dtype=float)


def mean_sd(values: pd.Series, decimals: int) -> str:
    values = values.dropna()
    if values.empty:
        return "n.m."
    if len(values) == 1:
        return f"{values.iloc[0]:.{decimals}f}"
    return f"{values.mean():.{decimals}f} $\\pm$ {values.std(ddof=1):.{decimals}f}"


def main() -> int:
    data = pd.read_csv(INPUT_TSV, sep="\t")
    data["aligner"] = data["aligner"].replace({"VACMap": "VACmap"})
    data = data.loc[data["aligner"].isin(ALIGNER_ORDER) & data["read_technology"].isin(TECHNOLOGY_ORDER)].copy()

    numeric = [
        "raw_total_sequences", "reads_mapped", "reads_unmapped", "total_length", "bases_mapped",
        "bases_mapped_cigar", "error_percent", "reads_mq0", "wallclock_runtime_hours", "peak_ram_gb", "threads",
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    raw_reads = verified_raw_input(data, "raw_total_sequences")
    raw_bases = verified_raw_input(data, "total_length")

    data["mapped_reads_percent"] = 100.0 * data["reads_mapped"] / raw_reads
    data["mapped_bases_percent"] = 100.0 * data["bases_mapped"] / raw_bases
    data["cigar_yield_percent"] = 100.0 * data["bases_mapped_cigar"] / raw_bases
    data["mq0_reads_percent"] = 100.0 * data["reads_mq0"] / data["reads_mapped"]
    data["runtime_hours"] = data["wallclock_runtime_hours"]

    rows = []
    for technology in TECHNOLOGY_ORDER:
        for aligner in ALIGNER_ORDER:
            group = data.loc[(data["read_technology"] == technology) & (data["aligner"] == aligner)]
            if len(group) != 3:
                raise ValueError(f"Expected 3 samples for {technology}/{aligner}, found {len(group)}.")
            threads = group["threads"].dropna().unique()
            row = {
                "technology": TECHNOLOGY_TITLES[technology],
                "aligner": aligner,
                "threads": int(threads[0]) if len(threads) == 1 else "/".join(str(int(t)) for t in threads),
            }
            for metric, decimals in METRICS.items():
                row[metric] = mean_sd(group[metric], decimals)
            rows.append(row)

    table = pd.DataFrame(rows)
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(table.to_string(index=False))
    print(f"\nWrote {OUTPUT_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
