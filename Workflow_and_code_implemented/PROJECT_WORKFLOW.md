# Project Workflow

## End-to-end flowchart

```mermaid
flowchart TD

    A["Raw long-read data<br/>6 FASTQ.gz files<br/><br/>HG002 / HG003 / HG004<br/>ONT + PacBio"]

    subgraph QC["1. FASTQ quality control"]
        B1["Count reads<br/><br/>Input: FASTQ.gz<br/>Output: fastq_read_counts.tsv"]
        B2["Calculate read lengths<br/><br/>Median, mean, maximum"]
        B3["Calculate N50"]
        B4["Calculate Q20 and Q30"]
        B5["Combine FASTQ QC results<br/><br/>Output: TSV summary tables"]
    end

    subgraph ALIGN["2. Reference and alignment"]
        C1["Prepare reference genome<br/><br/>Input: FASTA<br/>Output: reference indexes"]
        C2["Run one test alignment<br/><br/>Pipeline validation"]
        C3["Run complete benchmark<br/><br/>3 samples × 2 technologies × 4 configurations"]
        C4["Alignment configurations<br/><br/>mm2-ont<br/>mm2-pb<br/>pbmm2-ont<br/>pbmm2-pb"]
        C5["Alignment files<br/><br/>Output: SAM / BAM"]
    end

    subgraph STATS["3. Extract alignment statistics"]
        D1["samtools flagstat<br/><br/>Output: *.flagstat.txt"]
        D2["samtools idxstats<br/><br/>Output: *.idxstats.tsv"]
        D3["samtools stats<br/><br/>Output: *.samtools_stats.txt"]
    end

    subgraph SUMMARY["4. Build the master alignment table"]
        E1["build_alignment_summary.py"]
        E2["Validate file names and configurations"]
        E3["Check duplicates and missing combinations"]
        E4["alignment_summary.tsv<br/><br/>24 rows × 22 columns"]
    end

    subgraph ANALYSIS["5. Python data preparation"]
        F1["Read TSV with Pandas"]
        F2["Convert metrics to numeric"]
        F3["Filter technologies and configurations"]
        F4["Create aligner and input-file labels"]
        F5["Reshape with pivot()"]
        F6["Validate pairs, duplicates and missing values"]
        F7["Convert columns to NumPy arrays"]
    end

    subgraph FIGURES["6. Figures"]
        G1["FASTQ QC figures<br/><br/>Median, N50, Q20/Q30"]
        G2["Mapped bases figures"]
        G3["Mapped bases CIGAR figures"]
        G4["Error-rate box plots"]
        G5["ONT versus PacBio panels"]
        G6["Per-file horizontal bar plots"]
        G7["Outputs<br/><br/>PNG + PDF"]
    end

    subgraph REPORT["7. Statistics and reporting"]
        H1["Paired t-tests<br/><br/>scipy.stats.ttest_rel"]
        H2["Interpret aligner effects"]
        H3["Interpret preset effects"]
        H4["Compare ONT and PacBio"]
        H5["Add figures and conclusions to slides"]
        H6["Commit scripts, tables and figures with Git"]
    end

    A --> B1
    A --> B2
    B2 --> B3
    A --> B4
    B1 --> B5
    B2 --> B5
    B3 --> B5
    B4 --> B5

    A --> C3
    C1 --> C2
    C2 --> C3
    C4 --> C3
    C3 --> C5

    C5 --> D1
    C5 --> D2
    C5 --> D3

    D1 --> E1
    D2 --> E1
    D3 --> E1
    E1 --> E2
    E2 --> E3
    E3 --> E4

    E4 --> F1
    F1 --> F2
    F2 --> F3
    F3 --> F4
    F4 --> F5
    F5 --> F6
    F6 --> F7

    B5 --> G1
    F7 --> G2
    F7 --> G3
    F7 --> G4
    F7 --> G5
    F7 --> G6
    G1 --> G7
    G2 --> G7
    G3 --> G7
    G4 --> G7
    G5 --> G7
    G6 --> G7

    F7 --> H1
    G7 --> H2
    G7 --> H3
    G7 --> H4
    H1 --> H2
    H2 --> H5
    H3 --> H5
    H4 --> H5
    H5 --> H6

    class A input;
    class B1,B2,B3,B4,B5 qc;
    class C1,C2,C3,C4,C5 align;
    class D1,D2,D3 stats;
    class E1,E2,E3,E4 summary;
    class F1,F2,F3,F4,F5,F6,F7 analysis;
    class G1,G2,G3,G4,G5,G6,G7 figures;
    class H1,H2,H3,H4,H5,H6 report;

    classDef input fill:#1F3A5F,color:#FFFFFF,stroke:#132A45,stroke-width:2px;
    classDef qc fill:#E9F5F5,color:#14373D,stroke:#2A7F8E,stroke-width:1.5px;
    classDef align fill:#EAF0F8,color:#1C3557,stroke:#4C78A8,stroke-width:1.5px;
    classDef stats fill:#F2F4F7,color:#39424E,stroke:#7A8793,stroke-width:1.5px;
    classDef summary fill:#EEF6EC,color:#244B2A,stroke:#5B8C5A,stroke-width:1.5px;
    classDef analysis fill:#EDF6FB,color:#16324F,stroke:#5D87A1,stroke-width:1.5px;
    classDef figures fill:#FFF5E4,color:#5A3B00,stroke:#C89B3C,stroke-width:1.5px;
    classDef report fill:#FCEEEE,color:#5A1F1F,stroke:#B85C5C,stroke-width:1.5px;
```

## Stage 1 — Raw inputs

Six compressed FASTQ files:

```text
HG002.ont.1k.fastq.gz
HG002.pb.1k.fastq.gz
HG003.ont.1k.fastq.gz
HG003.pb.1k.fastq.gz
HG004.ont.1k.fastq.gz
HG004.pb.1k.fastq.gz
```

Each file contains 1,000 reads.

## Stage 2 — FASTQ quality control

The first analysis stage produced:

- read count
- read-length distribution
- median read length
- average read length
- maximum read length
- N50
- Q20 percentage
- Q30 percentage

## Stage 3 — Reference and alignment

The reference FASTA was prepared before alignment. A single test alignment was used to validate the pipeline before running the full benchmark.

### Tested configurations

| Project label | Aligner | Preset/settings |
|---|---|---|
| `mm2-ont` | minimap2 | `map-ont` |
| `mm2-pb` | minimap2 | `map-hifi` |
| `pbmm2-ont` | pbmm2 | `SUBREAD` |
| `pbmm2-pb` | pbmm2 | `CCS/HIFI` |

```text
3 samples × 2 technologies × 4 configurations = 24 alignments
```

## Stage 4 — Alignment statistics

For every alignment:

```text
*.flagstat.txt
*.idxstats.tsv
*.samtools_stats.txt
```

## Stage 5 — Master summary table

The script:

```text
alignment_analysis/scripts/build_alignment_summary.py
```

combines all reports into:

```text
alignment_analysis/tables/alignment_summary.tsv
```

The validated table contains 24 rows and 22 columns.

## Stage 6 — Python analysis and figures

The plotting scripts:

1. read the summary table,
2. convert metrics to numeric values,
3. filter technologies and configurations,
4. create labels,
5. reshape paired data,
6. validate duplicates and missing values,
7. generate PNG and PDF figures.

## Stage 7 — Statistics and reporting

Paired t-tests were used where the same biological samples were compared across configurations. Results were then interpreted and added to slides before being committed to Git.