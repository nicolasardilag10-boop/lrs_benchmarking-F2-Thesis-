"""Shared publication-style helpers for the 30x final figures.

Keeps one definition of the sample color mapping, font detection,
rcParams, and the small pieces of matplotlib boilerplate (spines,
gridlines, panel letters, legends, savefig) that every final figure
in alignment_analysis/scripts/30x/plots/final/ repeats. Nothing here
changes any plotted value — it only standardizes appearance.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ============================================================
# SAMPLE / ALIGNER IDENTITY (fixed across every figure)
# ============================================================

SAMPLE_ORDER = ["HG002", "HG003", "HG004"]

SAMPLE_COLORS = {
    "HG002": "#0072B2",
    "HG003": "#D55E00",
    "HG004": "#009E73",
}

ALIGNER_ORDER = ["minimap2", "pbmm2", "VACmap", "VG Giraffe"]

# Shape encodes aligner identity where it must be shown independently
# of sample color (e.g. the yield-vs-error scatter), instead of adding
# four more saturated hues.
ALIGNER_MARKERS = {
    "minimap2": "o",
    "pbmm2": "s",
    "VACmap": "^",
    "VG Giraffe": "D",
}

TECHNOLOGY_ORDER = ["ONT", "PacBio"]
TECHNOLOGY_TITLES = {"ONT": "ONT", "PacBio": "PacBio HiFi"}

NEUTRAL_GRAY = "#8C8C8C"
GRID_COLOR = "#E5E5E5"


# ============================================================
# FONT DETECTION
#
# Arial -> Helvetica -> Liberation Sans -> DejaVu Sans. Resolved once,
# programmatically, so an unavailable family in the chain never
# reaches matplotlib's findfont and never emits a warning.
# ============================================================

_FONT_PREFERENCE = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]


def detect_font_family() -> str:
    available = {font.name for font in font_manager.fontManager.ttflist}

    # The font manager's cache may predate fonts installed system-wide
    # (e.g. Liberation Sans on this host). Scan the OS font directories
    # directly and register any regular-weight match before falling
    # through the preference chain, so a real system font is used
    # instead of silently dropping to DejaVu Sans. DejaVu Sans itself
    # is always already available (it ships with matplotlib), so it
    # must not short-circuit this scan for the higher-priority fonts.
    higher_priority = set(_FONT_PREFERENCE) - {"DejaVu Sans"}
    if not available & higher_priority:
        for font_path in font_manager.findSystemFonts():
            lower_path = font_path.lower()
            if "italic" in lower_path or "oblique" in lower_path:
                continue
            try:
                name = font_manager.get_font(font_path).family_name
            except Exception:
                continue
            # Register both the regular and bold weights (panel letters
            # and titles need bold) so matplotlib never has to fall back
            # and emit a "Failed to find font weight bold" warning.
            if name in _FONT_PREFERENCE:
                font_manager.fontManager.addfont(font_path)
                available.add(name)

    for candidate in _FONT_PREFERENCE:
        if candidate in available:
            return candidate
    return "DejaVu Sans"


FONT_FAMILY = detect_font_family()


# ============================================================
# SIZE TARGETS
#
# ~182 mm (7.16 in) full-width figure, ~7 pt body text, ~8.5 pt bold
# panel letters, at final manuscript size.
# ============================================================

FULL_WIDTH_IN = 7.16

BASE_FONT_SIZE = 7.0
AXIS_LABEL_SIZE = 7.5
TITLE_SIZE = 8.0
TICK_LABEL_SIZE = 6.5
LEGEND_SIZE = 9.0
PANEL_LABEL_SIZE = 8.5


def apply_style() -> None:
    """Set the shared rcParams. Call once near the top of each script."""

    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": BASE_FONT_SIZE,
            "axes.labelsize": AXIS_LABEL_SIZE,
            "axes.titlesize": TITLE_SIZE,
            "xtick.labelsize": TICK_LABEL_SIZE,
            "ytick.labelsize": TICK_LABEL_SIZE,
            "legend.fontsize": LEGEND_SIZE,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": True,
        }
    )


# ============================================================
# AXIS / SPINE / GRID
# ============================================================


def clean_spines(axis: plt.Axes, keep=("left", "bottom")) -> None:
    for spine_name, spine in axis.spines.items():
        spine.set_visible(spine_name in keep)
    for spine_name in keep:
        axis.spines[spine_name].set_color("black")
        axis.spines[spine_name].set_linewidth(0.7)


def subtle_grid(axis: plt.Axes, axis_direction: str = "y") -> None:
    axis.grid(
        axis=axis_direction,
        color=GRID_COLOR,
        linewidth=0.45,
        linestyle="-",
        zorder=0,
    )
    axis.set_axisbelow(True)


def panel_letter(
    axis: plt.Axes,
    letter: str,
    x: float = -0.10,
    y: float = 1.05,
    fontsize: float = PANEL_LABEL_SIZE,
) -> None:
    axis.text(
        x,
        y,
        letter,
        transform=axis.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        fontstyle="normal",
        ha="left",
        va="bottom",
    )


# ============================================================
# LEGEND BUILDERS
# ============================================================


def sample_legend_handles(colors: dict[str, str] | None = None) -> list[Patch]:
    colors = colors or SAMPLE_COLORS
    return [
        Patch(facecolor=colors[sample], edgecolor="none", label=sample)
        for sample in SAMPLE_ORDER
    ]


def sample_marker_handles(marker: str = "o") -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="none",
            markerfacecolor=SAMPLE_COLORS[sample],
            markeredgecolor="black",
            markeredgewidth=0.4,
            markersize=5,
            label=sample,
        )
        for sample in SAMPLE_ORDER
    ]


def aligner_marker_handles(aligners: list[str] | None = None) -> list[Line2D]:
    aligners = aligners or ALIGNER_ORDER
    return [
        Line2D(
            [0],
            [0],
            marker=ALIGNER_MARKERS[aligner],
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.6,
            markersize=5,
            label=aligner,
        )
        for aligner in aligners
    ]


# ============================================================
# SAVE
# ============================================================

# Figures/tables this repo produces that the F2 thesis report includes,
# mapped to their destination path inside that sibling repo. Mirroring
# here means rerunning a plotting script after an edit also updates the
# report's copy immediately, instead of relying on a manual sync step.
REPORT_FIGURE_MAP = {
    "01_alignment_error_rate_30x.pdf": "figures/01_alignment_error_rate_30x.pdf",
    "04_mq0_reads_30x.pdf": "figures/04_mq0_reads_30x.pdf",
    "05_runtime_30x.pdf": "figures/05_runtime_30x.pdf",
    "06_cigar_composition_30x.pdf": "figures/06_cigar_composition_30x.pdf",
    "08_memory_30x.pdf": "figures/08_memory_30x.pdf",
    "11_mapped_unmapped_reads_30x.pdf": "figures/11_mapped_unmapped_reads_30x.pdf",
    "12_mapped_unmapped_bases_30x.pdf": "figures/12_mapped_unmapped_bases_30x.pdf",
    "13_cigar_yield_vs_error_30x.pdf": "figures/13_cigar_yield_vs_error_30x.pdf",
    "cigar_yield_error_diagnostics_combined.pdf": "figures/correlation_diagnostics/cigar_yield_error_diagnostics_combined.pdf",
    "cigar_yield_error_correlation_diagnostics.tsv": "tables/cigar_yield_error_correlation_diagnostics.tsv",
}


def mirror_to_report(output_path: Path) -> None:
    """Copy output_path into the f2-thesis-report repo, if applicable.

    No-op if output_path's filename isn't report-tracked, or if that repo
    (a sibling of lrs_benchmarking) isn't checked out on this machine.
    """
    dest_rel = REPORT_FIGURE_MAP.get(output_path.name)
    if dest_rel is None:
        return

    repo_root = next(
        (p for p in output_path.resolve().parents if (p / "alignment_analysis").is_dir()),
        None,
    )
    if repo_root is None:
        return

    report_root = repo_root.parent / "f2-thesis-report"
    if not report_root.is_dir():
        return

    dest = report_root / dest_rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(output_path.read_bytes())


def save_figure(figure: plt.Figure, output_png: Path, output_pdf: Path) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    mirror_to_report(output_pdf)
