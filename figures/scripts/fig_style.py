"""
Shared figure style for the LimbalMicroarrayAnalysis manuscript figures.

Aesthetic target: Nature Reviews / Molecular Plant conventions — thin axes, no
top/right spines, sans-serif type at 7-9 pt, restrained categorical palette,
white background, 600 dpi raster export alongside vector PDF.

Every figure script imports from here so that panel letters, colours and type
sizes are identical across the figure set.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

REPO = Path(__file__).resolve().parents[2]
DIR_DE = REPO / "results" / "de"
DIR_ENR = REPO / "results" / "enrichment"
DIR_PPI = REPO / "results" / "ppi"
DIR_MRA = REPO / "results" / "mra"
DIR_PROC = REPO / "data" / "processed"
DIR_OUT = REPO / "figures" / "output"
DIR_OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# Palette — kept in sync with PAL_* in R/00_config.R
# --------------------------------------------------------------------------- #

PAL_GROUP = {
    "Limbus": "#1B7A6E",
    "Cornea": "#C4622D",
    "Conjunctiva": "#3D5A8A",
    "Cultured": "#8A6BA8",
}

PAL_CONTRAST = {"C1": "#C4622D", "C2": "#3D5A8A", "C3": "#8A6BA8"}

PAL_DIR = {"Up": "#B03A2E", "Down": "#2471A3", "NS": "#C9C9C9"}

# Sequential ramp for module / enrichment shading.
PAL_SEQ = ["#F2F6F5", "#CFE3E0", "#9AC7C0", "#5FA69C", "#2E8378", "#1B7A6E"]

GREY = "#4D4D4D"
LIGHTGREY = "#BFBFBF"

CONTRAST_LABEL = {
    "C1": "Limbus vs cornea",
    "C2": "Limbus vs conjunctiva",
    "C3": "Limbus vs cultured cells",
}

# Fold-change thresholds used per contrast (see R/00_config.R).
LFC_CUT = {"C1": 1.0, "C2": 0.5, "C3": 1.0}
FDR_CUT = 0.05

# --------------------------------------------------------------------------- #
# rcParams
# --------------------------------------------------------------------------- #


def apply_style() -> None:
    """Install the shared rcParams. Call once at the top of each figure script."""
    mpl.rcParams.update(
        {
            # --- type ---
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 7.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "figure.titlesize": 9.5,
            # --- axes ---
            "axes.linewidth": 0.6,
            "axes.edgecolor": GREY,
            "axes.labelcolor": "black",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            # --- ticks ---
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.color": GREY,
            "ytick.color": GREY,
            "xtick.direction": "out",
            "ytick.direction": "out",
            # --- grid ---
            "grid.color": "#E8E8E8",
            "grid.linewidth": 0.5,
            # --- legend ---
            "legend.frameon": False,
            "legend.handlelength": 1.2,
            "legend.handletextpad": 0.5,
            "legend.columnspacing": 1.0,
            "legend.borderaxespad": 0.3,
            # --- figure / export ---
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,  # embed as TrueType, editable in Illustrator
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "lines.linewidth": 0.9,
            "patch.linewidth": 0.5,
        }
    )


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

# Journal single / 1.5 / double column widths in inches.
W1, W15, W2 = 3.35, 4.75, 6.85


def panel_label(ax, letter: str, dx: float = -0.16, dy: float = 1.05) -> None:
    """Place a bold lower-case panel letter in axes coordinates."""
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save(fig, stem: str, formats=("pdf", "png")) -> None:
    """Write a figure to figures/output/ in every requested format."""
    for ext in formats:
        path = DIR_OUT / f"{stem}.{ext}"
        fig.savefig(path)
        print(f"  wrote {path.relative_to(REPO)}")
    plt.close(fig)


def require(path: Path, what: str) -> Path:
    """Fail with an actionable message when an input table is missing."""
    if not Path(path).exists():
        raise SystemExit(
            f"MISSING INPUT: {what}\n"
            f"  expected at: {path}\n"
            f"  generate it by running the corresponding R/ step, or place the "
            f"existing file at that path."
        )
    return Path(path)


def load_membership() -> pd.DataFrame:
    """Cross-contrast DEG membership table (one row per gene DE in >= 1 contrast)."""
    return pd.read_csv(require(DIR_DE / "deg_membership.csv", "DEG membership table"))


def load_centrality(contrast: str) -> pd.DataFrame:
    """Per-node PPI topology for one contrast."""
    return pd.read_csv(
        require(DIR_PPI / f"centrality_{contrast}.csv", f"{contrast} PPI centrality table")
    )


def is_predicted(symbols: pd.Series) -> pd.Series:
    """Flag predicted-only RefSeq accessions, which carry no approved symbol."""
    return symbols.astype(str).str.match(r"^(XM_|XR_|NR_|NM_)")


def italicize(labels):
    """Render gene symbols in italics, as required for gene names."""
    return [f"$\\it{{{str(g).replace('_', chr(92) + '_')}}}$" for g in labels]
