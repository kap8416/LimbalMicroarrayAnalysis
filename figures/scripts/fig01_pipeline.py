"""
Figure 1 — Integrative systems-level pipeline.

Schematic of the analysis: dataset curation, cross-platform harmonization, the
three limbus-anchored contrasts, and the five analytical layers applied to each.

This is a conceptual figure and takes no data input; every count shown is a
design parameter defined in R/00_config.R and re-stated here so that the figure
and the configuration cannot drift apart silently.

Layout convention: stage labels occupy the left gutter (x < LEFT), all content
boxes start at x = LEFT. Coordinates are in an arbitrary 0-100 x 0-100 space.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from fig_style import GREY, LIGHTGREY, PAL_CONTRAST, PAL_GROUP, W2, apply_style, save

apply_style()

TEAL = PAL_GROUP["Limbus"]
ORANGE = "#C4622D"

LEFT = 15.0          # content left edge; the gutter to its left holds stage labels
RIGHT = 99.0
WIDTH = RIGHT - LEFT


def box(ax, x, y, w, h, *, title=None, body=None, fc="white", ec=GREY,
        title_fs=5.9, body_fs=5.4, lw=0.7, tc="black"):
    """Rounded box with an optional bold title above a lighter body block."""
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0,rounding_size=0.9",
            linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2,
        )
    )
    cx = x + w / 2
    if title and body:
        ax.text(cx, y + h * 0.70, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color=tc, zorder=3,
                linespacing=1.35)
        ax.text(cx, y + h * 0.29, body, ha="center", va="center",
                fontsize=body_fs, color=GREY, zorder=3, linespacing=1.4)
    else:
        ax.text(cx, y + h / 2, title or body, ha="center", va="center",
                fontsize=title_fs if title else body_fs,
                fontweight="bold" if title else "normal",
                color=tc, zorder=3, linespacing=1.4)


def arrow(ax, p0, p1, color=GREY, lw=0.8):
    ax.add_patch(
        FancyArrowPatch(
            p0, p1, arrowstyle="-|>", mutation_scale=6.5,
            linewidth=lw, color=color, zorder=1, shrinkA=1.0, shrinkB=1.0,
        )
    )


def main() -> None:
    fig, ax = plt.subplots(figsize=(W2, 7.1))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def stage(y, n, label):
        ax.text(1.0, y, f"{n}", fontsize=9.5, fontweight="bold", color=TEAL,
                ha="left", va="center")
        ax.text(4.6, y, label, fontsize=5.9, fontweight="bold", color=TEAL,
                ha="left", va="center", linespacing=1.3)

    # ================= 1. Data ============================================= #
    stage(92.5, "1", "DATA")
    hw = (WIDTH - 4) / 2
    box(ax, LEFT, 88, hw, 9,
        title="GSE38190 · GPL570",
        body="Affymetrix HG-U133 Plus 2.0\nlimbus 4 · cornea 3 · conjunctiva 3",
        fc="#F0F6F5", ec=TEAL)
    box(ax, LEFT + hw + 4, 88, hw, 9,
        title="GSE56421 · GPL6244",
        body="Affymetrix HuGene-1_0-st\ncornea 3 · cultured cells 3",
        fc="#F0F6F5", ec=TEAL)

    box(ax, LEFT, 80.5, WIDTH, 5.5,
        body="Curation: limbal fibroblasts (GSM1361191–93) and one conjunctival\n"
             "outlier (GSM724094, within-group $r$ = 0.92) excluded  →  "
             "16 samples · 4 groups",
        fc="white", ec=LIGHTGREY, body_fs=5.3)
    arrow(ax, (LEFT + hw / 2, 88), (LEFT + hw / 2, 86))
    arrow(ax, (LEFT + hw + 4 + hw / 2, 88), (LEFT + hw + 4 + hw / 2, 86))

    # ================= 2. Harmonization ==================================== #
    stage(72, "2", "BATCH\nCORRECTION")
    steps = [
        ("Quantile normalization", "within platform\nlimma::normalizeBetweenArrays"),
        ("Cornea-anchored shift", "the only group present\non both platforms"),
        ("Residual ComBat", "sva, mod = ~is_Cornea"),
    ]
    sw = (WIDTH - 2 * 3.5) / 3
    for i, (t, b) in enumerate(steps):
        x = LEFT + i * (sw + 3.5)
        box(ax, x, 67, sw, 9.5, title=t, body=b, fc="#FBF7F2", ec=ORANGE,
            title_fs=5.7, body_fs=5.0)
        if i:
            arrow(ax, (x - 3.5, 71.75), (x, 71.75), color=ORANGE)
    arrow(ax, (LEFT + WIDTH / 2, 80.5), (LEFT + WIDTH / 2, 76.5))

    box(ax, LEFT, 59.5, WIDTH, 6,
        body="45,378 RefSeq-level estimates  →  24,002 genes "
             "(highest-MAD isoform per symbol)\n"
             "QC gate: cross-platform corneal concordance $r$ = 0.953  "
             "(pre-specified $\\geq$ 0.95)",
        fc="white", ec=LIGHTGREY, body_fs=5.3)
    arrow(ax, (LEFT + WIDTH / 2, 67), (LEFT + WIDTH / 2, 65.5))

    # ================= 3. Contrasts ======================================== #
    stage(51, "3", "CONTRASTS")
    contrasts = [
        ("C1", "Limbus vs cornea", "|log$_2$FC| $\\geq$ 1.0\ncross-platform"),
        ("C2", "Limbus vs conjunctiva", "|log$_2$FC| $\\geq$ 0.5\nsame platform · most reliable"),
        ("C3", "Limbus vs cultured cells", "|log$_2$FC| $\\geq$ 1.0\nconfounded · exploratory"),
    ]
    cw = (WIDTH - 2 * 3.5) / 3
    cxs = []
    for i, (cid, title, note) in enumerate(contrasts):
        x = LEFT + i * (cw + 3.5)
        cx = x + cw / 2
        cxs.append(cx)
        col = PAL_CONTRAST[cid]
        box(ax, x, 45.5, cw, 10.5, title=f"{cid}\n{title}", body=note,
            fc="white", ec=col, lw=1.0, title_fs=5.9, body_fs=5.0, tc=col)
        arrow(ax, (cx, 59.5), (cx, 56), color=col)

    ax.text(RIGHT, 44.3,
            "limma · lmFit on a no-intercept design · eBayes(trend = TRUE) · "
            "FDR $\\leq$ 0.05 (Benjamini–Hochberg)",
            fontsize=5.0, ha="right", va="top", color=GREY, style="italic")

    # ================= 4. Analysis layers ================================== #
    stage(30, "4", "ANALYSIS\nLAYERS")

    # A single bus carries all three contrasts into every analytical layer,
    # rather than implying a one-to-one contrast-to-layer mapping.
    bus_y = 40.0
    for cx, (cid, _, _) in zip(cxs, contrasts):
        ax.plot([cx, cx], [45.5, bus_y], color=PAL_CONTRAST[cid], lw=0.8, zorder=1)
    ax.plot([LEFT + 2, RIGHT - 2], [bus_y, bus_y], color=GREY, lw=0.9, zorder=1)

    layers = [
        ("Differential\nexpression", "limma", ORANGE),
        ("Functional\nenrichment", "clusterProfiler\nGO-BP · KEGG", "#3D5A8A"),
        ("PPI network\ntopology", "STRING v12.0\nigraph · Louvain", "#8A6BA8"),
        ("Master regulator\nanalysis", "corto\n1,000 bootstraps", TEAL),
        ("Regulatory\nrewiring", "|$\\Delta\\rho$| $\\geq$ 0.6\nSpearman", "#7A6A3A"),
    ]
    gap = 2.0
    lw_ = (WIDTH - 4 * gap) / 5
    for i, (t, tool, col) in enumerate(layers):
        x = LEFT + i * (lw_ + gap)
        cx = x + lw_ / 2
        box(ax, x, 24.5, lw_, 11, title=t, body=tool, fc="white", ec=col,
            title_fs=5.5, body_fs=4.8)
        arrow(ax, (cx, bus_y), (cx, 35.5), color=GREY, lw=0.7)

    # ================= 5. Outputs ========================================== #
    stage(12, "5", "FINDINGS")
    ow = (WIDTH - 4) / 2
    box(ax, LEFT, 6.5, ow, 12,
        title="Conserved limbal identity",
        body="32-gene core signature\nFN1-anchored ECM interactome\n"
             "ECM organization and cell-substrate adhesion\nas dominant functional features",
        fc="#F0F6F5", ec=TEAL, lw=1.0, title_fs=5.9, body_fs=5.1)
    box(ax, LEFT + ow + 4, 6.5, ow, 12,
        title="Regulatory architecture",
        body="Master regulators ZEB1/2, ETS1, MEF2C\n"
             "culture-associated HMGA2, VIM, NF1\n"
             "module-specific rather than global rewiring\nacross limbal transitions",
        fc="#F0F6F5", ec=TEAL, lw=1.0, title_fs=5.9, body_fs=5.1)
    arrow(ax, (LEFT + ow / 2, 24.5), (LEFT + ow / 2, 18.5), color=TEAL)
    arrow(ax, (LEFT + ow + 4 + ow / 2, 24.5), (LEFT + ow + 4 + ow / 2, 18.5), color=TEAL)

    box(ax, LEFT, 0.8, WIDTH, 4.2,
        body="Code and intermediate results: github.com/kap8416/LimbalMicroarrayAnalysis",
        fc="#FAFAFA", ec=LIGHTGREY, body_fs=5.3, tc=GREY)

    save(fig, "Fig01_pipeline")


if __name__ == "__main__":
    main()
