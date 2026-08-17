"""
Figure 7 — Limbal identity as a niche-dependent regulatory attractor.

(a) Conceptual landscape: the limbal state occupies a deep basin stabilized by
    niche-derived ECM signalling; removal from the niche flattens the basin and
    the state drifts toward a culture-adapted configuration, while corneal
    differentiation is a separate, committed basin.
(b) Layered architecture of the attractor: the regulatory backbone, the ECM
    module that reinforces it, and the transitions that destabilize it.

Conceptual figure; no data input. Gene and regulator names shown are those
supported by the analyses in Figs. 2-6 and 8.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from fig_style import GREY, LIGHTGREY, PAL_GROUP, W2, apply_style, panel_label, save

apply_style()

TEAL = PAL_GROUP["Limbus"]
ORANGE = PAL_GROUP["Cornea"]
VIOLET = PAL_GROUP["Cultured"]
BLUE = PAL_GROUP["Conjunctiva"]


# --------------------------------------------------------------------------- #
# Panel a — attractor landscape
# --------------------------------------------------------------------------- #

def basin(x, centre, depth, width):
    """A single Gaussian well; summed wells give the potential landscape."""
    return -depth * np.exp(-((x - centre) ** 2) / (2 * width**2))


def panel_landscape(ax) -> None:
    x = np.linspace(0, 10, 1200)

    # Native niche: the limbal basin is deep and narrow, flanked by the
    # conjunctival and corneal basins.
    native = (
        basin(x, 2.0, 0.62, 0.72)      # conjunctiva
        + basin(x, 5.0, 1.00, 0.60)    # limbus
        + basin(x, 8.2, 0.78, 0.85)    # cornea
    )
    # Ex vivo: loss of niche ECM input shallows the limbal well and deepens a
    # culture-adapted well; the attractor is eroded rather than replaced.
    exvivo = (
        basin(x, 2.0, 0.55, 0.72)
        + basin(x, 5.0, 0.26, 0.78)    # limbal well strongly shallowed
        + basin(x, 8.2, 0.70, 0.85)
        + basin(x, 6.6, 0.80, 0.62)    # culture-adapted state, now the deeper well
    )

    ax.plot(x, native, color=TEAL, lw=1.6, zorder=3, label="Native niche")
    ax.plot(x, exvivo, color=VIOLET, lw=1.3, ls=(0, (4, 2)), zorder=3,
            label="Ex vivo culture")
    ax.fill_between(x, native, native.min() - 0.30, color=TEAL, alpha=0.07, zorder=1)

    # State markers, placed on the native curve.
    for cx, col, lab in [
        (2.0, BLUE, "Conjunctiva"),
        (5.0, TEAL, "Limbus"),
        (8.2, ORANGE, "Cornea"),
    ]:
        y = native[np.argmin(np.abs(x - cx))]
        ax.scatter([cx], [y], s=52, color=col, zorder=5, edgecolors="white",
                   linewidths=0.9)
        ax.annotate(lab, (cx, y), fontsize=6.2, fontweight="bold", color=col,
                    ha="center", va="top", xytext=(0, -9), textcoords="offset points")

    y_cult = exvivo[np.argmin(np.abs(x - 6.6))]
    ax.scatter([6.6], [y_cult], s=40, color=VIOLET, zorder=5, edgecolors="white",
               linewidths=0.9)
    ax.annotate("Culture-\nadapted", (6.6, y_cult), fontsize=5.8, color=VIOLET,
                ha="center", va="bottom", xytext=(0, 8), textcoords="offset points")

    # Destabilization arrow.
    y_limb = native[np.argmin(np.abs(x - 5.0))]
    ax.add_patch(
        FancyArrowPatch(
            (5.22, y_limb + 0.10), (6.42, y_cult + 0.06),
            arrowstyle="-|>", mutation_scale=7, lw=1.1, color=VIOLET,
            connectionstyle="arc3,rad=-0.42", zorder=6,
        )
    )
    ax.text(5.85, y_limb + 0.44, "niche removal\n(ex vivo expansion)", fontsize=5.4,
            color=VIOLET, ha="center", va="bottom", style="italic", zorder=6,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

    # Stabilizing input into the limbal well.
    ax.add_patch(
        FancyArrowPatch(
            (3.55, y_limb + 0.34), (4.80, y_limb + 0.04),
            arrowstyle="-|>", mutation_scale=7, lw=1.1, color=TEAL, zorder=6,
            connectionstyle="arc3,rad=0.22",
        )
    )
    ax.text(3.45, y_limb + 0.38, "FN1-centred ECM\nniche signalling", fontsize=5.6,
            color=TEAL, ha="right", va="bottom", fontweight="bold", zorder=6,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

    ax.set_xlim(0.3, 9.9)
    ax.set_ylim(native.min() - 0.34, 0.62)
    ax.set_xlabel("Transcriptional state space")
    ax.set_ylabel("Regulatory coherence\n(attractor depth)")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.legend(loc="upper right", fontsize=5.8)
    ax.set_title("Limbal identity as a regulatory attractor", pad=5)


# --------------------------------------------------------------------------- #
# Panel b — layered architecture
# --------------------------------------------------------------------------- #

def layer(ax, y, h, label, body, fc, ec, label_fs=6.0):
    ax.add_patch(
        FancyBboxPatch(
            (14, y), 84, h,
            boxstyle="round,pad=0,rounding_size=1.1",
            facecolor=fc, edgecolor=ec, lw=0.9, zorder=2,
        )
    )
    ax.text(11.5, y + h / 2, label, fontsize=label_fs, fontweight="bold",
            color=ec, ha="right", va="center", linespacing=1.25)
    ax.text(56, y + h / 2, body, fontsize=5.4, color="black",
            ha="center", va="center", linespacing=1.5, zorder=3)


def panel_architecture(ax) -> None:
    ax.set_xlim(0, 104)
    ax.set_ylim(0, 100)
    ax.axis("off")

    layer(
        ax, 74, 20,
        "NICHE\nINPUT",
        "Limbal stromal and vascular microenvironment\n"
        "$\\bf{FN1}$ · COL8A1/2 · THBS4 · LAMB3 basement membrane\n"
        "paracrine IL6 · CXCL12 · FGF2 · BMP4",
        "#F0F6F5", TEAL,
    )
    layer(
        ax, 50, 20,
        "REGULATORY\nBACKBONE",
        "Master regulators: $\\bf{ZEB1}$ · $\\bf{ZEB2}$ · $\\bf{ETS1}$ · MEF2C · THRB · PPARG\n"
        "PAX6 constitutively expressed, permissive rather than instructive\n"
        "repressed regulon activity: PTCH1 · OVOL1 · KLF4",
        "#FBF7F2", "#C4622D",
    )
    layer(
        ax, 26, 20,
        "CORE\nSIGNATURE",
        "32 genes upregulated in limbus across all three contrasts\n"
        "ECM (CDH11, COL8A1/2, THBS4, ITGBL1, ITGA10, ABI3BP, FIBIN)\n"
        "signalling (BMP4, FZD7, FGF18, ANGPT1, NTRK3) · protection (SOD3, CYTL1)",
        "#F4F1F7", "#8A6BA8",
    )
    layer(
        ax, 2, 20,
        "DESTABILIZ-\nATION",
        "Ex vivo expansion → integrin-cytoskeletal rewiring (ITGA2/5/6, ITGB4/6, PTK2, MMP9)\n"
        "HMGA2 · VIM · NF1 regulon activity rises; ECM-instructive input is lost\n"
        "module-specific rewiring (MEF2C, MITF, ZEB2, IKZF1, FLI1), not global collapse",
        "#FAFAFA", GREY,
    )

    # Reinforcement / erosion arrows between layers.
    for y0, y1, col, lab in [
        (74, 70, TEAL, "instructs"),
        (50, 46, "#C4622D", "maintains"),
    ]:
        ax.add_patch(
            FancyArrowPatch((56, y0), (56, y1), arrowstyle="-|>",
                            mutation_scale=7, lw=1.0, color=col, zorder=1)
        )
        ax.text(58, (y0 + y1) / 2, lab, fontsize=5.2, color=col,
                ha="left", va="center", style="italic")

    ax.add_patch(
        FancyArrowPatch((56, 26), (56, 22), arrowstyle="-|>",
                        mutation_scale=7, lw=1.0, color=GREY, zorder=1)
    )
    ax.text(58, 24, "erodes", fontsize=5.2, color=GREY, ha="left", va="center",
            style="italic")

    # Positive feedback loop from the core signature back to the niche.
    ax.add_patch(
        FancyArrowPatch((98, 36), (98, 84), arrowstyle="-|>", mutation_scale=7,
                        lw=1.0, color=TEAL, connectionstyle="arc3,rad=-0.55",
                        zorder=1)
    )
    ax.text(100.5, 60, "self-reinforcing loop", fontsize=5.2, color=TEAL,
            ha="center", va="center", style="italic", rotation=90)

    ax.set_title("Layered architecture of the attractor", pad=5)


# --------------------------------------------------------------------------- #

def main() -> None:
    fig = plt.figure(figsize=(W2, 6.3))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[1.0, 1.22], hspace=0.30,
        left=0.115, right=0.975, top=0.935, bottom=0.055,
    )

    ax_a = fig.add_subplot(gs[0])
    panel_landscape(ax_a)
    panel_label(ax_a, "a", dx=-0.105, dy=1.10)

    ax_b = fig.add_subplot(gs[1])
    panel_architecture(ax_b)
    panel_label(ax_b, "b", dx=-0.105, dy=1.055)

    save(fig, "Fig07_attractor_model")


if __name__ == "__main__":
    main()
