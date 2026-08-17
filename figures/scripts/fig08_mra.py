"""
Figure 8 — Master regulator architecture and regulatory rewiring.

(a-c) Normalized enrichment score (NES) of the most significant regulators in
      each contrast, ordered by |NES|.
(d)   Regulators ranked by the number of rewired co-expression edges
      (|Delta rho| >= 0.6 between conditions).

NES quantifies REGULON ACTIVITY, not transcript abundance. PAX6 illustrates why
the distinction matters: its transcript is constitutively expressed and
non-differential in every contrast, yet its regulon activity is significantly
reduced in limbus relative to cornea. That distinction belongs in the caption,
not inside the figure: explanatory boxes drawn on the axes consume panel area and
duplicate text the reader already has beneath the figure.

Inputs
  results/mra/MRA_all.csv                produced by R/06_mra_rewiring.R
  results/mra/rewiring_tf_summary.csv    produced by R/06_mra_rewiring.R
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import (
    DIR_MRA,
    GREY,
    PAL_CONTRAST,
    PAL_DIR,
    CONTRAST_LABEL,
    W2,
    apply_style,
    panel_label,
    require,
    save,
)

apply_style()

CONTRASTS = ["C1", "C2", "C3"]
N_TF = 10  # regulators shown per contrast (top |NES|)
N_REWIRED = 18  # regulators shown in the rewiring panel


def panel_nes(ax, mra: pd.DataFrame, contrast: str) -> None:
    sub = mra[mra["contrast"] == contrast].copy()
    if "padj" in sub.columns:
        sub = sub[sub["padj"] <= 0.05]
    sub = sub.reindex(sub["NES"].abs().sort_values(ascending=False).index).head(N_TF)
    sub = sub.sort_values("NES")

    y = np.arange(len(sub))
    colors = [PAL_DIR["Up"] if v > 0 else PAL_DIR["Down"] for v in sub["NES"]]
    ax.barh(y, sub["NES"], 0.66, color=colors)
    ax.axvline(0, color=GREY, lw=0.7)

    for yi, (_, r) in zip(y, sub.iterrows()):
        right = r["NES"] > 0
        ax.annotate(
            f"{r['NES']:+.2f}",
            (r["NES"], yi),
            va="center", ha="left" if right else "right",
            fontsize=5.2, color=GREY,
            xytext=(2.5 if right else -2.5, 0), textcoords="offset points",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(sub["tf"], fontstyle="italic", fontsize=6.2)
    # Extra margin so the numeric labels beyond the bar ends stay inside the frame.
    lim = max(abs(sub["NES"])) * 1.48
    ax.set_xlim(-lim, lim)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True)
    ax.set_xlabel("NES (regulon activity in limbus)")
    ax.set_title(f"{contrast} · {CONTRAST_LABEL[contrast]}",
                 color=PAL_CONTRAST[contrast], pad=4, fontsize=7.2)

    if contrast == "C3":
        ax.text(0.98, 0.03, "exploratory", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5.2, style="italic", color=GREY)


def panel_rewiring(ax) -> None:
    rw = pd.read_csv(require(DIR_MRA / "rewiring_tf_summary.csv", "rewiring summary"))
    col = "n_rewired_edges" if "n_rewired_edges" in rw.columns else rw.columns[1]
    d = rw.nlargest(N_REWIRED, col).sort_values(col)

    y = np.arange(len(d))
    ax.barh(y, d[col], 0.66, color="#7A6A3A")
    for yi, v in zip(y, d[col]):
        ax.annotate(f"{int(v)}", (v, yi), va="center", ha="left",
                    fontsize=5.4, color=GREY, xytext=(2.5, 0),
                    textcoords="offset points")

    ax.set_yticks(y)
    ax.set_yticklabels(d["tf"], fontstyle="italic", fontsize=6.2)
    ax.set_xlabel("Rewired co-expression edges (|$\\Delta\\rho$| $\\geq$ 0.6)")
    ax.set_xlim(0, d[col].max() * 1.16)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True)
    ax.set_title("Context-dependent edge turnover", pad=4)


def main() -> None:
    mra = pd.read_csv(require(DIR_MRA / "MRA_all.csv", "master regulator results"))

    fig = plt.figure(figsize=(W2, 5.8))
    gs = fig.add_gridspec(
        2, 3, height_ratios=[1.0, 1.05], hspace=0.52, wspace=0.62,
        left=0.085, right=0.98, top=0.925, bottom=0.11,
    )

    for j, (c, letter) in enumerate(zip(CONTRASTS, "abc")):
        ax = fig.add_subplot(gs[0, j])
        panel_nes(ax, mra, c)
        panel_label(ax, letter, dx=-0.42, dy=1.16)

    ax_d = fig.add_subplot(gs[1, :])
    panel_rewiring(ax_d)
    panel_label(ax_d, "d", dx=-0.075, dy=1.09)

    save(fig, "Fig08_master_regulators")


if __name__ == "__main__":
    main()
