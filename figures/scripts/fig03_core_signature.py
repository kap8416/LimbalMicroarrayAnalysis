"""
Figure 3 — Contrast-specific and shared limbal programs converge on a 32-gene core.

(a) UpSet-style intersection plot of genes upregulated in limbus across C1-C3.
(b) log2 fold change of the 32-gene core signature in all three contrasts.
(c) Functional annotation of the core signature.

Inputs
  results/de/deg_membership.csv

The core signature is defined as genes upregulated in limbus in ALL THREE
contrasts. Predicted-only RefSeq accessions (XM_/XR_) meet the same statistical
criterion but carry no approved symbol; they are counted separately in panel (a)
and excluded from panel (b), which is why the core is 32 named genes out of 49
significant features.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from fig_style import (
    DIR_DE,
    GREY,
    LIGHTGREY,
    PAL_CONTRAST,
    W2,
    apply_style,
    is_predicted,
    load_membership,
    panel_label,
    save,
)

apply_style()

CONTRASTS = ["C1", "C2", "C3"]

# Functional assignment of the core signature, used for panel (c) row grouping.
# Categories follow the interpretation given in the Discussion.
CORE_FUNCTION = {
    "ECM structure and regulation": [
        "COL8A1", "COL8A2", "THBS4", "ABI3BP", "FIBIN", "ITGBL1", "ITGA10", "CDH11",
    ],
    "Developmental / progenitor signalling": [
        "BMP4", "FZD7", "FGF18", "ANGPT1", "NTRK3", "GRP", "PTGFR", "SHC3",
    ],
    "Transcriptional regulation": [
        "ZEB1", "EYA4", "POU3F3", "SHOX", "ADARB1", "PRR16",
    ],
    "Metabolic and oxidative protection": [
        "SOD3", "CYTL1", "PLPP3", "HACD1",
    ],
    "Other / uncharacterized": [
        "C14orf132", "CCDC102B", "CEP126", "JCAD", "KCNE4", "LRRC3B",
    ],
}


def core_table(mm: pd.DataFrame) -> pd.DataFrame:
    """Named genes upregulated in limbus in all three contrasts."""
    up = np.logical_and.reduce([mm[f"{c}_status"] == "Up" for c in CONTRASTS])
    core = mm[up].copy()
    core = core[~is_predicted(core["Symbol"])]
    return core.sort_values("C2_logFC", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Panel a — intersections of the upregulated sets
# --------------------------------------------------------------------------- #

def panel_upset(ax_bar, ax_mat, mm: pd.DataFrame) -> None:
    up = {c: set(mm.loc[mm[f"{c}_status"] == "Up", "Symbol"]) for c in CONTRASTS}

    combos = [
        (("C1",), "C1 only"),
        (("C3",), "C3 only"),
        (("C1", "C3"), "C1 & C3"),
        (("C2",), "C2 only"),
        (("C2", "C3"), "C2 & C3"),
        (("C1", "C2"), "C1 & C2"),
        (("C1", "C2", "C3"), "core"),
    ]

    rows = []
    for members, _ in combos:
        inc = set.intersection(*[up[c] for c in members])
        exc = set().union(*[up[c] for c in CONTRASTS if c not in members]) if len(members) < 3 else set()
        genes = inc - exc
        named = {g for g in genes if not str(g).startswith(("XM_", "XR_", "NR_"))}
        rows.append({"members": members, "n": len(genes), "n_named": len(named)})

    df = pd.DataFrame(rows).sort_values("n", ascending=False).reset_index(drop=True)
    x = np.arange(len(df))

    # Bars: total significant features, with the named-gene fraction overlaid.
    ax_bar.bar(x, df["n"], 0.62, color=LIGHTGREY, label="all significant features")
    ax_bar.bar(x, df["n_named"], 0.62, color="#1B7A6E", label="named genes")

    for i, r in df.iterrows():
        ax_bar.annotate(
            f"{r['n_named']}" if r["n_named"] == r["n"] else f"{r['n_named']}/{r['n']}",
            (i, r["n"]),
            ha="center", va="bottom", fontsize=6, color=GREY,
            xytext=(0, 1.5), textcoords="offset points",
        )

    ax_bar.set_ylabel("Genes upregulated\nin limbus")
    ax_bar.set_xlim(-0.6, len(df) - 0.4)
    ax_bar.set_ylim(0, df["n"].max() * 1.16)
    ax_bar.set_xticks([])
    ax_bar.yaxis.grid(True)
    ax_bar.legend(loc="upper right", fontsize=6)
    ax_bar.set_title("Overlap of limbal-upregulated gene sets", pad=5)

    # Membership matrix beneath the bars.
    for yi, c in enumerate(CONTRASTS):
        y = len(CONTRASTS) - 1 - yi
        ax_mat.scatter(x, [y] * len(x), s=26, color="#E4E4E4", zorder=1, linewidths=0)
        sel = [i for i, r in df.iterrows() if c in r["members"]]
        ax_mat.scatter(sel, [y] * len(sel), s=26, color=PAL_CONTRAST[c], zorder=3, linewidths=0)

    for i, r in df.iterrows():
        ys = [len(CONTRASTS) - 1 - CONTRASTS.index(c) for c in r["members"]]
        if len(ys) > 1:
            ax_mat.plot([i, i], [min(ys), max(ys)], color=GREY, lw=1.0, zorder=2)

    # Mark the three-way intersection directly on its bar.
    core_i = int(df.index[df["members"].apply(len) == 3][0])
    ax_bar.annotate(
        "core\nsignature",
        (core_i, df.loc[core_i, "n"]),
        ha="center", va="bottom", fontsize=6, style="italic", color="#1B7A6E",
        xytext=(0, 12), textcoords="offset points",
        arrowprops=dict(arrowstyle="-", lw=0.6, color="#1B7A6E",
                        shrinkA=0.5, shrinkB=1.5),
    )

    ax_mat.set_xlim(-0.6, len(df) - 0.4)
    ax_mat.set_ylim(-0.6, len(CONTRASTS) - 0.4)
    ax_mat.set_yticks(range(len(CONTRASTS)))
    ax_mat.set_yticklabels(CONTRASTS[::-1])
    ax_mat.set_xticks([])
    for s in ax_mat.spines.values():
        s.set_visible(False)
    ax_mat.tick_params(length=0)


# --------------------------------------------------------------------------- #
# Panel b — core signature effect sizes
# --------------------------------------------------------------------------- #

def panel_core_heatmap(ax, core: pd.DataFrame) -> None:
    # Order rows by functional category, then by mean effect size.
    fn_of = {g: k for k, v in CORE_FUNCTION.items() for g in v}
    core = core.assign(
        func=core["Symbol"].map(fn_of).fillna("Other / uncharacterized"),
        mean_lfc=core[[f"{c}_logFC" for c in CONTRASTS]].mean(axis=1),
    )
    cat_order = list(CORE_FUNCTION)
    core["func"] = pd.Categorical(core["func"], categories=cat_order, ordered=True)
    core = core.sort_values(["func", "mean_lfc"], ascending=[True, False]).reset_index(drop=True)

    mat = core[[f"{c}_logFC" for c in CONTRASTS]].to_numpy()
    norm = TwoSlopeNorm(vmin=0, vcenter=2.0, vmax=max(5.5, np.nanmax(mat)))

    im = ax.imshow(mat, cmap="YlGnBu", norm=norm, aspect="auto")

    ax.set_xticks(range(len(CONTRASTS)))
    ax.set_xticklabels(CONTRASTS)
    ax.set_yticks(range(len(core)))
    ax.set_yticklabels(core["Symbol"], fontsize=5.6, fontstyle="italic")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # Cell values.
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(
                j, i, f"{v:.1f}",
                ha="center", va="center", fontsize=4.9,
                color="white" if v > 3.2 else "#333333",
            )

    # Functional category brackets on the right.
    y0 = 0
    for cat in cat_order:
        n = int((core["func"] == cat).sum())
        if n == 0:
            continue
        ax.plot([2.62, 2.62], [y0 - 0.4, y0 + n - 0.6], color=GREY, lw=0.8,
                clip_on=False)
        ax.text(
            2.78, y0 + (n - 1) / 2,
            cat.replace(" and ", " &\n").replace(" / ", " /\n"),
            fontsize=5.4, va="center", ha="left", color=GREY, clip_on=False,
        )
        y0 += n

    ax.set_title("32-gene core limbal signature", pad=5)

    cb = ax.figure.colorbar(im, ax=ax, fraction=0.05, pad=0.95, aspect=16)
    cb.set_label("log$_2$FC vs comparator", fontsize=6)
    cb.ax.tick_params(labelsize=5.5, length=1.8)
    cb.outline.set_visible(False)


# --------------------------------------------------------------------------- #

def main() -> None:
    mm = load_membership()
    core = core_table(mm)
    assert len(core) == 32, f"Expected 32 named core genes, found {len(core)}"

    fig = plt.figure(figsize=(W2, 6.1))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.30, 1.0], height_ratios=[1.0, 0.24],
        hspace=0.06, wspace=0.52,
        left=0.075, right=0.735, top=0.945, bottom=0.055,
    )

    ax_bar = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bar)
    panel_upset(ax_bar, ax_mat, mm)
    panel_label(ax_bar, "a", dx=-0.150, dy=1.135)

    ax_b = fig.add_subplot(gs[:, 1])
    panel_core_heatmap(ax_b, core)
    panel_label(ax_b, "b", dx=-0.66, dy=1.055)

    save(fig, "Fig03_core_signature")


if __name__ == "__main__":
    main()
