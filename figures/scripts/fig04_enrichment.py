"""
Figure 4 — Functional reconfiguration of limbal identity across epithelial contexts.

Dot plot of GO-BP and KEGG over-representation results: one column per
contrast x direction, terms on the y axis, dot size = gene count, colour =
-log10 FDR.

Terms whose GO names invite over-reading are marked with a dagger and explained
in the legend, so the annotation rather than the label drives interpretation:
  - "ossification" (C2 up)              : shared ECM and BMP-pathway genes, not
                                          osteogenic differentiation.
  - "epidermis development" (C1/C3 down): keratinization programs shared by
                                          stratified epithelia, not epidermal
                                          specification.

Input
  results/enrichment/enrichment_all.csv     produced by R/04_functional_enrichment.R
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import (
    DIR_ENR,
    GREY,
    LIGHTGREY,
    CONTRAST_LABEL,
    W2,
    apply_style,
    panel_label,
    require,
    save,
)

apply_style()

CONTRASTS = ["C1", "C2", "C3"]
N_TERMS = 8  # terms retained per contrast x direction

CAVEAT_TERMS = {
    "ossification":
        "shared ECM / BMP-pathway genes, not osteogenic differentiation",
    "epidermis development":
        "keratinization shared by stratified epithelia, not epidermal specification",
    "keratinocyte differentiation":
        "keratinization shared by stratified epithelia, not epidermal specification",
}


def load() -> pd.DataFrame:
    df = pd.read_csv(require(DIR_ENR / "enrichment_all.csv", "enrichment results"))
    need = {"contrast", "direction", "source", "term", "padj", "count"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"enrichment_all.csv is missing columns: {sorted(missing)}")
    return df


def select_terms(df: pd.DataFrame) -> pd.DataFrame:
    """Top N_TERMS terms by FDR within each contrast x direction x source."""
    return (
        df.sort_values("padj")
        .groupby(["contrast", "direction", "source"], as_index=False)
        .head(N_TERMS)
    )


def panel_dotplot(ax, df: pd.DataFrame, source: str) -> None:
    sub = select_terms(df[df["source"] == source]).copy()
    if sub.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No significant {source} terms", ha="center", va="center",
                fontsize=6.5, color=GREY)
        return

    sub["col"] = sub["contrast"] + " " + sub["direction"]
    col_order = [f"{c} {d}" for c in CONTRASTS for d in ("Up", "Down")]
    col_order = [c for c in col_order if c in set(sub["col"])]

    # Order terms by the contrast in which each is most significant, so related
    # programs sit together rather than interleaving.
    key = sub.groupby("term").agg(best=("padj", "min"),
                                  first_col=("col", lambda s: col_order.index(s.iloc[0])))
    term_order = key.sort_values(["first_col", "best"], ascending=[True, False]).index.tolist()

    sub["xi"] = sub["col"].map({c: i for i, c in enumerate(col_order)})
    sub["yi"] = sub["term"].map({t: i for i, t in enumerate(term_order)})
    sub["nlp"] = -np.log10(sub["padj"].clip(lower=1e-300))

    sizes = 8 + 52 * (sub["count"] / sub["count"].max())
    sc = ax.scatter(
        sub["xi"], sub["yi"], s=sizes, c=sub["nlp"], cmap="YlOrRd",
        edgecolors="white", linewidths=0.4, zorder=3,
    )

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels([c.replace(" ", "\n") for c in col_order], fontsize=6)
    ax.set_yticks(range(len(term_order)))
    ax.set_yticklabels(
        [t + (" †" if t.lower() in CAVEAT_TERMS else "") for t in term_order],
        fontsize=5.8,
    )
    ax.set_xlim(-0.6, len(col_order) - 0.4)
    ax.set_ylim(-0.7, len(term_order) - 0.3)
    ax.grid(True, axis="both", color="#F0F0F0", lw=0.5, zorder=1)
    ax.tick_params(length=0)
    ax.set_title(f"{source.replace('_', '-')} over-representation", pad=5)

    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.035, pad=0.02, aspect=18)
    cb.set_label("-log$_{10}$ FDR", fontsize=6)
    cb.ax.tick_params(labelsize=5.5, length=1.8)
    cb.outline.set_visible(False)

    # Size legend.
    for i, n in enumerate([sub["count"].min(), sub["count"].median(), sub["count"].max()]):
        ax.scatter([], [], s=8 + 52 * (n / sub["count"].max()), c=LIGHTGREY,
                   label=f"{int(n)} genes")
    ax.legend(loc="lower right", fontsize=5.4, labelspacing=0.9,
              borderpad=0.5, handletextpad=0.8)


def main() -> None:
    df = load()

    fig = plt.figure(figsize=(W2, 6.5))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[1.35, 1.0], hspace=0.30,
        left=0.30, right=0.90, top=0.945, bottom=0.075,
    )

    ax_a = fig.add_subplot(gs[0])
    panel_dotplot(ax_a, df, "GO_BP")
    panel_label(ax_a, "a", dx=-0.42, dy=1.07)

    ax_b = fig.add_subplot(gs[1])
    panel_dotplot(ax_b, df, "KEGG")
    panel_label(ax_b, "b", dx=-0.42, dy=1.09)

    caveats = sorted(set(CAVEAT_TERMS.values()))
    fig.text(
        0.30, 0.012,
        "† " + "\n† ".join(caveats),
        fontsize=5.2, ha="left", va="bottom", color=GREY, style="italic",
        linespacing=1.5,
    )

    save(fig, "Fig04_functional_enrichment")


if __name__ == "__main__":
    main()
