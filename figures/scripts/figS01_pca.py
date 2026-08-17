"""
Supplementary Figure 1 — Integration quality of the harmonized matrix.

(a) PCA of the batch-corrected gene-level matrix, coloured by tissue group and
    marked by platform of origin. Successful harmonization means samples separate
    by tissue, not by platform; corneal samples from both platforms should
    co-localize.
(b) Sample-sample Pearson correlation heatmap, which is what identified
    GSM724094 as a conjunctival outlier (within-group r = 0.92 versus 0.98-0.99).

Inputs
  data/processed/pca_coordinates.csv   preferred, written by R/02_integration.R
  data/processed/expr_gene.csv         used to compute PCA if the above is absent
  data/raw/sample_metadata.csv
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import (
    DIR_PROC,
    GREY,
    LIGHTGREY,
    PAL_GROUP,
    REPO,
    W2,
    apply_style,
    panel_label,
    require,
    save,
)

apply_style()

MARKER = {"GPL570": "o", "GPL6244": "^"}


def load_meta() -> pd.DataFrame:
    return pd.read_csv(require(REPO / "data" / "raw" / "sample_metadata.csv",
                               "sample metadata"))


def load_pca(meta: pd.DataFrame):
    """Return (dataframe with PC1-PC3, variance-explained triple)."""
    p = DIR_PROC / "pca_coordinates.csv"
    if p.exists():
        df = pd.read_csv(p)
        qc = REPO / "results" / "de" / "integration_qc.csv"
        vexp = (None, None, None)
        if qc.exists():
            q = pd.read_csv(qc).set_index("metric")["value"]
            vexp = tuple(
                float(q.get(f"PC{i}_var_explained", np.nan)) for i in (1, 2, 3)
            )
        return df, vexp

    expr = pd.read_csv(require(DIR_PROC / "expr_gene.csv", "gene-level matrix"),
                       index_col=0)
    x = expr.T.to_numpy(float)
    x = (x - x.mean(0)) / np.where(x.std(0) == 0, 1, x.std(0))
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u * s
    vexp = tuple((s**2 / (s**2).sum())[:3])
    df = pd.DataFrame(
        {"gsm": expr.columns, "PC1": scores[:, 0], "PC2": scores[:, 1],
         "PC3": scores[:, 2]}
    ).merge(meta[["gsm", "group", "platform"]], on="gsm", how="left")
    return df, vexp


def panel_pca(ax, df: pd.DataFrame, vexp) -> None:
    for (grp, plat), sub in df.groupby(["group", "platform"], sort=False):
        ax.scatter(
            sub["PC1"], sub["PC2"],
            s=42, color=PAL_GROUP.get(grp, LIGHTGREY),
            marker=MARKER.get(plat, "o"),
            edgecolors="white", linewidths=0.7, zorder=3,
            label=f"{grp} ({plat})",
        )

    def lab(i, v):
        return f"PC{i}" + (f"  ({100 * v:.1f}%)" if v and not np.isnan(v) else "")

    ax.set_xlabel(lab(1, vexp[0]))
    ax.set_ylabel(lab(2, vexp[1]))
    ax.axhline(0, color=LIGHTGREY, lw=0.5, zorder=1)
    ax.axvline(0, color=LIGHTGREY, lw=0.5, zorder=1)
    ax.legend(loc="best", fontsize=5.2, ncol=1)
    ax.set_title("PCA of the batch-corrected matrix", pad=5)


def panel_corr(ax, meta: pd.DataFrame) -> None:
    p = DIR_PROC / "expr_gene.csv"
    if not p.exists():
        ax.axis("off")
        ax.text(0.5, 0.5,
                "requires data/processed/expr_gene.csv",
                ha="center", va="center", fontsize=6.0, color=GREY)
        return

    expr = pd.read_csv(p, index_col=0)
    order = (
        meta.assign(g=pd.Categorical(meta["group"], categories=list(PAL_GROUP)))
        .sort_values(["g", "gsm"])["gsm"]
    )
    order = [g for g in order if g in expr.columns]
    c = expr[order].corr(method="pearson")

    im = ax.imshow(c, cmap="RdYlBu_r", vmin=c.to_numpy().min(), vmax=1.0)
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=4.8)
    ax.set_yticklabels(order, fontsize=4.8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # Group separators.
    grp = meta.set_index("gsm")["group"][order].tolist()
    bounds = [i for i in range(1, len(grp)) if grp[i] != grp[i - 1]]
    for b in bounds:
        ax.axhline(b - 0.5, color="white", lw=1.2)
        ax.axvline(b - 0.5, color="white", lw=1.2)

    ax.set_title("Sample-sample Pearson correlation", pad=5)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03, aspect=16)
    cb.set_label("$r$", fontsize=6)
    cb.ax.tick_params(labelsize=5.5, length=1.8)
    cb.outline.set_visible(False)


def main() -> None:
    meta = load_meta()
    df, vexp = load_pca(meta)

    fig = plt.figure(figsize=(W2, 3.2))
    gs = fig.add_gridspec(
        1, 2, wspace=0.34, left=0.075, right=0.94, top=0.90, bottom=0.15
    )

    ax_a = fig.add_subplot(gs[0])
    panel_pca(ax_a, df, vexp)
    panel_label(ax_a, "a", dx=-0.20, dy=1.13)

    ax_b = fig.add_subplot(gs[1])
    panel_corr(ax_b, meta)
    panel_label(ax_b, "b", dx=-0.24, dy=1.13)

    save(fig, "FigS01_integration_QC")


if __name__ == "__main__":
    main()
