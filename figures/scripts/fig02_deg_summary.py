"""
Figure 2 — Differential expression across the three limbus-anchored contrasts.

(a) Number of differentially expressed genes per contrast, split by direction.
(b-d) Volcano plots for C1, C2 and C3 with contrast-specific fold-change
      thresholds and the genes named in the manuscript labelled.

Inputs
  results/de/DEA_counts_summary.csv
  results/de/deg_membership.csv                     (fallback volcano source)
  results/de/DEA_C{1,2,3}_*.csv                     (preferred: full background)

Full topTables give a complete volcano cloud including non-significant genes.
When only the membership table is available the volcano is restricted to genes
significant in at least one contrast; the panel is annotated accordingly so the
restriction is never silently hidden.
"""

from __future__ import annotations

import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import (
    DIR_DE,
    FDR_CUT,
    LFC_CUT,
    PAL_CONTRAST,
    PAL_DIR,
    CONTRAST_LABEL,
    GREY,
    W2,
    apply_style,
    is_predicted,
    load_membership,
    panel_label,
    require,
    save,
)

apply_style()

CONTRASTS = ["C1", "C2", "C3"]

# Genes labelled on each volcano.
#
# IMPORTANT: these lists are contrast-specific and were verified to be
# significant IN THAT CONTRAST. The genes named in the current manuscript text
# for C2 (GEM, FMOD, SPARCL1, CRTAC1, GJB6) carry C1 fold changes and are NOT
# significant in C2; they are therefore not labelled here. See docs/ERRATA.md.
HIGHLIGHT = {
    "C1": ["POSTN", "CXCL12", "FN1", "IL6", "CDH11", "KRT3", "MAL", "KIT"],
    "C2": ["MYOC", "CYTL1", "KERA", "COL8A1", "ITGBL1", "CDH11", "PTHLH", "MMP7"],
    "C3": ["COL1A1", "MMP2", "FGF2", "JUN", "MMP9", "ITGB4", "LAMB3", "PTHLH"],
}


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #

def load_volcano(contrast: str):
    """Return (dataframe, is_full_background) for one contrast.

    Prefers the full limma topTable; falls back to the membership table.
    """
    hits = sorted(glob.glob(str(DIR_DE / f"DEA_{contrast}_*.csv")))
    hits = [h for h in hits if not h.endswith(("_UP.csv", "_DOWN.csv"))]
    if hits:
        df = pd.read_csv(hits[0])
        df = df.rename(columns={df.columns[0]: "row_id"}) if df.columns[0] == "" else df
        need = {"Symbol", "logFC", "adj.P.Val"}
        if need.issubset(df.columns):
            return df[["Symbol", "logFC", "adj.P.Val"]].copy(), True

    mm = load_membership()
    df = mm[["Symbol", f"{contrast}_logFC", f"{contrast}_adj.P.Val"]].copy()
    df.columns = ["Symbol", "logFC", "adj.P.Val"]
    return df, False


def classify(df: pd.DataFrame, contrast: str) -> pd.DataFrame:
    lfc = LFC_CUT[contrast]
    df = df.dropna(subset=["logFC", "adj.P.Val"]).copy()
    df["status"] = "NS"
    df.loc[(df["adj.P.Val"] <= FDR_CUT) & (df["logFC"] >= lfc), "status"] = "Up"
    df.loc[(df["adj.P.Val"] <= FDR_CUT) & (df["logFC"] <= -lfc), "status"] = "Down"
    # -log10(0) is undefined; clamp at the smallest representable FDR.
    floor = df.loc[df["adj.P.Val"] > 0, "adj.P.Val"].min()
    df["y"] = -np.log10(df["adj.P.Val"].clip(lower=floor))
    return df


# --------------------------------------------------------------------------- #
# Panels
# --------------------------------------------------------------------------- #

def panel_counts(ax) -> None:
    cs = pd.read_csv(require(DIR_DE / "DEA_counts_summary.csv", "DEG counts summary"))

    # Resolve the contrast identifier from whichever column carries it. The
    # current pipeline writes a "contrast" column holding C1/C2/C3 directly;
    # legacy files name rows by comparison ("Limbus_vs_AMLE"). Trying the label
    # column first fails silently on the current format, because the long label
    # contains no bare C1/C2/C3 token.
    cid = pd.Series(pd.NA, index=cs.index, dtype="object")
    for col in ("contrast", "label"):
        if col in cs.columns:
            hit = cs[col].astype(str).str.extract(r"\b(C[123])\b")[0]
            cid = cid.fillna(hit)
    if cid.isna().any() and "contrast" in cs.columns:
        order = {"Cornea": "C1", "Conjuctiva": "C2", "Conjunctiva": "C2",
                 "AMLE": "C3", "CulturedEpi": "C3", "Cultured": "C3"}
        cid = cid.fillna(cs["contrast"].astype(str).str.split("_vs_").str[-1].map(order))

    cs = cs.assign(cid=cid)
    unresolved = cs.loc[cs["cid"].isna(), :]
    if len(unresolved):
        raise SystemExit(
            "Could not map these rows of DEA_counts_summary.csv to C1/C2/C3:\n"
            + unresolved.to_string()
        )
    cs = cs.drop_duplicates(subset="cid").set_index("cid").reindex(CONTRASTS)

    x = np.arange(len(CONTRASTS))
    w = 0.36
    b1 = ax.bar(x - w / 2, cs["up"], w, color=PAL_DIR["Up"], label="Higher in limbus")
    b2 = ax.bar(x + w / 2, cs["down"], w, color=PAL_DIR["Down"], label="Lower in limbus")

    for bars in (b1, b2):
        for b in bars:
            ax.annotate(
                f"{int(b.get_height()):,}",
                (b.get_x() + b.get_width() / 2, b.get_height()),
                ha="center",
                va="bottom",
                fontsize=6.5,
                color=GREY,
                xytext=(0, 1.5),
                textcoords="offset points",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{c}\n{CONTRAST_LABEL[c]}\n|log$_2$FC| $\\geq$ {LFC_CUT[c]:g}" for c in CONTRASTS]
    )
    ax.set_ylabel("Differentially expressed genes")
    ax.set_ylim(0, cs["up"].max() * 1.18)
    ax.legend(loc="upper right")
    ax.yaxis.grid(True)
    ax.set_title("Differential expression burden per contrast", pad=5)


def panel_volcano(ax, contrast: str) -> None:
    df, full = load_volcano(contrast)
    df = classify(df, contrast)
    lfc = LFC_CUT[contrast]

    for st in ("NS", "Down", "Up"):
        sub = df[df["status"] == st]
        ax.scatter(
            sub["logFC"],
            sub["y"],
            s=2.0 if st == "NS" else 3.2,
            c=PAL_DIR[st],
            alpha=0.35 if st == "NS" else 0.75,
            linewidths=0,
            rasterized=True,
        )

    for v in (-lfc, lfc):
        ax.axvline(v, color=GREY, lw=0.5, ls=(0, (3, 2)))
    ax.axhline(-np.log10(FDR_CUT), color=GREY, lw=0.5, ls=(0, (3, 2)))

    # Headroom so labels and annotations do not collide with the frame.
    ax.set_ylim(-0.10 * df["y"].max(), df["y"].max() * 1.30)
    xr = max(abs(df["logFC"].min()), abs(df["logFC"].max())) * 1.30
    ax.set_xlim(-xr, xr)

    # Label the named genes, alternating vertical offsets to reduce overlap.
    lab = (
        df[df["Symbol"].isin(HIGHLIGHT[contrast]) & ~is_predicted(df["Symbol"])]
        .sort_values("logFC")
        .reset_index(drop=True)
    )
    for i, r in lab.iterrows():
        right = r["logFC"] > 0
        ax.annotate(
            r["Symbol"],
            (r["logFC"], r["y"]),
            fontsize=5.8,
            fontstyle="italic",
            ha="left" if right else "right",
            va="center",
            xytext=(3.5 if right else -3.5, 5.5 * (1 if i % 2 else -1)),
            textcoords="offset points",
            color="black",
            annotation_clip=False,
        )

    n_up = int((df["status"] == "Up").sum())
    n_dn = int((df["status"] == "Down").sum())
    ax.set_title(f"{contrast} · {CONTRAST_LABEL[contrast]}", color=PAL_CONTRAST[contrast], pad=3)
    ax.set_ylabel("-log$_{10}$ FDR")
    # Counts inside the frame, in the empty corner above the NS valley.
    ax.text(
        0.5, 0.985,
        f"{n_up:,} $\\uparrow$      {n_dn:,} $\\downarrow$",
        transform=ax.transAxes, va="top", ha="center",
        fontsize=6.5, color=GREY,
    )
    # Solo se conserva la advertencia sobre la restricción del universo, que es
    # una limitación de lo que el panel muestra. El carácter exploratorio de C3 se
    # declara en el pie de figura, no dentro del eje.
    notes = []
    if not full:
        notes.append("plotted cloud restricted to the DEG union")
    if notes:
        ax.text(
            0.5, -0.185,
            "\n".join(notes),
            transform=ax.transAxes, va="top", ha="center",
            fontsize=5.4, style="italic", color=GREY,
        )


# --------------------------------------------------------------------------- #
# Assemble
# --------------------------------------------------------------------------- #

def main() -> None:
    fig = plt.figure(figsize=(W2, 5.6))
    gs = fig.add_gridspec(
        2, 3, height_ratios=[1.0, 1.15], hspace=0.62, wspace=0.34,
        left=0.085, right=0.985, top=0.945, bottom=0.135,
    )

    ax_a = fig.add_subplot(gs[0, :])
    panel_counts(ax_a)
    panel_label(ax_a, "a", dx=-0.062, dy=1.12)

    for j, (c, letter) in enumerate(zip(CONTRASTS, "bcd")):
        ax = fig.add_subplot(gs[1, j])
        panel_volcano(ax, c)
        panel_label(ax, letter, dx=-0.30, dy=1.19)

    # One shared x-axis label instead of three overlapping ones.
    fig.supxlabel(
        "log$_2$ fold change (limbus / comparator)", fontsize=7.5, y=0.005
    )

    save(fig, "Fig02_DEG_summary")


if __name__ == "__main__":
    main()
