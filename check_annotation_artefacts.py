"""Sanity check for GO annotation artefacts in the limbal enrichment results.

Rationale
---------
Over-representation analysis of an ECM-rich DEG list reliably returns skeletal
GO terms ("ossification", "bone development", "osteoblast differentiation"),
because GO annotates most collagens, small leucine-rich proteoglycans and
BMP/TGF-beta components to both extracellular matrix organization and skeletal
development. The term name invites the reader to infer an osteogenic program
that the data do not support.

This script tests that inference directly: if the limbal epithelium were
executing an osteogenic program, the terminal osteoblast effectors would be
induced. They are not.

Usage:  python check_ann.py            (run from the repository root)
Exit status is non-zero if any osteoblast effector is genuinely upregulated,
so the conclusion is re-checked automatically whenever the pipeline is re-run.
"""
from __future__ import annotations
import sys, pathlib, pandas as pd

REPO = pathlib.Path(__file__).resolve().parent
DEA = {"C1": "results/de/DEA_C1_Limbus_vs_Cornea.csv",
       "C2": "results/de/DEA_C2_Limbus_vs_Conjunctiva.csv",
       "C3": "results/de/DEA_C3_Limbus_vs_Cultured.csv"}
# Terminal osteoblast / osteocyte / osteoclast effectors. None of these is an
# ECM-generic gene; all are specific to bone cell differentiation.
OSTEOBLAST = ["SP7", "BGLAP", "ALPL", "IBSP", "DMP1", "SPP1", "COL10A1",
              "PHEX", "MEPE", "SOST", "ACP5", "CALCR"]
LFC_CUT, FDR_CUT = 1.0, 0.05
# "cytoskeleton" contains "skelet": exclude it explicitly.
SKELETAL = r"ossif|bone|osteo|(?<!cyto)skelet|cartilage|chondro"


def main() -> int:
    dea = {c: pd.read_csv(REPO / p).dropna(subset=["Symbol"])
                 .drop_duplicates("Symbol").set_index("Symbol")
           for c, p in DEA.items()}

    print("Osteoblast effectors, limbus versus each comparator")
    print(f"{'gene':<10}" + "".join(f"{c:>22}" for c in dea))
    induced = []
    for g in OSTEOBLAST:
        row = f"{g:<10}"
        for c, df in dea.items():
            if g not in df.index:
                row += f"{'absent':>22}"; continue
            lfc, fdr = float(df.loc[g, "logFC"]), float(df.loc[g, "adj.P.Val"])
            row += f"{lfc:>+10.2f} (FDR {fdr:>7.2g})"
            if lfc >= LFC_CUT and fdr <= FDR_CUT:
                induced.append((g, c, lfc, fdr))
        print(row)

    enr = pd.read_csv(REPO / "results/enrichment/enrichment_all.csv")
    skel = enr[enr.term.str.contains(SKELETAL, case=False, na=False)]
    print()
    for c in ("C1", "C2", "C3"):
        sub = skel[(skel.contrast == c) & (skel.direction == "Up")]
        if not len(sub):
            continue
        genes = set()
        for s in sub.genes:
            genes |= set(str(s).split("/"))
        print(f"  {c}: {len(sub):>2} skeletal terms enriched, driven by "
              f"{len(genes)} distinct genes "
              f"({len(sub) / max(len(genes), 1):.2f} terms per gene)")

    print()
    if induced:
        print("FAIL: an osteoblast effector is genuinely upregulated:")
        for g, c, lfc, fdr in induced:
            print(f"  {g} in {c}: log2FC = {lfc:+.2f}, FDR = {fdr:.2g}")
        return 1
    print("PASS: no terminal osteoblast effector is upregulated in limbus in "
          "any contrast. The skeletal GO terms reflect shared extracellular "
          "matrix and BMP annotation, not an osteogenic program.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
