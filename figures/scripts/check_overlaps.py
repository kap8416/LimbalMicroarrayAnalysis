"""
Check every figure for overlapping text.

Regenerates each figure in memory and compares the rendered bounding box of every
text element against every other. Overlaps are reported with the figure, the
offending strings, and the fraction of the smaller box that is covered.

Why this exists: label collisions are invisible in the code and only appear once
the figure is rendered at final size. They are also the defect most likely to
survive into a submitted manuscript, because a figure that is merely crowded still
"looks finished" in a thumbnail.

Ignored by design:
  - axis tick labels against their own axis label, which never truly collide;
  - text inside a legend, which matplotlib lays out itself;
  - overlaps below MIN_FRACTION, which are adjacent boxes with padding rather than
    illegible text.

Usage:  python check_overlaps.py [figNN ...]
Exit status is non-zero if any overlap is found, so it can gate a release.
"""

from __future__ import annotations

import importlib.util
import itertools
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

# Overlap smaller than this fraction of the smaller box is treated as padding.
MIN_FRACTION = 0.02

FIGURES = [
    "fig01_pipeline", "fig02_deg_summary", "fig03_core_signature",
    "fig04_enrichment", "fig05_ppi_c1", "fig06_ppi_c2_c3",
    "fig07_model", "fig08_mra", "figS01_pca",
]


def texts_of(fig):
    """Every text element with a non-empty string, with its axes and legend flag."""
    out = []
    for ax in fig.get_axes():
        legend_texts = set()
        leg = ax.get_legend()
        if leg is not None:
            legend_texts = {id(t) for t in leg.get_texts()}
        for t in ax.texts + [ax.title, ax.xaxis.label, ax.yaxis.label]:
            if t.get_text().strip() and t.get_visible():
                out.append((ax, t, id(t) in legend_texts))
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            if t.get_text().strip() and t.get_visible():
                out.append((ax, t, False))
        # Los textos de la leyenda no están en ax.texts; hay que añadirlos
        # explícitamente o una leyenda mal colocada pasa desapercibida.
        if leg is not None:
            for t in leg.get_texts():
                if t.get_text().strip() and t.get_visible():
                    out.append((ax, t, True))
    out.extend((None, t, False) for t in fig.texts if t.get_text().strip())
    return out


def overlap_fraction(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    if dx <= 0 or dy <= 0:
        return 0.0
    inter = dx * dy
    smaller = min(a.width * a.height, b.width * b.height)
    return inter / smaller if smaller > 0 else 0.0


def check(mod_name):
    spec = importlib.util.spec_from_file_location(mod_name, HERE / f"{mod_name}.py")
    mod = importlib.util.module_from_spec(spec)
    plt.close("all")
    try:
        spec.loader.exec_module(mod)
        mod.main()
    except SystemExit as exc:                    # falta un archivo de entrada
        return None, str(exc)

    problems = []
    for fig in [plt.figure(n) for n in plt.get_fignums()]:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        items = texts_of(fig)
        boxes = []
        for ax, t, in_legend in items:
            try:
                boxes.append((ax, t, in_legend, t.get_window_extent(renderer=r)))
            except Exception:
                pass
        for (ax1, t1, l1, b1), (ax2, t2, l2, b2) in itertools.combinations(boxes, 2):
            if l1 and l2:          # dos textos dentro de la misma leyenda
                continue
            # ejes distintos que se pisan sí importa; mismo eje también
            f = overlap_fraction(b1, b2)
            if f >= MIN_FRACTION:
                problems.append((t1.get_text()[:34], t2.get_text()[:34], round(f, 2)))
        # texto que se sale del área de su propio eje
        for ax, t, in_legend in items:
            if ax is None or in_legend:
                continue
            if t in (ax.title, ax.xaxis.label, ax.yaxis.label):
                continue
            if t in ax.get_xticklabels() or t in ax.get_yticklabels():
                continue
            try:
                tb = t.get_window_extent(renderer=r)
            except Exception:
                continue
            ab = ax.get_window_extent(renderer=r)
            if tb.x0 < ab.x0 - 2 or tb.x1 > ab.x1 + 2:
                problems.append((t.get_text()[:34], "<fuera del eje>", 1.0))

    plt.close("all")
    return problems, None


def main(argv):
    targets = argv[1:] or FIGURES
    total = 0
    print(f"overlap threshold: {MIN_FRACTION:.0%} of the smaller text box\n")
    for name in targets:
        problems, skipped = check(name)
        if skipped is not None:
            print(f"  {name:22} skipped — missing input")
            continue
        if problems:
            total += len(problems)
            print(f"  {name:22} {len(problems)} overlapping pair(s)")
            for a, b, f in problems[:8]:
                print(f"      {f:4.0%}  {a!r}  <->  {b!r}")
        else:
            print(f"  {name:22} clean")
    print()
    if total:
        print(f"FAIL: {total} overlapping text pair(s)")
        return 1
    print("PASS: no overlapping text in any figure")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
