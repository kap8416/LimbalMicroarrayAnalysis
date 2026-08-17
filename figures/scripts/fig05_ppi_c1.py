"""
Figure 5 — The limbus-versus-cornea (C1) interactome is anchored on FN1.

(a) Composite hub ranking of the 15 top-scoring nodes, coloured by direction and
    annotated with degree.
(b) Module landscape: module size versus the proportion of limbus-upregulated
    nodes, with the named functional modules highlighted.
(c) Node-link diagram of the network. Drawn only when an edge list is available
    (results/ppi/edges_C1.csv); otherwise the panel reports what is missing
    rather than substituting a decorative layout.

Inputs
  results/ppi/centrality_C1.csv     (required)
  results/ppi/edges_C1.csv          (optional, enables panel c)
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import (
    DIR_PPI,
    GREY,
    LIGHTGREY,
    PAL_DIR,
    W2,
    apply_style,
    load_centrality,
    panel_label,
    save,
)
from ppi_common import (
    draw_network,
    have_edges,
    missing_edges_note,
    panel_hub_ranking,
    panel_module_landscape,
)

apply_style()

CONTRAST = "C1"

# Modules named in Results 3.3.1, keyed by a representative member so the
# annotation survives changes in Louvain module numbering between runs.
NAMED_MODULES = {
    "FN1": "ECM organization\nand remodelling",
    "IL6": "Immune surveillance\nand chemotaxis",
    "BMP4": "BMP / Wnt\ndevelopmental regulators",
    "NT5E": "cAMP-calcium\nsignalling",
}


def main() -> None:
    cen = load_centrality(CONTRAST)

    fig = plt.figure(figsize=(W2, 5.9))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.0],
        hspace=0.42, wspace=0.34,
        left=0.10, right=0.975, top=0.93, bottom=0.075,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    panel_hub_ranking(ax_a, cen, title="Top 15 hubs by composite centrality")
    panel_label(ax_a, "a", dx=-0.30, dy=1.14)

    ax_b = fig.add_subplot(gs[0, 1])
    panel_module_landscape(ax_b, cen, NAMED_MODULES)
    panel_label(ax_b, "b", dx=-0.24, dy=1.14)

    ax_c = fig.add_subplot(gs[1, :])
    if have_edges(CONTRAST):
        draw_network(
            ax_c,
            CONTRAST,
            cen,
            label_genes=[],
            focus_hubs=True,
            title="C1 interactome: neighbourhood of the leading hubs",
        )
    else:
        missing_edges_note(ax_c, CONTRAST)
    panel_label(ax_c, "c", dx=-0.055, dy=1.06)

    save(fig, "Fig05_PPI_C1")


if __name__ == "__main__":
    main()
