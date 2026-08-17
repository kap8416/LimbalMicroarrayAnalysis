"""
Figure 6 — A compact ECM-adhesion circuit distinguishes limbus from conjunctiva
(C2), whereas ex vivo culture reorganizes the interactome broadly (C3).

(a) C2 node-link diagram. The full topology of this 5-node network, with LAMB3
    as the central node.
(b) C3 composite hub ranking, separating limbus-upregulated ECM hubs from
    culture-associated integrin hubs.
(c) C3 module landscape.
(d) C3 node-link diagram, when the edge list is available.

Inputs
  results/ppi/centrality_C2.csv, results/ppi/edges_C2.csv
  results/ppi/centrality_C3.csv, results/ppi/edges_C3.csv  (optional for panel d)

Run reconstruct_c2_edges.py first if results/ppi/edges_C2.csv is absent.
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from fig_style import (
    GREY,
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

# Modules named in Results 3.3.3.
C3_NAMED_MODULES = {
    "COL1A1": "ECM remodelling",
    "FGF2": "Immune and paracrine\nsignalling",
    "JUN": "Stress response and\ntranscriptional regulation",
    "OCLN": "Barrier and differentiation\n(up in cultured cells)",
}


def main() -> None:
    c2 = load_centrality("C2")
    c3 = load_centrality("C3")

    fig = plt.figure(figsize=(W2, 6.4))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.0, 1.20], hspace=0.40, wspace=0.34,
        left=0.10, right=0.975, top=0.905, bottom=0.085,
    )

    # --- a: C2, complete topology ------------------------------------------ #
    ax_a = fig.add_subplot(gs[0, 0])
    if have_edges("C2"):
        draw_network(
            ax_a,
            "C2",
            c2,
            label_genes=list(c2["gene"]),
            title="C2 · limbus vs conjunctiva\nLAMB3-centred ECM-adhesion circuit",
            reconstructed=True,
            layout="manual",
        )
    else:
        missing_edges_note(ax_a, "C2")
    panel_label(ax_a, "a", dx=-0.10, dy=1.10)

    # --- b: C3 hubs --------------------------------------------------------- #
    ax_b = fig.add_subplot(gs[0, 1])
    panel_hub_ranking(ax_b, c3, title="C3 · top 15 hubs")
    panel_label(ax_b, "b", dx=-0.30, dy=1.10)

    # --- c: C3 modules ------------------------------------------------------ #
    ax_c = fig.add_subplot(gs[1, 0])
    panel_module_landscape(ax_c, c3, C3_NAMED_MODULES)
    panel_label(ax_c, "c", dx=-0.26, dy=1.10)

    # --- d: C3 network ------------------------------------------------------ #
    ax_d = fig.add_subplot(gs[1, 1])
    if have_edges("C3"):
        draw_network(
            ax_d,
            "C3",
            c3,
            label_genes=[],
            focus_hubs=True,
            title="C3 interactome: neighbourhood of the leading hubs",
        )
    else:
        missing_edges_note(ax_d, "C3")
    panel_label(ax_d, "d", dx=-0.08, dy=1.10)

    fig.text(
        0.5, 0.005,
        "C3 is exploratory: the contrast is fully confounded with platform of origin.",
        ha="center", va="bottom", fontsize=5.8, style="italic", color=GREY,
    )

    save(fig, "Fig06_PPI_C2_C3")


if __name__ == "__main__":
    main()
