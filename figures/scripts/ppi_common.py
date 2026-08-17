"""
Shared panel builders for the PPI figures (Fig. 5 and Fig. 6).

Two design decisions are enforced here:

1. Node and edge counts are always read from the same object. Quoting a node
   count from the retained connected component alongside an edge count from the
   full STRING response is how network sizes get misreported; every count shown
   by these panels is derived from the centrality table (nodes) and the edge
   list (edges) of that same component.

2. When an edge list is absent, the node-link panel is NOT faked with a random
   or force-directed layout over invented edges. It states what is missing.
   For small networks whose degree sequence uniquely determines the topology
   (C2), the edge list is reconstructed and the reconstruction is labelled.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fig_style import DIR_PPI, GREY, LIGHTGREY, PAL_DIR


# --------------------------------------------------------------------------- #
# Edge-list availability
# --------------------------------------------------------------------------- #

def edge_path(contrast: str) -> Path:
    return DIR_PPI / f"edges_{contrast}.csv"


def have_edges(contrast: str) -> bool:
    return edge_path(contrast).exists()


def load_edges(contrast: str) -> pd.DataFrame:
    df = pd.read_csv(edge_path(contrast))
    cols = {c.lower(): c for c in df.columns}
    a = cols.get("from") or cols.get("node1") or cols.get("preferredname_a") or df.columns[0]
    b = cols.get("to") or cols.get("node2") or cols.get("preferredname_b") or df.columns[1]
    out = df[[a, b]].copy()
    out.columns = ["from", "to"]
    return out


def missing_edges_note(ax, contrast: str) -> None:
    """Explicit placeholder when the edge list has not been exported."""
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        f"Node-link diagram for {contrast}\nrequires the edge list.\n\n"
        f"Export it from R/05_ppi_networks.R:\n\n"
        f"write.csv(\n"
        f"  igraph::as_data_frame(g, \"edges\"),\n"
        f"  \"results/ppi/edges_{contrast}.csv\",\n"
        f"  row.names = FALSE)\n\n"
        f"The other panels are complete and\nuse the centrality table only.",
        ha="center",
        va="center",
        fontsize=6.0,
        color=GREY,
        linespacing=1.6,
        family="monospace",
    )


# --------------------------------------------------------------------------- #
# Panel: hub ranking
# --------------------------------------------------------------------------- #

def panel_hub_ranking(ax, cen: pd.DataFrame, title: str, n: int = 15) -> None:
    """Horizontal bars of composite hub score for the top-n nodes."""
    d = cen.nlargest(n, "hub_score").sort_values("hub_score")
    y = np.arange(len(d))
    colors = [PAL_DIR.get(s, LIGHTGREY) for s in d["direction"]]

    ax.barh(y, d["hub_score"], 0.66, color=colors)

    for yi, (_, r) in zip(y, d.iterrows()):
        ax.annotate(
            f"k={int(r['degree'])}   log$_2$FC {r['logFC']:+.2f}",
            (r["hub_score"], yi),
            va="center", ha="left", fontsize=5.4, color=GREY,
            xytext=(2.5, 0), textcoords="offset points",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(d["gene"], fontstyle="italic", fontsize=6.4)
    ax.set_xlabel("Composite hub score")
    ax.set_xlim(0, 1.42)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.grid(True)
    ax.set_title(title, pad=7, loc="left")

    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=PAL_DIR["Up"]),
        plt.Rectangle((0, 0), 1, 1, fc=PAL_DIR["Down"]),
    ]
    # La leyenda va sobre el eje: dentro chocaba con las barras inferiores y con
    # sus anotaciones de grado y log2FC.
    # Título a la izquierda y leyenda a la derecha, ambos sobre el eje: dentro del
    # área de barras la leyenda tapaba las filas inferiores y sus anotaciones.
    ax.legend(handles, ["Up in limbus", "Down in limbus"],
              loc="lower right", bbox_to_anchor=(1.0, 1.0),
              ncol=2, fontsize=5.6, borderaxespad=0.0, handlelength=0.9,
              columnspacing=0.9, handletextpad=0.4)


# --------------------------------------------------------------------------- #
# Panel: module landscape
# --------------------------------------------------------------------------- #

def panel_module_landscape(ax, cen: pd.DataFrame, named: dict[str, str]) -> None:
    """Module size versus the proportion of limbus-upregulated nodes.

    Reading modules on these two axes separates modules that carry the limbal
    program (large and near-uniformly upregulated) from modules that mark the
    comparator state (predominantly downregulated), which a node-link diagram
    alone does not make quantitative.
    """
    mod = (
        cen.groupby("module")
        .agg(
            n=("gene", "size"),
            pct_up=("direction", lambda s: 100 * (s == "Up").mean()),
            top=("hub_score", "max"),
        )
        .reset_index()
    )
    # Representative gene per module, used to attach the functional labels.
    rep = cen.sort_values("hub_score", ascending=False).groupby("module")["gene"].first()
    mod["rep"] = mod["module"].map(rep)

    sizes = 14 + 46 * (mod["n"] / mod["n"].max())
    sc = ax.scatter(
        mod["n"], mod["pct_up"], s=sizes,
        c=mod["pct_up"], cmap="RdBu_r", vmin=0, vmax=100,
        edgecolors="white", linewidths=0.5, zorder=3,
    )
    ax.axhline(50, color=GREY, lw=0.5, ls=(0, (3, 2)), zorder=1)

    # Label the modules named in the Results by their representative member.
    # Offsets are staggered because the largest, most upregulated modules cluster
    # in the same corner of this space and their labels would otherwise collide.
    xmax = mod["n"].max() * 1.42
    stagger = [(9, 14), (9, -16), (9, 20), (9, -22), (9, 0)]
    for oi, (key, lab) in enumerate(named.items()):
        hit = cen.loc[cen["gene"] == key, "module"]
        if hit.empty:
            continue
        m = mod[mod["module"] == hit.iloc[0]]
        if m.empty:
            continue
        r = m.iloc[0]
        dx, dy = stagger[oi % len(stagger)]
        # Si no hay espacio a la derecha, etiquetar hacia la izquierda.
        if r["n"] > 0.62 * xmax:
            dx, ha = -abs(dx), "right"
        else:
            ha = "left"
        ax.annotate(
            lab, (r["n"], r["pct_up"]),
            fontsize=5.3, ha=ha, va="center", color="black",
            xytext=(dx, dy), textcoords="offset points", annotation_clip=False,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY, shrinkA=0, shrinkB=3),
        )

    ax.set_xlabel("Nodes per Louvain module")
    ax.set_ylabel("Nodes upregulated in limbus (%)")
    ax.set_ylim(-6, 106)
    ax.set_xlim(-mod["n"].max() * 0.04, mod["n"].max() * 1.42)
    ax.yaxis.grid(True)
    ax.set_title(f"Module landscape ({len(mod)} modules)", pad=5)


# --------------------------------------------------------------------------- #
# Panel: node-link diagram
# --------------------------------------------------------------------------- #

def _spring_layout(nodes, edges, seed=42, iters=600):
    """Minimal Fruchterman-Reingold layout, so networkx is not a dependency."""
    rng = np.random.default_rng(seed)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    pos = rng.uniform(-1, 1, size=(n, 2))
    ei = np.array([[idx[a], idx[b]] for a, b in edges if a in idx and b in idx])
    k = 1.6 / np.sqrt(max(n, 1))   # mayor separación entre nodos
    t = 0.15

    for _ in range(iters):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(delta, axis=-1)
        np.fill_diagonal(dist, np.inf)
        rep = (k**2 / dist**2)[:, :, None] * delta
        disp = rep.sum(axis=1)

        if len(ei):
            d = pos[ei[:, 0]] - pos[ei[:, 1]]
            dd = np.linalg.norm(d, axis=1, keepdims=True)
            dd[dd == 0] = 1e-9
            att = (dd / k) * d
            np.add.at(disp, ei[:, 0], -att)
            np.add.at(disp, ei[:, 1], att)

        norm = np.linalg.norm(disp, axis=1, keepdims=True)
        norm[norm == 0] = 1e-9
        pos += disp / norm * np.minimum(norm, t)
        t *= 0.994

    pos -= pos.mean(axis=0)
    span = np.abs(pos).max()
    return pos / (span if span else 1.0)




def hub_neighbourhood(edges, cen, n_hubs=8, max_nodes=140):
    """Subgraph induced by the top hubs and their immediate neighbours.

    A 1,285-node network drawn at figure scale is a hairball: it shows density but
    no structure, and no label can be placed legibly in its core. Showing the
    neighbourhood of the leading hubs instead displays the part of the network the
    text actually discusses, at a scale where nodes and labels are readable. The
    full node and edge counts remain in the panel annotation, so nothing is hidden.

    Neighbours are ranked by degree so the subgraph keeps the well-connected
    context rather than an arbitrary sample of pendant nodes.
    """
    hubs = list(cen.nlargest(n_hubs, "hub_score")["gene"])
    deg = cen.set_index("gene")["degree"]
    hubset = set(hubs)

    nbrs = set()
    for a, b in edges[["from", "to"]].itertuples(index=False, name=None):
        if a in hubset: nbrs.add(b)
        if b in hubset: nbrs.add(a)
    nbrs -= hubset

    room = max(0, max_nodes - len(hubs))
    keep = set(hubs) | set(sorted(nbrs, key=lambda g: -deg.get(g, 0))[:room])
    sub = edges[edges["from"].isin(keep) & edges["to"].isin(keep)]
    return sub, hubs


def _declutter_labels(ax, anchors, texts, fontsize, seed=0, iters=1200):
    """Place labels near their anchors without overlapping each other.

    Labels start at their node and are then pushed apart by a short force
    relaxation: each label repels every other label whose box would overlap, and
    is pulled back toward its own anchor. A thin leader line connects the final
    label position to the node, so displacement never creates ambiguity about
    which node a label belongs to.

    Crowded node-link panels otherwise stack several gene names on top of one
    another in the dense core of the graph, which is unreadable and is the single
    most common defect in published network figures.
    """
    if not len(anchors):
        return
    fig = ax.figure
    fig.canvas.draw()                       # necesario para medir texto
    inv = ax.transData.inverted()

    # tamaño aproximado de cada caja de texto, en coordenadas de datos
    r = fig.canvas.get_renderer()
    sizes = []
    for t in texts:
        tp = ax.text(0, 0, t, fontsize=fontsize, fontstyle="italic")
        bb = tp.get_window_extent(renderer=r)
        p0 = inv.transform((0, 0)); p1 = inv.transform((bb.width, bb.height))
        sizes.append((abs(p1[0] - p0[0]) * 1.45, abs(p1[1] - p0[1]) * 2.8))
        tp.remove()
    sizes = np.array(sizes)

    pos = np.array(anchors, dtype=float)
    rng = np.random.default_rng(seed)
    pos += rng.normal(0, 1e-3, pos.shape)   # romper empates exactos
    anch = np.array(anchors, dtype=float)

    for _ in range(iters):
        disp = np.zeros_like(pos)
        d = pos[:, None, :] - pos[None, :, :]
        need = (sizes[:, None, :] + sizes[None, :, :]) / 2
        overlap = (np.abs(d) < need).all(axis=-1)
        np.fill_diagonal(overlap, False)
        for i, j in zip(*np.where(overlap)):
            v = d[i, j]
            n = np.linalg.norm(v)
            v = v / n if n > 1e-9 else rng.normal(0, 1, 2)
            disp[i] += v * need[i, j].mean() * 0.85
        disp += (anch - pos) * 0.05          # atracción al nodo, débil
        pos += disp * 0.5

    for (x, y), (ax_, ay), t in zip(pos, anch, texts):
        if np.hypot(x - ax_, y - ay) > 1e-3:
            ax.plot([ax_, x], [ay, y], color=GREY, lw=0.35, alpha=0.65, zorder=3.5)
        ax.text(x, y, t, fontsize=fontsize, fontstyle="italic",
                ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.88))

def draw_network(
    ax,
    contrast: str,
    cen: pd.DataFrame,
    label_genes: list[str],
    title: str,
    reconstructed: bool = False,
    layout: str | None = None,
    focus_hubs: bool = False,
) -> None:
    all_edges = load_edges(contrast)
    keep_all = set(cen["gene"])
    all_edges = all_edges[all_edges["from"].isin(keep_all) & all_edges["to"].isin(keep_all)]
    n_full, e_full = len(keep_all), len(all_edges)

    if focus_hubs:
        edges, label_genes = hub_neighbourhood(all_edges, cen)
    else:
        edges = all_edges

    nodes = sorted({x for e in edges[["from", "to"]].itertuples(index=False, name=None) for x in e})

    if layout == "manual" and len(nodes) <= 8:
        # Small networks read better on a fixed radial layout than on a
        # stochastic one.
        centre = cen.nlargest(1, "degree")["gene"].iloc[0]
        others = [g for g in nodes if g != centre]
        ang = np.linspace(0, 2 * np.pi, len(others), endpoint=False) + np.pi / 6
        pos = {centre: np.array([0.0, 0.0])}
        pos.update({g: np.array([np.cos(a), np.sin(a)]) for g, a in zip(others, ang)})
        xy = np.array([pos[g] for g in nodes])
    else:
        xy = _spring_layout(nodes, list(edges.itertuples(index=False, name=None)))

    pos = {g: xy[i] for i, g in enumerate(nodes)}

    for a, b in edges.itertuples(index=False, name=None):
        pa, pb = pos[a], pos[b]
        ax.plot(
            [pa[0], pb[0]], [pa[1], pb[1]],
            color="#B9C4C2", lw=0.45 if len(nodes) > 60 else 1.1,
            alpha=0.5 if len(nodes) > 60 else 0.9, zorder=1,
            solid_capstyle="round",
        )

    deg = cen.set_index("gene")["degree"]
    smax = deg.max()
    sizes = (14 + 190 * (deg[nodes].to_numpy() / smax) if len(nodes) > 60
             else 110 + 320 * (deg[nodes].to_numpy() / smax))
    colors = [PAL_DIR.get(s, LIGHTGREY) for s in cen.set_index("gene")["direction"][nodes]]

    ax.scatter(
        xy[:, 0], xy[:, 1], s=sizes, c=colors,
        edgecolors="white", linewidths=0.5 if len(nodes) > 60 else 0.8,
        zorder=3,
    )

    small = len(nodes) <= 60
    lab = [g for g in dict.fromkeys(label_genes) if g in pos]

    if small:
        # Redes pequeñas: desplazamiento radial, que es predecible y legible.
        for g in lab:
            q = pos[g]
            rr = np.linalg.norm(q)
            u = q / rr if rr > 1e-9 else np.array([0.0, 1.0])
            ha = "left" if u[0] > 0.25 else ("right" if u[0] < -0.25 else "center")
            va = "bottom" if u[1] > 0.25 else ("top" if u[1] < -0.25 else "center")
            ax.annotate(g, q, fontsize=7.0, fontstyle="italic", ha=ha, va=va,
                        xytext=tuple(u * 20), textcoords="offset points", zorder=5,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.9))
    else:
        # Redes densas: relajación de etiquetas con línea guía.
        _declutter_labels(ax, [pos[g] for g in lab], lab, fontsize=6.0)

    ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
    # Margin so radial labels stay inside the panel.
    m = 0.62 if small else 0.30
    ax.set_xlim(xy[:, 0].min() - m, xy[:, 0].max() + m)
    ax.set_ylim(xy[:, 1].min() - m, xy[:, 1].max() + m)
    ax.set_title(title, pad=4)

    if focus_hubs:
        note = (f"neighbourhood of the {len(label_genes)} leading hubs: "
                f"{len(nodes)} nodes · {len(edges)} edges\n"
                f"full network {n_full} nodes · {e_full} edges · "
                f"{cen['module'].nunique()} modules")
    else:
        note = f"{len(nodes)} nodes · {len(edges)} edges · {cen['module'].nunique()} modules"
    if reconstructed:
        note += "\ntopology reconstructed from the degree sequence"
    ax.text(
        0.5, -0.02, note, transform=ax.transAxes, ha="center", va="top",
        fontsize=5.8, color=GREY, style="italic" if reconstructed else "normal",
    )
