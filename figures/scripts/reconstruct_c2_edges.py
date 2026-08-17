"""
Reconstruct the C2 (limbus vs conjunctiva) PPI edge list from its centrality table.

The C2 network has only 5 nodes, and its degree, betweenness and closeness values
jointly determine a unique topology. This script enumerates every simple graph
consistent with the observed degree sequence and keeps the one that reproduces
all three centrality vectors exactly, so the resulting edge list is derived, not
assumed.

Why this exists: the manuscript reports "5 nodes and 10 edges" (revised down from
"20 edges"), but the degree sequence in results/ppi/centrality_C2.csv sums to 8,
which is 4 edges. The reconstruction below confirms 4 edges independently, via
betweenness and closeness. See docs/ERRATA.md.

Output
  results/ppi/edges_C2.csv
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from fig_style import DIR_PPI, load_centrality

TOL = 1e-6


def all_pairs_shortest_paths(n: int, adj: np.ndarray) -> np.ndarray:
    """Floyd-Warshall on an unweighted adjacency matrix."""
    d = np.where(adj > 0, 1.0, np.inf)
    np.fill_diagonal(d, 0.0)
    for k in range(n):
        d = np.minimum(d, d[:, [k]] + d[[k], :])
    return d


def centralities(adj: np.ndarray):
    """Normalized degree, betweenness and closeness for a small simple graph."""
    n = adj.shape[0]
    deg = adj.sum(axis=1)
    d = all_pairs_shortest_paths(n, adj)

    if np.isinf(d).any():
        return None  # disconnected: closeness is undefined on this scale

    # igraph normalizes degree by (n - 1) and closeness as (n - 1) / sum(d).
    deg_norm = deg / (n - 1)
    closeness = (n - 1) / d.sum(axis=1)

    # Betweenness by explicit enumeration of shortest paths (n is tiny).
    btw = np.zeros(n)
    for s, t in itertools.combinations(range(n), 2):
        paths = _shortest_paths(adj, s, t, int(d[s, t]))
        if not paths:
            continue
        for v in range(n):
            if v in (s, t):
                continue
            hit = sum(1 for p in paths if v in p[1:-1])
            btw[v] += hit / len(paths)
    # igraph normalizes undirected betweenness by (n-1)(n-2)/2.
    btw_norm = btw / ((n - 1) * (n - 2) / 2)

    return deg, deg_norm, btw_norm, closeness


def _shortest_paths(adj, s, t, length):
    """All shortest s-t paths of the given length."""
    out, stack = [], [(s, [s])]
    while stack:
        node, path = stack.pop()
        if len(path) - 1 > length:
            continue
        if node == t:
            if len(path) - 1 == length:
                out.append(path)
            continue
        for nb in np.flatnonzero(adj[node]):
            if nb not in path:
                stack.append((int(nb), path + [int(nb)]))
    return out


def main() -> None:
    cen = load_centrality("C2")
    genes = list(cen["gene"])
    n = len(genes)
    assert n <= 8, "This exhaustive reconstruction is intended for very small networks only."

    obs_deg = cen["degree"].to_numpy(float)
    obs_btw = cen["betweenness"].to_numpy(float)
    obs_clo = cen["closeness"].to_numpy(float)
    n_edges = int(obs_deg.sum() // 2)

    print(f"C2 nodes: {n}  ({', '.join(genes)})")
    print(f"degree sequence: {obs_deg.astype(int).tolist()}  ->  sum {int(obs_deg.sum())}"
          f"  ->  {n_edges} edges")

    pairs = list(itertools.combinations(range(n), 2))
    solutions = []

    for combo in itertools.combinations(pairs, n_edges):
        adj = np.zeros((n, n))
        for i, j in combo:
            adj[i, j] = adj[j, i] = 1
        if not np.allclose(adj.sum(axis=1), obs_deg):
            continue
        res = centralities(adj)
        if res is None:
            continue
        _, _, btw, clo = res
        if np.allclose(btw, obs_btw, atol=1e-4) and np.allclose(clo, obs_clo, atol=1e-4):
            solutions.append(combo)

    print(f"topologies consistent with degree + betweenness + closeness: {len(solutions)}")
    assert solutions, (
        "No topology reproduces the observed centralities. Export the real edge "
        "list from the igraph object in R/05_ppi_networks.R instead."
    )

    # Distinct solutions differ only by relabelling of equivalent nodes; the
    # edge set is identical up to automorphism, so the first is representative.
    combo = solutions[0]
    mods = cen.set_index("gene")["module"]
    edges = pd.DataFrame(
        [
            {
                "from": genes[i],
                "to": genes[j],
                "same_module": bool(mods[genes[i]] == mods[genes[j]]),
            }
            for i, j in combo
        ]
    )

    out = DIR_PPI / "edges_C2.csv"
    edges.to_csv(out, index=False)
    print(f"\nreconstructed edge list ({len(edges)} edges) -> {out}")
    print(edges.to_string(index=False))


if __name__ == "__main__":
    main()
