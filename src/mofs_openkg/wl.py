"""Small pure-Python Weisfeiler–Lehman label refinement utility."""

from collections.abc import Iterable, Mapping
import hashlib


def _digest(value: str) -> str:
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:16]


def wl_labels(nodes: Iterable[object], edges: Iterable[tuple[object, object]],
              iterations: int = 3, node_labels: Mapping | None = None) -> list[dict[object, str]]:
    """Return node labels for iteration 0 through ``iterations``.

    Edges are treated as undirected, matching the local-neighborhood use in
    the MOFs-OpenKG WL fingerprint.  The function has no third-party
    dependency so the release package remains easy to install.
    """
    nodes = list(nodes)
    neighbors = {node: set() for node in nodes}
    for left, right in edges:
        neighbors.setdefault(left, set()).add(right)
        neighbors.setdefault(right, set()).add(left)
    labels = {node: str((node_labels or {}).get(node, node)) for node in neighbors}
    history = [dict(labels)]
    for _ in range(iterations):
        labels = {
            node: _digest(labels[node] + "|" + "|".join(sorted(labels[n] for n in neighbors[node])))
            for node in neighbors
        }
        history.append(dict(labels))
    return history
