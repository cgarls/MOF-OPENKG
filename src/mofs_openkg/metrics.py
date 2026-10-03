"""Ranking metrics used for relation-specific tail prediction."""

from collections.abc import Iterable


def ranking_metrics(ranks: Iterable[int], ks=(1, 3, 5, 10)) -> dict[str, float]:
    """Compute MR, MRR and Hits@K from one-based target ranks."""
    values = [int(rank) for rank in ranks]
    if not values or any(rank < 1 for rank in values):
        raise ValueError("ranks must be a non-empty iterable of positive integers")
    result = {"MR": sum(values) / len(values), "MRR": sum(1 / r for r in values) / len(values)}
    result.update({f"H@{k}": sum(r <= k for r in values) / len(values) for k in ks})
    return result
