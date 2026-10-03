"""Core, dependency-light utilities for the MOFs-OpenKG project."""

from .schema import RELATIONS, relation_schema
from .mapping import record_to_triples, records_to_tsv
from .metrics import ranking_metrics
from .wl import wl_labels

__all__ = [
    "RELATIONS",
    "relation_schema",
    "record_to_triples",
    "records_to_tsv",
    "ranking_metrics",
    "wl_labels",
]
