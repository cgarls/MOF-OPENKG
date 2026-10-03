from mofs_openkg.mapping import record_to_triples
from mofs_openkg.metrics import ranking_metrics
from mofs_openkg.wl import wl_labels


def test_mapping_and_metrics():
    triples = record_to_triples({"identifier": "x", "M_precursor": ["CuCl2"]})
    assert triples[0][1] == "hasMetal"
    assert ranking_metrics([1, 2, 10])["H@1"] == 1 / 3


def test_wl_is_deterministic():
    first = wl_labels([1, 2, 3], [(1, 2), (2, 3)], iterations=2)
    second = wl_labels([1, 2, 3], [(1, 2), (2, 3)], iterations=2)
    assert first == second
