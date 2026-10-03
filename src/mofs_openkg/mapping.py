"""Convert synthesis records into deterministic MOFs-OpenKG triples.

The mapper intentionally does not normalize synonyms: evaluation must use the
same labels as the source data unless a separately documented normalization
table is supplied.
"""

import csv
import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path

FIELD_RELATIONS = {
    "M_precursor": "hasMetal",
    "O_precursor": "hasLinker",
    "S_precursor": "hasSolvent",
    "method": "hasMethod",
    "operation": "hasOperation",
}


def stable_id(kind: str, value: object) -> str:
    """Create a reproducible identifier without exposing raw long labels."""
    text = str(value).strip()
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]
    return f"{kind}:{digest}"


def _values(value: object) -> Iterable[object]:
    if value is None or value == "":
        return ()
    if isinstance(value, (list, tuple, set)):
        return value
    return (value,)


def record_to_triples(record: Mapping, mof_key: str = "identifier") -> list[tuple[str, str, str]]:
    """Map one synthesis record to ``(head, relation, tail)`` triples."""
    mof_value = record.get(mof_key) or record.get("name") or record.get("mof")
    if mof_value is None:
        raise ValueError(f"record has no {mof_key!r}, 'name', or 'mof' field")
    head = stable_id("MOF", mof_value)
    triples = []
    for field, relation in FIELD_RELATIONS.items():
        for value in _values(record.get(field)):
            if isinstance(value, Mapping):
                value = value.get("name") or value.get("label") or value.get("id")
            if value is not None and str(value).strip():
                triples.append((head, relation, stable_id(relation, value)))
    return triples


def records_to_tsv(records: Iterable[Mapping], output: str | Path) -> int:
    """Write triples as a compact TSV file and return the number of rows."""
    count = 0
    with Path(output).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("head", "relation", "tail"))
        for record in records:
            for triple in record_to_triples(record):
                writer.writerow(triple)
                count += 1
    return count


def json_to_tsv(input_path: str | Path, output_path: str | Path) -> int:
    """Convert a JSON list of synthesis records to triples."""
    with Path(input_path).open(encoding="utf-8") as handle:
        records = json.load(handle)
    if not isinstance(records, list):
        raise ValueError("input JSON must contain a list of records")
    return records_to_tsv(records, output_path)
