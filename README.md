# MOFs-OpenKG

Compact, reproducible utilities accompanying the MOFs-OpenKG study.

This release contains only the small, reusable parts of the workflow:

- the relation schema used by the synthesis knowledge graph;
- conversion of synthesis records into deterministic KG triples;
- a dependency-light Weisfeiler–Lehman local-label refinement utility;
- MR, MRR and Hits@K ranking metrics;
- a command-line converter for preparing triples.

The original KAIST/CSD/PubChem files, CIF/SDF structures, LLM outputs,
training checkpoints and generated figures are intentionally excluded. They
are too large and may have separate redistribution terms. Put local data in a
separate ignored `data/` directory and do not commit credentials or API keys.

## Install

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -e .
```

## Quick start

```bash
mofs-openkg json-to-tsv examples/records.json triples.tsv
```

or:

```python
from mofs_openkg import record_to_triples, ranking_metrics, wl_labels
```

## Reproducibility note

The released package does not include the CompGCN training implementation,
because the supplied project directory currently contains preprocessing,
mapping and structure scripts but no standalone CompGCN trainer. Add that
trainer under `src/mofs_openkg/models/` once its source location is confirmed.
