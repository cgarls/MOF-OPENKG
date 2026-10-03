"""The compact relation schema used in the MOFs-OpenKG release."""

RELATIONS = {
    "hasMetal": ("M_precursor", "Precursor"),
    "hasLinker": ("O_precursor", "Precursor"),
    "hasSolvent": ("S_precursor", "Precursor"),
    "hasMethod": ("method", "Method"),
    "hasOperation": ("operation", "Operation"),
    "hasPaper": ("publication identifier", "Paper"),
    "hasAuthor": ("paper authorship", "Author"),
    "hasJournal": ("paper journal", "Journal"),
    "hasColor": ("crystal descriptor", "Color"),
    "hasHabit": ("crystal descriptor", "Habit"),
    "hasBond": ("atom–bond structure", "Bond"),
    "hasKernel": ("local environment", "Kernel"),
    "subKernel": ("kernel hierarchy", "Kernel"),
}


def relation_schema():
    """Return rows suitable for exporting the relation table."""
    return [
        {"source_field_or_object": source, "graph_relation": relation,
         "target_type": target}
        for relation, (source, target) in RELATIONS.items()
    ]
