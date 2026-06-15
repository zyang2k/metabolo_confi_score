"""Set Transformer MVP for evidence fusion in LC-MS/MS confidence scoring.

Modules:
    modules.py  — MAB / ISAB / PMA building blocks (masked).
    dataset.py  — (B, K, d) tensor builder over feature_table_v2.csv.
    model.py    — EvidenceSetTransformer with equivariant + set-pooled heads.
"""
