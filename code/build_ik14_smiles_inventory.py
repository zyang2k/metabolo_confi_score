"""Stage 1: build IK14 -> canonical SMILES inventory for labeled-bin scope.

Strategy:
  feature_table_v2.csv has hit_ik14 but no smiles.
  orbitrap_hits_v2.csv has smiles but no hit_ik14.
  Both share (wiki_id, library_wiki_id) -- one row per library match.

  Inner-join on that key, restrict to labeled (non-blank) bins, group by
  hit_ik14, pick most-common canonical SMILES per IK14.

Output:
  data/ik14_to_smiles_labeled.csv  [hit_ik14, smiles, n_occurrences]
  data/ik14_missing_smiles.csv     IK14s in scope but with no SMILES found
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

FT_PATH = Path("data/feature_table_v2.csv")
OH_PATH = Path("data/orbitrap_hits_v2.csv")
OUT_PATH = Path("data/ik14_to_smiles_labeled.csv")
MISS_PATH = Path("data/ik14_missing_smiles.csv")


def canonicalize(smi: str) -> str | None:
    if not isinstance(smi, str) or not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def main() -> None:
    print(f"reading {FT_PATH} ...")
    ft = pd.read_csv(
        FT_PATH,
        usecols=["wiki_id", "library_wiki_id", "hit_ik14",
                 "anno_ik14", "spectrum_label"],
        dtype={"library_wiki_id": "string", "wiki_id": "string"},
    )
    labeled = ft[ft.spectrum_label != "blank"].copy()
    target_iks = set(labeled.hit_ik14.dropna().unique()) | \
                 set(labeled.anno_ik14.dropna().unique())
    print(f"  labeled rows: {len(labeled):,}")
    print(f"  target IK14s (hit_ik14 ∪ anno_ik14): {len(target_iks):,}")

    join_keys = labeled[["wiki_id", "library_wiki_id", "hit_ik14"]].dropna()
    join_keys = join_keys.drop_duplicates()
    print(f"  unique (wiki_id, library_wiki_id) join keys: {len(join_keys):,}")

    print(f"\nstreaming {OH_PATH} in chunks ...")
    ik14_to_smiles: dict[str, Counter] = {}
    chunk_iter = pd.read_csv(
        OH_PATH,
        usecols=["wiki_id", "library_wiki_id", "smiles"],
        dtype={"library_wiki_id": "string", "wiki_id": "string",
               "smiles": "string"},
        chunksize=200_000,
    )
    n_rows = 0
    for ci, chunk in enumerate(chunk_iter):
        n_rows += len(chunk)
        merged = chunk.merge(join_keys, on=["wiki_id", "library_wiki_id"],
                             how="inner")
        for ik14, smi in zip(merged.hit_ik14, merged.smiles):
            csmi = canonicalize(smi)
            if csmi is None:
                continue
            ik14_to_smiles.setdefault(ik14, Counter())[csmi] += 1
        print(f"  chunk {ci:>2d}: rows={n_rows:>10,d}   "
              f"ik14s_seen={len(ik14_to_smiles):>6,d}")

    print(f"\nresolving canonical SMILES per IK14 ...")
    rows: list[dict] = []
    for ik14, ctr in ik14_to_smiles.items():
        smi, n = ctr.most_common(1)[0]
        rows.append({"hit_ik14": ik14, "smiles": smi, "n_occurrences": n})
    df = pd.DataFrame(rows).sort_values("n_occurrences", ascending=False)

    df.to_csv(OUT_PATH, index=False)
    print(f"  wrote {OUT_PATH}  ({len(df):,} rows)")

    missing = sorted(target_iks - set(df.hit_ik14))
    pd.DataFrame({"hit_ik14": missing}).to_csv(MISS_PATH, index=False)
    print(f"  wrote {MISS_PATH}  ({len(missing):,} IK14s with no SMILES)")

    print(f"\nsummary:")
    print(f"  target IK14s          : {len(target_iks):>7,d}")
    print(f"  resolved with SMILES  : {len(df):>7,d}  "
          f"({len(df)/len(target_iks):.1%})")
    print(f"  missing               : {len(missing):>7,d}  "
          f"({len(missing)/len(target_iks):.1%})")


if __name__ == "__main__":
    main()
