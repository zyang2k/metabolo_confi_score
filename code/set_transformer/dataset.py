"""Build the (B, K, d) candidate-set tensor dataset from feature_table_v2.csv.

Per spectrum (wiki_id), we collect every candidate hit, sort by entropy_similarity
desc, truncate to K_max, and pad shorter sets to K_max with zeros.

Per-hit features are deliberately chosen to NOT include sim_gap or n_candidates —
those are exactly the set-level summaries we want the Set Transformer to learn.
Other set-level GBM features (`compound_has_ok_adduct`, `n_candidate_adducts`)
are kept because they encode multi-library reinforcement at the COMPOUND level,
which is not redundant with cross-candidate attention.

Outputs:
    SetDataset.X       (n_spectra, K, d_in)   per-hit features (zero-padded)
    SetDataset.M       (n_spectra, K)         True for valid candidates
    SetDataset.y_hit   (n_spectra, K)         1=TP candidate, 0=FP, -1=pad
    SetDataset.y_set   (n_spectra,)           1 if any TP candidate, else 0
    SetDataset.groups  (n_spectra,)           anno_ik14 string for GroupKFold
    SetDataset.wiki_id (n_spectra,)           bin id, for output joins
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset


# ── Feature schema ────────────────────────────────────────────────────────────
# Per-hit numeric features. Kept aligned with score_gbm_v2.py minus the set-summary
# features (sim_gap, n_candidates) that the Set Transformer is meant to learn.
NUMERIC_FEATURES: List[str] = [
    'entropy_similarity', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy',
    'n_candidate_adducts', 'compound_has_ok_adduct',
    'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]

# Categoricals are integer-coded; the model embeds them.
CATEGORICAL_FEATURES: List[str] = ['hit_adduct_cat', 'db', 'polarity']


@dataclass
class SetDataset(Dataset):
    """In-memory tensor dataset for the Set Transformer."""

    X: torch.Tensor                  # (N, K, d_num)
    C: torch.Tensor                  # (N, K, n_cat)  integer-coded categoricals
    M: torch.Tensor                  # (N, K) bool
    y_hit: torch.Tensor              # (N, K) long; -1 for padding
    y_set: torch.Tensor              # (N,)   float
    groups: np.ndarray               # (N,)   str
    wiki_id: np.ndarray              # (N,)   str
    cat_vocab: Dict[str, List[str]] = field(default_factory=dict)
    feature_names: List[str] = field(default_factory=list)

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        return {
            'X': self.X[idx],
            'C': self.C[idx],
            'M': self.M[idx],
            'y_hit': self.y_hit[idx],
            'y_set': self.y_set[idx],
        }

    @property
    def d_num(self) -> int:
        return self.X.shape[-1]

    @property
    def cat_cardinalities(self) -> List[int]:
        return [len(self.cat_vocab[c]) for c in CATEGORICAL_FEATURES]


def _encode_categoricals(df: pd.DataFrame) -> Tuple[np.ndarray, Dict[str, List[str]]]:
    """Integer-encode each categorical column. Reserved index 0 = '<pad/missing>'."""
    out = np.zeros((len(df), len(CATEGORICAL_FEATURES)), dtype=np.int64)
    vocab: Dict[str, List[str]] = {}
    for j, c in enumerate(CATEGORICAL_FEATURES):
        s = df[c].astype(str).fillna('missing')
        cats = ['<pad/missing>'] + sorted(set(s.unique()) - {'<pad/missing>'})
        idx = {v: i for i, v in enumerate(cats)}
        out[:, j] = s.map(lambda v: idx.get(v, 0)).values
        vocab[c] = cats
    return out, vocab


def build_dataset(feature_table_csv: str,
                  K_max: int = 25,
                  labeled_only: bool = True) -> SetDataset:
    """Construct the set tensor dataset.

    Args:
        feature_table_csv: path to data/feature_table_v2.csv
        K_max: max set size per spectrum; sort by esim desc and truncate
        labeled_only: keep only spectra with spectrum_label in {TP, FP}

    Returns:
        SetDataset ready for DataLoader.
    """
    df = pd.read_csv(feature_table_csv, low_memory=False)
    if labeled_only:
        df = df[df['spectrum_label'].isin(['TP', 'FP'])].reset_index(drop=True)

    # Coerce numerics
    for c in NUMERIC_FEATURES:
        if c not in df.columns:
            df[c] = np.nan
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # Encode categoricals once over the whole table for a stable vocabulary
    cat_codes, cat_vocab = _encode_categoricals(df)
    df['_cat_codes_idx'] = np.arange(len(df))   # row index for fast gather

    # Per-spectrum sort by esim desc, truncate to K_max
    df = df.sort_values(['wiki_id', 'entropy_similarity'],
                        ascending=[True, False]).reset_index(drop=True)
    df['_within_rank'] = df.groupby('wiki_id').cumcount()
    df = df[df['_within_rank'] < K_max].reset_index(drop=True)
    cat_codes = cat_codes[df['_cat_codes_idx'].values]   # realign after sort/filter

    # Group meta (one row per spectrum)
    spec_meta = (
        df.groupby('wiki_id', sort=False)
          .agg(anno_ik14=('anno_ik14', 'first'),
               spectrum_label=('spectrum_label', 'first'))
          .reset_index()
    )
    n_spec = len(spec_meta)
    wiki_to_idx = {w: i for i, w in enumerate(spec_meta['wiki_id'].values)}

    d_num = len(NUMERIC_FEATURES)
    X = np.zeros((n_spec, K_max, d_num), dtype=np.float32)
    C = np.zeros((n_spec, K_max, len(CATEGORICAL_FEATURES)), dtype=np.int64)
    M = np.zeros((n_spec, K_max), dtype=bool)
    y_hit = np.full((n_spec, K_max), -1, dtype=np.int64)

    spec_idx = df['wiki_id'].map(wiki_to_idx).values
    rank = df['_within_rank'].values
    Xnum = df[NUMERIC_FEATURES].astype(np.float32).values
    # Zero-fill NaNs (the per-hit numerics are bounded; missing → neutral 0).
    Xnum = np.nan_to_num(Xnum, nan=0.0, posinf=0.0, neginf=0.0)
    yh = df['hit_label'].astype(np.int64).values

    X[spec_idx, rank] = Xnum
    C[spec_idx, rank] = cat_codes
    M[spec_idx, rank] = True
    y_hit[spec_idx, rank] = yh

    # Set-level label: 1 if any candidate is TP, else 0
    y_set = (y_hit == 1).any(axis=1).astype(np.float32)

    # Group key for GroupKFold; empty → unique placeholder so it doesn't pool spectra
    groups = spec_meta['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    return SetDataset(
        X=torch.from_numpy(X),
        C=torch.from_numpy(C),
        M=torch.from_numpy(M),
        y_hit=torch.from_numpy(y_hit),
        y_set=torch.from_numpy(y_set),
        groups=groups,
        wiki_id=spec_meta['wiki_id'].values,
        cat_vocab=cat_vocab,
        feature_names=NUMERIC_FEATURES,
    )


def collate_passthrough(batch):
    """Default collate that stacks dict tensors."""
    out = {}
    for k in batch[0]:
        out[k] = torch.stack([b[k] for b in batch], dim=0)
    return out


__all__ = ['SetDataset', 'build_dataset', 'collate_passthrough',
           'NUMERIC_FEATURES', 'CATEGORICAL_FEATURES']
