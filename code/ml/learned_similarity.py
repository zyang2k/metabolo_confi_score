"""learned_similarity.py — MS2DeepScore-style Siamese encoder for evidence fusion.

Pretraining objective (mechanism 3 from project_lcms_neural_survey_pause):
  given a pair of MS2 spectra (A, B) of compounds with known SMILES,
  predict Tanimoto(Morgan_FP(A), Morgan_FP(B)) — a continuous chemistry-derived
  similarity. The label is FREE (computed from RDKit), so we get richer
  supervision than the binary TP/FP we've been training on.

Trained per-fold matching the GBM's GroupKFold(5) on anno_ik14, so the resulting
`learned_similarity` feature column is OOF and safe to feed into the GBM.

Pipeline:
  - Bin each spectrum into 2000 bins of 0.5 Da from 0 to 1000 m/z, sqrt-intensity.
  - Siamese encoder (MLP): 2000 → 1000 → 500 → 200, L2-normalized output.
  - Cosine similarity between embeddings → MSE against Tanimoto.
  - Per-fold inference: for every (query_wiki_id, library_wiki_id) row in
    feature_table_v2.csv where query is in the held-out fold, encode both
    spectra and write the cosine to data/learned_similarity_oof.csv.

Output:
  data/learned_similarity_oof.csv  (wiki_id, library_wiki_id, learned_similarity)

Usage:
  python code/learned_similarity.py            # full 5-fold OOF run
  python code/learned_similarity.py --smoke    # 1 fold, 5 epochs
"""

from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
from sklearn.model_selection import GroupKFold

RDLogger.DisableLog('rdApp.*')


# ── Config ───────────────────────────────────────────────────────────────────

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUERY_PEAKS = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')
LIB_PEAKS = os.path.join(ROOT, 'data', 'library_peaks_cache.json')
NEG_CUR = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv')
POS_CUR = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
OUT_OOF = os.path.join(ROOT, 'data', 'learned_similarity_oof.csv')

MZ_MIN, MZ_MAX = 0.0, 1000.0
BIN_WIDTH = 0.5
N_BINS = int((MZ_MAX - MZ_MIN) / BIN_WIDTH)   # 2000

EMB_DIM = 200
FP_RADIUS = 2
FP_NBITS = 2048


# ── Spectrum binning ─────────────────────────────────────────────────────────

def bin_spectrum(peaks: List[List[float]]) -> np.ndarray:
    """Bin peaks into a fixed-width vector with sqrt-intensity weighting."""
    v = np.zeros(N_BINS, dtype=np.float32)
    if not peaks:
        return v
    for mz, inten in peaks:
        if mz < MZ_MIN or mz >= MZ_MAX:
            continue
        idx = int((mz - MZ_MIN) / BIN_WIDTH)
        if 0 <= idx < N_BINS:
            v[idx] += float(inten)
    # Sqrt intensity, then L2 normalize
    v = np.sqrt(np.clip(v, 0.0, None))
    n = np.linalg.norm(v)
    if n > 0:
        v /= n
    return v


# ── Tanimoto on Morgan fingerprints ──────────────────────────────────────────

def smiles_to_fp(smiles: str):
    """Return Morgan ExplicitBitVect or None if unparseable."""
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_NBITS)


def tanimoto(fp_a, fp_b) -> float:
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


# ── Model ────────────────────────────────────────────────────────────────────

class SiameseEncoder(nn.Module):
    def __init__(self, n_bins: int = N_BINS, emb_dim: int = EMB_DIM, dropout: float = 0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_bins, 1000),
            nn.GELU(), nn.Dropout(dropout), nn.LayerNorm(1000),
            nn.Linear(1000, 500),
            nn.GELU(), nn.Dropout(dropout), nn.LayerNorm(500),
            nn.Linear(500, emb_dim),
        )

    def forward(self, x):
        z = self.net(x)
        return F.normalize(z, p=2, dim=-1)   # so cosine = dot product


def cosine_pair(za: torch.Tensor, zb: torch.Tensor) -> torch.Tensor:
    return (za * zb).sum(dim=-1)


# ── Pair dataset ─────────────────────────────────────────────────────────────

@dataclass
class SpectraIndex:
    """Holds binned spectra + per-spectrum metadata for fast pair sampling."""
    binned: torch.Tensor              # (N, N_BINS) float32
    ik14: np.ndarray                  # (N,) str
    fp_by_ik14: Dict[str, object]     # IK14 → Morgan FP
    spec_idx_by_ik14: Dict[str, List[int]]   # IK14 → list of spectrum row indices
    wiki_id: np.ndarray               # (N,) str

    def __len__(self):
        return self.binned.shape[0]


def build_train_index(peaks_cache: Dict[str, list],
                      anno_meta: pd.DataFrame) -> SpectraIndex:
    """Build SpectraIndex from labeled spectra with valid SMILES + cached peaks.

    anno_meta must have columns: wiki_id, anno_ik14, anno_smiles
    """
    rows = []
    fp_cache: Dict[str, object] = {}
    print(f'Building spectra index from {len(anno_meta):,} candidates...')
    n_no_peaks = 0
    n_no_fp = 0
    for _, r in anno_meta.iterrows():
        wid = str(r['wiki_id']).strip()
        ik = str(r.get('anno_ik14', '')).strip()
        smiles = str(r.get('anno_smiles', '')).strip()
        if not ik or not smiles:
            continue
        peaks = peaks_cache.get(wid)
        if peaks is None:
            n_no_peaks += 1
            continue
        if ik not in fp_cache:
            fp = smiles_to_fp(smiles)
            if fp is None:
                n_no_fp += 1
                continue
            fp_cache[ik] = fp
        v = bin_spectrum(peaks)
        rows.append((wid, ik, v))

    binned = np.stack([r[2] for r in rows], axis=0).astype(np.float32)
    wids = np.array([r[0] for r in rows])
    iks = np.array([r[1] for r in rows])
    spec_idx_by_ik14: Dict[str, List[int]] = {}
    for i, ik in enumerate(iks):
        spec_idx_by_ik14.setdefault(ik, []).append(i)
    print(f'  built: {len(rows):,} spectra, {len(spec_idx_by_ik14):,} unique IK14, '
          f'skipped {n_no_peaks:,} no-peaks + {n_no_fp:,} no-FP')
    return SpectraIndex(
        binned=torch.from_numpy(binned),
        ik14=iks, fp_by_ik14=fp_cache,
        spec_idx_by_ik14=spec_idx_by_ik14,
        wiki_id=wids,
    )


def sample_pair_batch(idx: SpectraIndex,
                      train_mask: np.ndarray,
                      batch_size: int,
                      same_compound_frac: float = 0.4,
                      rng: np.random.Generator = None):
    """Sample a batch of pairs balanced toward high-Tanimoto.

    same_compound_frac of the batch are pairs of spectra from the same IK14
    (Tanimoto = 1). The remaining are random pairs (mostly low Tanimoto).
    Returns (idx_a, idx_b, tanimoto) all length batch_size.
    """
    rng = rng or np.random.default_rng()
    n_same = int(batch_size * same_compound_frac)
    n_rand = batch_size - n_same

    # Allowed indices
    allowed = np.where(train_mask)[0]

    a_list, b_list, t_list = [], [], []

    # Same-compound pairs
    iks_with_multiple = [ik for ik, ids in idx.spec_idx_by_ik14.items()
                         if sum(1 for i in ids if train_mask[i]) >= 2]
    if iks_with_multiple:
        for _ in range(n_same):
            ik = rng.choice(iks_with_multiple)
            ids = [i for i in idx.spec_idx_by_ik14[ik] if train_mask[i]]
            ai, bi = rng.choice(ids, size=2, replace=False)
            a_list.append(ai); b_list.append(bi); t_list.append(1.0)

    # Random pairs (likely low Tanimoto)
    for _ in range(n_rand):
        ai, bi = rng.choice(allowed, size=2, replace=False)
        ik_a, ik_b = idx.ik14[ai], idx.ik14[bi]
        if ik_a == ik_b:
            t = 1.0
        else:
            t = tanimoto(idx.fp_by_ik14[ik_a], idx.fp_by_ik14[ik_b])
        a_list.append(ai); b_list.append(bi); t_list.append(t)

    return (np.array(a_list, dtype=np.int64),
            np.array(b_list, dtype=np.int64),
            np.array(t_list, dtype=np.float32))


# ── Training ─────────────────────────────────────────────────────────────────

def train_one_fold(idx: SpectraIndex,
                   train_mask: np.ndarray,
                   *,
                   n_epochs: int = 30,
                   batches_per_epoch: int = 200,
                   batch_size: int = 256,
                   lr: float = 3e-4,
                   wd: float = 1e-4,
                   seed: int = 42,
                   verbose: bool = True) -> SiameseEncoder:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = SiameseEncoder()
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_epochs)
    binned = idx.binned

    for epoch in range(n_epochs):
        model.train()
        loss_sum, n_seen = 0.0, 0
        for _ in range(batches_per_epoch):
            ai, bi, t = sample_pair_batch(idx, train_mask, batch_size, rng=rng)
            xa = binned[ai]
            xb = binned[bi]
            target = torch.from_numpy(t)

            opt.zero_grad()
            za = model(xa)
            zb = model(xb)
            sim = cosine_pair(za, zb)
            loss = F.mse_loss(sim, target)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            loss_sum += loss.item() * batch_size
            n_seen += batch_size
        sched.step()

        # Light per-epoch diagnostic
        if verbose:
            model.eval()
            with torch.no_grad():
                ai, bi, t = sample_pair_batch(idx, train_mask, 1024, rng=rng)
                za = model(binned[ai]); zb = model(binned[bi])
                sim = cosine_pair(za, zb).numpy()
                corr = float(np.corrcoef(sim, t)[0, 1])
            print(f'    epoch {epoch:3d}  loss={loss_sum/n_seen:.4f}  '
                  f'eval-pair corr(cos, Tanimoto)={corr:.4f}')
    return model


# ── Inference ────────────────────────────────────────────────────────────────

def encode_spectrum(model: SiameseEncoder, peaks: list) -> np.ndarray:
    """Encode a single spectrum to embedding."""
    v = bin_spectrum(peaks)
    with torch.no_grad():
        z = model(torch.from_numpy(v).unsqueeze(0))
    return z.squeeze(0).numpy()


def infer_fold(model: SiameseEncoder,
               test_wids: np.ndarray,
               feature_table: pd.DataFrame,
               query_cache: Dict[str, list],
               lib_cache: Dict[str, list]) -> pd.DataFrame:
    """For every (wiki_id, library_wiki_id) row in feature_table where
    wiki_id is in test_wids, compute learned_similarity = cos(encode(query), encode(library))."""
    sub = feature_table[feature_table['wiki_id'].isin(test_wids)][
        ['wiki_id', 'library_wiki_id']].dropna().drop_duplicates()
    print(f'    test rows to score: {len(sub):,}  (unique wiki_ids: {sub.wiki_id.nunique():,})')

    # Pre-encode unique queries and libraries to avoid re-encoding
    model.eval()
    q_emb: Dict[str, np.ndarray] = {}
    for wid in sub['wiki_id'].unique():
        peaks = query_cache.get(str(wid))
        if peaks is None:
            continue
        q_emb[str(wid)] = encode_spectrum(model, peaks)

    l_emb: Dict[str, np.ndarray] = {}
    for lid in sub['library_wiki_id'].unique():
        peaks = lib_cache.get(str(lid))
        if peaks is None:
            continue
        l_emb[str(lid)] = encode_spectrum(model, peaks)

    rows = []
    for _, r in sub.iterrows():
        wid = str(r['wiki_id']); lid = str(r['library_wiki_id'])
        if wid in q_emb and lid in l_emb:
            sim = float(np.dot(q_emb[wid], l_emb[lid]))
        else:
            sim = np.nan
        rows.append({'wiki_id': wid, 'library_wiki_id': lid, 'learned_similarity': sim})
    return pd.DataFrame(rows)


# ── Main ─────────────────────────────────────────────────────────────────────

def main(*, smoke: bool = False, n_epochs: int = 30,
         batches_per_epoch: int = 200, batch_size: int = 256,
         seed: int = 42):
    if smoke:
        n_epochs = 5
        batches_per_epoch = 50

    print(f'[1/4] Loading spectra metadata...')
    # Build anno_meta from curated CSVs (wiki_id, anno_ik14, anno_smiles)
    rows = []
    for path in [NEG_CUR, POS_CUR]:
        df = pd.read_csv(path, low_memory=False)
        # Annotated spectra only
        df = df[df['is_manual_annotated'].astype(bool) & df['name'].notna()
                & (df['name'].astype(str).str.strip() != '')]
        df['anno_smiles'] = df['annotation-smiles'].fillna(df['smiles'])
        rows.append(df[['wiki_id', 'anno_smiles']])
    cur = pd.concat(rows, ignore_index=True)

    # Pull anno_ik14 from feature_table for consistency with GBM splits
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False,
                     usecols=['wiki_id', 'spectrum_label', 'anno_ik14',
                              'library_wiki_id'])
    labeled = (ft[ft['spectrum_label'].isin(['TP', 'FP'])]
               .drop_duplicates('wiki_id')[['wiki_id', 'anno_ik14']])
    anno_meta = labeled.merge(cur, on='wiki_id', how='left')
    anno_meta = anno_meta.dropna(subset=['anno_smiles'])
    print(f'  anno_meta: {len(anno_meta):,} labeled spectra with SMILES')

    print(f'[2/4] Loading peak caches...')
    with open(QUERY_PEAKS) as f:
        query_cache = json.load(f)
    with open(LIB_PEAKS) as f:
        lib_cache = json.load(f)
    print(f'  query peaks: {len(query_cache):,}  library peaks: {len(lib_cache):,}')

    spec_idx = build_train_index(query_cache, anno_meta)

    print(f'[3/4] 5-fold GroupKFold OOF training (n_epochs={n_epochs})...')
    # Build group key + index alignment to ft labeled rows
    wid_to_pos = {w: i for i, w in enumerate(spec_idx.wiki_id)}
    aligned = anno_meta[anno_meta['wiki_id'].astype(str).isin(wid_to_pos)].copy()
    aligned['_pos'] = aligned['wiki_id'].astype(str).map(wid_to_pos)
    aligned = aligned.dropna(subset=['_pos']).drop_duplicates('wiki_id')

    groups = aligned['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    spec_positions = aligned['_pos'].astype(int).values
    n_total = len(spec_idx)
    gkf = GroupKFold(n_splits=5)
    splits = list(gkf.split(spec_positions, groups=groups))
    if smoke:
        splits = splits[:1]

    all_results: List[pd.DataFrame] = []
    for fold, (tr_pos_idx, te_pos_idx) in enumerate(splits):
        t0 = time.time()
        train_mask = np.zeros(n_total, dtype=bool)
        train_mask[spec_positions[tr_pos_idx]] = True
        test_wids = aligned.iloc[te_pos_idx]['wiki_id'].astype(str).values

        print(f'\n=== Fold {fold} ===  train={train_mask.sum()}  test={len(test_wids)}')
        model = train_one_fold(
            spec_idx, train_mask,
            n_epochs=n_epochs, batches_per_epoch=batches_per_epoch,
            batch_size=batch_size, seed=seed + fold,
        )
        # Inference on test fold's full feature_table rows
        ft_full = pd.read_csv(FEATURE_TABLE, low_memory=False,
                              usecols=['wiki_id', 'library_wiki_id', 'spectrum_label'])
        ft_full = ft_full[ft_full['spectrum_label'].isin(['TP', 'FP'])]
        df_fold = infer_fold(model, test_wids, ft_full, query_cache, lib_cache)
        df_fold['fold'] = fold
        all_results.append(df_fold)
        print(f'  fold time: {time.time() - t0:.1f}s')

    out = pd.concat(all_results, ignore_index=True)
    out.to_csv(OUT_OOF, index=False)
    print(f'\nWrote {OUT_OOF}: {len(out):,} rows  '
          f'(unique wiki_ids: {out["wiki_id"].nunique():,})')

    # Quick summary
    valid = out['learned_similarity'].dropna()
    print(f'  learned_similarity: median={valid.median():.4f}  '
          f'p05={valid.quantile(0.05):.4f}  p95={valid.quantile(0.95):.4f}  '
          f'NaN={out["learned_similarity"].isna().sum()}')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--smoke', action='store_true', help='1 fold, 5 epochs')
    p.add_argument('--epochs', type=int, default=30)
    p.add_argument('--batches', type=int, default=200)
    p.add_argument('--batch_size', type=int, default=256)
    args = p.parse_args()
    main(smoke=args.smoke, n_epochs=args.epochs,
         batches_per_epoch=args.batches, batch_size=args.batch_size)
