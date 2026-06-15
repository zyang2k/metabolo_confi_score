"""bench_label_propagate.py — Route 2: graph-native label-propagation score.

The honest realization of "convert the molecular network into a score" (vs Route 1 =
relational features into the GBM). Build a spectral-similarity graph over the labeled
bins, seed curated TP=1 / FP=0, diffuse, read each bin's soft label, calibrate.

Graph: nodes = labeled top-1 bins; edges = entropy similarity >= SIM_EDGE between bins
in the same polarity within precursor +-10 mDa, weight = similarity. (This is the GNPS
edge; we do NOT precursor-lock to same-compound only — analog edges carry label mass,
which is the point of propagation.)

Propagation (personalized / PageRank-style):  f <- a * S f + (1-a) * y0,
S = row-normalized W, y0 = seed (train labels, test = 0), ~50 iters, a = 0.85.

LEAKAGE: GroupKFold(anno_ik14). In each fold only TRAINING bins are seeded; test bins
(and their same-IK14 siblings, kept together by GroupKFold) carry no seed, so a bin can
never propagate its own held-out label. OOF score = test-fold f. Then isotonic-calibrate.

Reports: graph coverage, OOF AUC vs GBM 0.9131, Brier/ECE (raw+cal), and AUC restricted
to the connected subset (where the graph actually carries signal).
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import ms_entropy as me
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bench_harness import FEATURE_TABLE, build_top1, make_groups, ece

ROOT = Path(__file__).resolve().parent.parent
QP = ROOT / 'data' / 'query_peaks_cache_v2.json'
MS1_TOL, MS2_TOL = 0.01, 0.02
SIM_EDGE = 0.50       # edge threshold (lower than 0.70 to give propagation reach)
ALPHA, N_ITER = 0.85, 50


def build_edges(top1, peaks):
    """Symmetric weighted edge list (i, j, sim) among labeled bins."""
    t = top1.reset_index(drop=True)
    mz = pd.to_numeric(t['precursor_mz'], errors='coerce').values
    wid = t['wiki_id'].values
    pol = t['polarity'].astype(str).values
    n = len(t)
    edges = []
    for p in np.unique(pol):
        idx = np.where(pol == p)[0]
        order = idx[np.argsort(mz[idx])]
        mzs = mz[order]
        for a in range(len(order)):
            i = order[a]
            if wid[i] not in peaks or not np.isfinite(mz[i]):
                continue
            qp = np.asarray(peaks[wid[i]], np.float64)
            b = a + 1
            while b < len(order) and mzs[b] - mzs[a] <= MS1_TOL:
                j = order[b]
                if wid[j] in peaks:
                    s = me.calculate_entropy_similarity(
                        qp, np.asarray(peaks[wid[j]], np.float64),
                        ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
                    if s >= SIM_EDGE:
                        edges.append((i, j, float(s)))
                b += 1
    return edges, n


def propagate(edges, n, seed, seed_mask):
    """f <- a S f + (1-a) y0, S row-normalized over the weighted graph."""
    deg = np.zeros(n)
    for i, j, w in edges:
        deg[i] += w; deg[j] += w
    deg[deg == 0] = 1.0
    y0 = np.where(seed_mask, seed, 0.0)
    f = y0.copy()
    for _ in range(N_ITER):
        agg = np.zeros(n)
        for i, j, w in edges:
            agg[i] += w * f[j]
            agg[j] += w * f[i]
        f = ALPHA * (agg / deg) + (1 - ALPHA) * y0
    return f


def main():
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    t = build_top1(ft).reset_index(drop=True)
    peaks = {k: v for k, v in json.load(open(QP)).items() if v}
    print(f'labeled top-1 bins: {len(t):,}')

    edges, n = build_edges(t, peaks)
    deg = np.zeros(n)
    for i, j, w in edges:
        deg[i] += 1; deg[j] += 1
    cov = (deg > 0).mean()
    print(f'edges (sim>={SIM_EDGE}): {len(edges):,}  |  nodes with >=1 edge: {cov:.1%}')

    y = t['hit_label'].values.astype(float)
    groups = make_groups(t)
    oof = np.full(n, np.nan)
    gkf = GroupKFold(n_splits=5)
    for tr, te in gkf.split(t, y, groups):
        seed_mask = np.zeros(n, bool); seed_mask[tr] = True
        f = propagate(edges, n, y, seed_mask)
        oof[te] = f[te]

    auc = roc_auc_score(y, oof)
    iso = IsotonicRegression(out_of_bounds='clip'); iso.fit(oof, y)
    oc = iso.transform(oof)
    print('\n=== Route 2: graph-native propagated score (OOF) ===')
    print(f'  AUC           {auc:.4f}   (GBM baseline 0.9131)')
    print(f'  Brier raw/cal {brier_score_loss(y, oof):.4f} / {brier_score_loss(y, oc):.4f}')
    print(f'  ECE   raw/cal {ece(oof, y):.4f} / {ece(oc, y):.4f}')

    # restricted to connected nodes (where the graph carries signal)
    conn = deg > 0
    if conn.sum() > 30 and len(np.unique(y[conn])) == 2:
        print(f'\n  connected subset ({conn.sum():,} bins, '
              f'TP rate {y[conn].mean():.3f}):')
        print(f'    AUC {roc_auc_score(y[conn], oof[conn]):.4f}')
        iso2 = IsotonicRegression(out_of_bounds='clip'); iso2.fit(oof[conn], y[conn])
        print(f'    Brier cal {brier_score_loss(y[conn], iso2.transform(oof[conn])):.4f}')
    iso_iso = oof[~conn]
    print(f'  isolated nodes: {(~conn).sum():,} (score=0, no graph signal)')

    pd.DataFrame({'wiki_id': t['wiki_id'], 'prop_score': oof,
                  'prop_cal': oc, 'degree': deg, 'hit_label': y}).to_csv(
        ROOT / 'data' / 'label_propagate_scores.csv', index=False)
    print('\nWrote data/label_propagate_scores.csv')


if __name__ == '__main__':
    main()
