"""bench_ttof_relational.py — Cross-platform check for the relational KG features.

Auditor gate (project_relational_kg_features_20260604): the +0.0070 lift from
rel_nn_sim / rel_n_confirmed_nbr was measured on Orbitrap only. A confusability graph
built from curated Orbitrap bins may not transfer to TTOF. This bench:

  1. Reuses build_ttof_table (golden TP / yy_ FP labels) from bench_ttof_golden.
  2. Computes the relational features WITHIN TTOF (each TTOF labeled top-1 bin searched
     against the other TTOF labeled bins — the deployment-realistic in-platform graph),
     same leakage guard (neighbours blocked from sharing anno_ik14).
  3. Trains the Orbitrap model on BASE_NUMERIC (mirrors the +0.0070 bench — no struct_logit)
     with vs without the 2 relational features, predicts on TTOF, reports per-polarity AUC.

Pass = relational delta on TTOF is >= 0 and same sign/direction as Orbitrap (FP more
confusable). Fail = inverts or vanishes -> recast as triage flag, do not ship as feature.
"""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import ms_entropy as me
import xgboost as xgb
from sklearn.metrics import roc_auc_score, brier_score_loss

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bench_harness import BASE_NUMERIC, CATEGORICAL, FEATURE_TABLE, build_top1, prep, train
from bench_ttof_golden import (
    build_ttof_table, load_adduct_taxonomy, LIB_PEAKS,
    TTOF_NEG_HITS, TTOF_POS_HITS, TTOF_NEG_SPEC, TTOF_POS_SPEC,
    TTOF_NEG_PEAKS, TTOF_POS_PEAKS,
)

ROOT = Path(__file__).resolve().parent.parent
REL = ['rel_nn_sim', 'rel_n_confirmed_nbr']
MS1_TOL, MS2_TOL, SIM_THRESH = 0.01, 0.02, 0.70


def compute_relational(top1: pd.DataFrame, peaks: dict) -> pd.DataFrame:
    """Within-platform confusability features on a labeled top-1 frame.

    Same definition + leakage guard as build_relational_features.py: neighbours blocked
    from sharing anno_ik14 (empty -> own block via wiki_id). Returns wiki_id + REL cols.
    """
    t = top1.copy()
    t['precursor_mz'] = pd.to_numeric(t['precursor_mz'], errors='coerce')
    t['anno_ik14'] = t['anno_ik14'].fillna('').astype(str)
    t['wiki_id'] = t['wiki_id'].astype(str)
    t['block'] = np.where(t['anno_ik14'] != '', t['anno_ik14'], 'self::' + t['wiki_id'])
    rows = []
    for pol, sub in t.groupby('polarity'):
        sub = sub.reset_index(drop=True)
        mz = sub['precursor_mz'].values
        order = np.argsort(mz)
        mz_sorted = mz[order]
        wid = sub['wiki_id'].values
        blk = sub['block'].values
        lab = sub['hit_label'].values.astype(float)
        have = np.array([w in peaks for w in wid])
        for i in range(len(sub)):
            if not have[i] or not np.isfinite(mz[i]):
                rows.append((wid[i], 0.0, 0)); continue
            lo = np.searchsorted(mz_sorted, mz[i] - MS1_TOL, 'left')
            hi = np.searchsorted(mz_sorted, mz[i] + MS1_TOL, 'right')
            cand = [order[p] for p in range(lo, hi)
                    if order[p] != i and blk[order[p]] != blk[i] and have[order[p]]]
            qp = np.asarray(peaks[wid[i]], dtype=np.float64)
            nn_sim, n_conf = 0.0, 0
            for j in cand:
                s = me.calculate_entropy_similarity(
                    qp, np.asarray(peaks[wid[j]], dtype=np.float64),
                    ms2_tolerance_in_da=MS2_TOL, clean_spectra=True)
                if s > nn_sim:
                    nn_sim = s
                if s >= SIM_THRESH and lab[j] == 1:
                    n_conf += 1
            rows.append((wid[i], float(nn_sim), int(n_conf)))
    return pd.DataFrame(rows, columns=['wiki_id'] + REL)


def predict(model, frame, extra, n_num, n_cat):
    X = prep(frame, extra)
    d = xgb.DMatrix(X, enable_categorical=True,
                    feature_types=['q'] * n_num + ['c'] * n_cat)
    return model.predict(d)


def main():
    tax = load_adduct_taxonomy()
    print('Loading library peaks cache (320 MB) ...')
    lib_peaks = json.load(open(LIB_PEAKS))

    neg = build_ttof_table(TTOF_NEG_HITS, TTOF_NEG_SPEC, TTOF_NEG_PEAKS, lib_peaks, 0, tax)
    pos = build_ttof_table(TTOF_POS_HITS, TTOF_POS_SPEC, TTOF_POS_PEAKS, lib_peaks, 1, tax)
    ttof = pd.concat([neg, pos], ignore_index=True)
    ttof_top1 = build_top1(ttof).reset_index(drop=True)
    print(f'\nTTOF top-1 labeled: TP={(ttof_top1.spectrum_label=="TP").sum():,} '
          f'FP={(ttof_top1.spectrum_label=="FP").sum():,}')

    # within-TTOF relational features
    print('Computing within-TTOF relational features ...')
    tp = {**json.load(open(TTOF_NEG_PEAKS)), **json.load(open(TTOF_POS_PEAKS))}
    relt = compute_relational(ttof_top1, tp)
    ttof_top1['wiki_id'] = ttof_top1['wiki_id'].astype(str)
    ttof_top1 = ttof_top1.merge(relt, on='wiki_id', how='left')
    # direction check on TTOF
    for c in REL:
        gTP = ttof_top1.loc[ttof_top1.hit_label == 1, c].mean()
        gFP = ttof_top1.loc[ttof_top1.hit_label == 0, c].mean()
        print(f'  {c:22s} TP {gTP:.3f}  FP {gFP:.3f}  (FP>TP confusability: {gFP>gTP})')

    # Orbitrap training top-1 + relational (from the shipped Orbitrap side table)
    orb = pd.read_csv(FEATURE_TABLE, low_memory=False)
    orb_top1 = build_top1(orb).reset_index(drop=True)
    orb_top1['wiki_id'] = orb_top1['wiki_id'].astype(str)
    rel_orb = pd.read_csv(ROOT / 'data' / 'relational_features.csv')
    rel_orb['wiki_id'] = rel_orb['wiki_id'].astype(str)
    orb_top1 = orb_top1.merge(rel_orb[['wiki_id'] + REL], on='wiki_id', how='left')

    y_tr = orb_top1['hit_label'].values
    n_cat = len(CATEGORICAL)
    results = {}
    for arm, extra in [('baseline', []), ('+relational', REL)]:
        n_num = len(BASE_NUMERIC) + len(extra)
        model = train(prep(orb_top1, extra), y_tr, n_num, n_cat)
        ttof_top1[f'pred_{arm}'] = predict(model, ttof_top1, extra, n_num, n_cat)
        results[arm] = model

    # per-polarity + combined AUC for each arm
    print('\n' + '=' * 64)
    print('TTOF cross-platform AUC  (Orbitrap-trained, BASE_NUMERIC)')
    print('=' * 64)
    rows = []
    lab = ttof_top1.spectrum_label.isin(['TP', 'FP'])
    for pol_val, pol in [(0, 'neg'), (1, 'pos'), (None, 'all')]:
        m = lab if pol_val is None else (lab & (ttof_top1.polarity == pol_val))
        y = ttof_top1.loc[m, 'hit_label'].values
        if len(np.unique(y)) < 2:
            continue
        ab = roc_auc_score(y, ttof_top1.loc[m, 'pred_baseline'].values)
        ar = roc_auc_score(y, ttof_top1.loc[m, 'pred_+relational'].values)
        rows.append(dict(polarity=pol, n=len(y), n_TP=int((y == 1).sum()),
                         AUC_base=ab, AUC_rel=ar, delta=ar - ab))
        print(f'  {pol:4s} n={len(y):4d}  base {ab:.4f}  +rel {ar:.4f}  Δ {ar-ab:+.4f}')
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / 'data' / 'bench_ttof_relational_summary.csv', index=False)
    print(f"\nWrote data/bench_ttof_relational_summary.csv")
    print('\nReference: Orbitrap in-platform Δ(+relational) = +0.0070 (0.9131→0.9201)')


if __name__ == '__main__':
    main()
