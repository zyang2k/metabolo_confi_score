"""bench_st_hybrid.py — Head C: GBM stack on Set Transformer OOF embeddings.

Pipeline:
  1. Build the same top-1-by-esim training set as score_gbm_v2.py
  2. Run a baseline GBM on original features alone → record OOF AUC
  3. Join in OOF (z_set, z_top, set_logit, top1_logit) from data/st_oof_embeddings.npz
  4. Run a hybrid GBM on [original features + ST embeddings] on identical splits
  5. Report Δ AUC; both with isotonic calibration

The OOF embeddings are properly held-out (each spectrum's embedding was generated
by an ST model that didn't see that spectrum's anno_ik14 group during training),
so feeding them into a downstream GBM with the same GroupKFold splits is a clean
stacking pattern.

Outputs:
  data/bench_st_hybrid_summary.csv   per-fold and aggregate AUCs for both arms
"""

from __future__ import annotations

import os
import warnings
from typing import List

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
ST_EMB = os.path.join(ROOT, 'data', 'st_oof_embeddings.npz')
SUMMARY_OUT = os.path.join(ROOT, 'data', 'bench_st_hybrid_summary.csv')

# Mirrors score_gbm_v2.py exactly — same NUMERIC + CATEGORICAL split.
NUMERIC_FEATURES = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEGORICAL_FEATURES = ['hit_adduct_cat', 'db', 'polarity']


def prep_features(df: pd.DataFrame, extra_numeric: List[str]) -> pd.DataFrame:
    X = df.copy()
    for c in NUMERIC_FEATURES:
        if c not in X.columns:
            X[c] = np.nan
        X[c] = pd.to_numeric(X[c], errors='coerce')
    for c in CATEGORICAL_FEATURES:
        s = X[c].astype(str).fillna('missing')
        cats = sorted(s.unique().tolist())
        idx = {v: i for i, v in enumerate(cats)}
        X[c] = s.map(idx).astype('int32')
    cols = NUMERIC_FEATURES + extra_numeric + CATEGORICAL_FEATURES
    return X[cols]


def train_xgb(X_train: pd.DataFrame, y_train: np.ndarray, n_num: int, n_cat: int):
    dtrain = xgb.DMatrix(X_train, label=y_train,
                         enable_categorical=True,
                         feature_types=['q'] * n_num + ['c'] * n_cat)
    params = {
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'tree_method': 'hist',
        'max_depth': 5,
        'learning_rate': 0.05,
        'subsample': 0.85,
        'colsample_bytree': 0.85,
        'min_child_weight': 5,
        'reg_alpha': 0.1,
        'reg_lambda': 1.0,
        'seed': 42,
        'verbosity': 0,
    }
    return xgb.train(params, dtrain, num_boost_round=500)


def ece(p: np.ndarray, y: np.ndarray, nbins: int = 10) -> float:
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0:
            continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / max(n, 1)


def run_oof(top1: pd.DataFrame, extra_numeric: List[str], groups: np.ndarray,
            tag: str) -> np.ndarray:
    """Return OOF predictions on top1 rows. Same 5-fold GroupKFold as score_gbm_v2.py."""
    X_df = prep_features(top1, extra_numeric)
    n_num = len(NUMERIC_FEATURES) + len(extra_numeric)
    n_cat = len(CATEGORICAL_FEATURES)
    y = top1['hit_label'].values

    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        model = train_xgb(X_df.iloc[tr], y[tr], n_num, n_cat)
        dte = xgb.DMatrix(X_df.iloc[te],
                          enable_categorical=True,
                          feature_types=['q'] * n_num + ['c'] * n_cat)
        oof[te] = model.predict(dte)
        fold_auc = roc_auc_score(y[te], oof[te])
        fold_aucs.append(fold_auc)
        print(f'  [{tag}] Fold {fold}: AUC={fold_auc:.4f}')
    auc_oof = roc_auc_score(y, oof)
    print(f'  [{tag}] OOF AUC: {auc_oof:.4f}   per-fold: '
          f'[{", ".join(f"{a:.4f}" for a in fold_aucs)}]')
    return oof, fold_aucs, auc_oof


def main():
    print(f'Loading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1_idx = labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    print(f'  {len(top1):,} top-1 rows; prior_TP={top1["hit_label"].mean():.3f}')

    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    print(f'\nLoading ST OOF embeddings: {ST_EMB}')
    st = np.load(ST_EMB, allow_pickle=True)
    st_df = pd.DataFrame({
        'wiki_id': st['wiki_id'],
        'set_logit': st['set_logit'],
        'top1_logit': st['top1_logit'],
    })
    z_set = pd.DataFrame(st['z_set'], columns=[f'z_set_{i}' for i in range(st['z_set'].shape[1])])
    z_top = pd.DataFrame(st['z_top'], columns=[f'z_top_{i}' for i in range(st['z_top'].shape[1])])
    z_set['wiki_id'] = st['wiki_id']
    z_top['wiki_id'] = st['wiki_id']

    n_emb = z_set.shape[1] - 1   # minus wiki_id
    print(f'  embeddings: z_set ({n_emb}-d), z_top ({n_emb}-d), 2 logits')

    # Merge ST features into top1 (left-join; spectra without ST embedding stay NaN)
    top1 = top1.merge(st_df, on='wiki_id', how='left')
    top1 = top1.merge(z_set, on='wiki_id', how='left')
    top1 = top1.merge(z_top, on='wiki_id', how='left')
    coverage = top1['set_logit'].notna().mean()
    print(f'  ST embedding coverage: {coverage*100:.1f}% of top-1 rows')

    extra_baseline: List[str] = []
    extra_logits = ['set_logit', 'top1_logit']
    extra_zset = [f'z_set_{i}' for i in range(n_emb)]
    extra_ztop = [f'z_top_{i}' for i in range(n_emb)]

    print('\n=== Baseline GBM (no ST features) ===')
    oof_b, folds_b, auc_b = run_oof(top1, extra_baseline, groups, tag='baseline')

    print('\n=== Hybrid GBM (+ ST logits) ===')
    oof_l, folds_l, auc_l = run_oof(top1, extra_logits, groups, tag='+logits')

    print('\n=== Hybrid GBM (+ ST logits + z_set) ===')
    oof_zs, folds_zs, auc_zs = run_oof(top1, extra_logits + extra_zset, groups, tag='+logits+z_set')

    print('\n=== Hybrid GBM (+ ST logits + z_set + z_top) ===')
    oof_full, folds_full, auc_full = run_oof(top1, extra_logits + extra_zset + extra_ztop,
                                              groups, tag='+full')

    # Summary
    rows = []
    for name, oof, folds, auc in [('baseline', oof_b, folds_b, auc_b),
                                   ('+st_logits', oof_l, folds_l, auc_l),
                                   ('+st_logits+z_set', oof_zs, folds_zs, auc_zs),
                                   ('+st_logits+z_set+z_top', oof_full, folds_full, auc_full)]:
        y = top1['hit_label'].values
        iso = IsotonicRegression(out_of_bounds='clip')
        iso.fit(oof, y)
        oof_cal = iso.transform(oof)
        rows.append({
            'arm': name,
            'auc_oof': auc,
            'auc_fold0': folds[0], 'auc_fold1': folds[1], 'auc_fold2': folds[2],
            'auc_fold3': folds[3], 'auc_fold4': folds[4],
            'brier_raw': brier_score_loss(y, oof),
            'brier_cal': brier_score_loss(y, oof_cal),
            'ece_raw': ece(oof, y),
            'ece_cal': ece(oof_cal, y),
            'delta_vs_baseline': auc - auc_b,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_OUT, index=False)

    print(f'\n=== Final summary ===')
    cols_show = ['arm', 'auc_oof', 'delta_vs_baseline', 'brier_cal', 'ece_cal']
    print(summary[cols_show].to_string(index=False, float_format=lambda x: f'{x:.4f}'))
    print(f'\nWrote {SUMMARY_OUT}')


if __name__ == '__main__':
    main()
