"""bench_noise_features.py — Test Oliver's 5 noise features as GBM additions.

Two arms on identical 5-fold GroupKFold splits, both isotonic-calibrated:
  baseline      : current 18-feature GBM (matches score_gbm_v2.py)
  +noise        : baseline + n_ions, normalized_entropy, top5_pct,
                                ion_density, median_to_base

Prerequisite: data/feature_table_v2.csv must contain the 5 noise feature
columns. Run code/add_noise_features.py first.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
SUMMARY_OUT = os.path.join(ROOT, 'data', 'bench_noise_features_summary.csv')

# Same 18 features as score_gbm_v2.py
BASE_NUMERIC = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEGORICAL = ['hit_adduct_cat', 'db', 'polarity']

NOISE_FEATURES = ['n_ions', 'normalized_entropy', 'top5_pct',
                  'ion_density', 'median_to_base']


def prep(df: pd.DataFrame, extra_num: list[str]) -> pd.DataFrame:
    X = df.copy()
    for c in BASE_NUMERIC + extra_num:
        if c not in X.columns:
            X[c] = np.nan
        X[c] = pd.to_numeric(X[c], errors='coerce')
    for c in CATEGORICAL:
        s = X[c].astype(str).fillna('missing')
        cats = sorted(s.unique().tolist())
        idx = {v: i for i, v in enumerate(cats)}
        X[c] = s.map(idx).astype('int32')
    return X[BASE_NUMERIC + extra_num + CATEGORICAL]


def train(X, y, n_num, n_cat):
    dtr = xgb.DMatrix(X, label=y, enable_categorical=True,
                      feature_types=['q'] * n_num + ['c'] * n_cat)
    params = {'objective': 'binary:logistic', 'eval_metric': 'auc',
              'tree_method': 'hist', 'max_depth': 5, 'learning_rate': 0.05,
              'subsample': 0.85, 'colsample_bytree': 0.85, 'min_child_weight': 5,
              'reg_alpha': 0.1, 'reg_lambda': 1.0, 'seed': 42, 'verbosity': 0}
    return xgb.train(params, dtr, num_boost_round=500)


def ece(p, y, nbins=10):
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


def run(top1, extra, groups, tag):
    X = prep(top1, extra)
    y = top1['hit_label'].values
    n_num = len(BASE_NUMERIC) + len(extra)
    n_cat = len(CATEGORICAL)
    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        m = train(X.iloc[tr], y[tr], n_num, n_cat)
        dte = xgb.DMatrix(X.iloc[te], enable_categorical=True,
                          feature_types=['q'] * n_num + ['c'] * n_cat)
        oof[te] = m.predict(dte)
        a = roc_auc_score(y[te], oof[te])
        fold_aucs.append(a)
        print(f'  [{tag}] Fold {fold}: AUC={a:.4f}')
    auc = roc_auc_score(y, oof)
    print(f'  [{tag}] OOF AUC: {auc:.4f}   per-fold: '
          f'[{", ".join(f"{a:.4f}" for a in fold_aucs)}]')

    # Final model for feature importance
    m_final = train(X, y, n_num, n_cat)
    imp = m_final.get_score(importance_type='gain')
    return oof, fold_aucs, auc, imp


def main():
    print(f'Loading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows, {ft["wiki_id"].nunique():,} unique spectra')

    missing = [c for c in NOISE_FEATURES if c not in ft.columns]
    if missing:
        raise SystemExit(f'Missing noise feature columns: {missing}. '
                          f'Run code/add_noise_features.py first.')

    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1_idx = labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    print(f'  {len(top1):,} top-1 rows; prior_TP={top1["hit_label"].mean():.3f}')

    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    # NaN audit on noise features
    for c in NOISE_FEATURES:
        n_nan = top1[c].isna().sum()
        print(f'  noise feature {c:22s} NaN in top-1: {n_nan}')

    print('\n=== Baseline GBM (no noise features) ===')
    oof_b, folds_b, auc_b, imp_b = run(top1, [], groups, 'baseline')

    print('\n=== Extended GBM (+ Oliver noise features) ===')
    oof_e, folds_e, auc_e, imp_e = run(top1, NOISE_FEATURES, groups, '+noise')

    # Calibrate + report
    rows = []
    for name, oof, folds, auc, imp in [('baseline', oof_b, folds_b, auc_b, imp_b),
                                         ('+noise', oof_e, folds_e, auc_e, imp_e)]:
        y = top1['hit_label'].values
        iso = IsotonicRegression(out_of_bounds='clip')
        iso.fit(oof, y)
        oc = iso.transform(oof)
        rows.append({
            'arm': name, 'auc_oof': auc,
            **{f'auc_fold{i}': folds[i] for i in range(5)},
            'brier_raw': brier_score_loss(y, oof),
            'brier_cal': brier_score_loss(y, oc),
            'ece_raw': ece(oof, y), 'ece_cal': ece(oc, y),
            'delta_vs_baseline': auc - auc_b,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_OUT, index=False)

    print(f'\n=== Final summary ===')
    cols_show = ['arm', 'auc_oof', 'delta_vs_baseline', 'brier_cal', 'ece_cal']
    print(summary[cols_show].to_string(index=False, float_format=lambda x: f'{x:.4f}'))
    print(f'\nWrote {SUMMARY_OUT}')

    # Feature importance comparison
    print('\n=== Top-15 feature importance (gain) — extended model ===')
    imp_e_sorted = sorted(imp_e.items(), key=lambda x: -x[1])
    for f, g in imp_e_sorted[:15]:
        flag = '  *' if f in NOISE_FEATURES else ''
        print(f'  {f:28s}  {g:>12.2f}{flag}')


if __name__ == '__main__':
    main()
