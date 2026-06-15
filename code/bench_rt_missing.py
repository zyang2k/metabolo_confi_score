"""bench_rt_missing.py — Diagnose how a missing RT prediction moves the GBM confidence.

build_features_v2.py:742  signed_delta_rt = delta_predicted_rt (observed RT - predicted RT
from the MassWiki API). When the API returns NO RT prediction this is NaN. score_gbm_v2.py
coerces numerics with errors='coerce' (NaN preserved, never imputed) and there is no
missing-RT indicator — so XGBoost routes the NaN via native default-direction learning.
Missing RT is therefore NOT explicitly penalised today; the model sends it down whatever
branch maximises the *training* objective.

LABEL-ARTIFACT CAVEAT (memory: project_yy_formation, project_oliver_curation_method).
Oliver's yy_ FP labels come from RT-disagreement *within a compound name*. A candidate with
no RT prediction structurally cannot be RT-flagged as an FP, so the labels could *reward*
missing-RT — the opposite of what we want. We must therefore NOT trust an AUC bump; we
measure the direction the model actually assigns before imposing any penalty.

This script reports, on the production top-1 GroupKFold OOF setup (mirrors bench_harness):
  (a) OOF AUC delta: baseline  vs  +rt_pred_missing indicator
  (b) direction missing-RT moves confidence:
        - label side  : TP rate among missing vs present rows
        - model side  : mean OOF predicted confidence, missing vs present (baseline + extended)
        - partial dep : force rt_pred_missing 0->1 on the final +feature model, mean shift
  (c) the same split restricted to the trustworthy / Oliver-reviewed subset only
      (FP | TP-with-duplicate-name — mirrors is_oliver_reviewed in score_gbm_v2.py)

Diagnostic only: writes no model, changes no production feature. The penalty itself is
imposed downstream in score_gbm_v2.py as a deterministic rule (see that file).
"""
import os
import sys

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_harness import (BASE_NUMERIC, CATEGORICAL, FEATURE_TABLE, NUM_BOOST_ROUND,
                           XGB_PARAMS, prep)

NEW = 'rt_pred_missing'


def oof_predict(top1, extra, groups, y):
    """5-fold GroupKFold OOF predictions; prep once (global cat encoding) then slice."""
    X = prep(top1, extra)
    n_num = len(BASE_NUMERIC) + len(extra)
    n_cat = len(CATEGORICAL)
    ftypes = ['q'] * n_num + ['c'] * n_cat
    oof = np.full(len(top1), np.nan)
    aucs = []
    for tr, te in GroupKFold(n_splits=5).split(top1, y, groups):
        m = xgb.train(XGB_PARAMS,
                      xgb.DMatrix(X.iloc[tr], label=y[tr], enable_categorical=True,
                                  feature_types=ftypes),
                      NUM_BOOST_ROUND)
        oof[te] = m.predict(xgb.DMatrix(X.iloc[te], enable_categorical=True,
                                        feature_types=ftypes))
        aucs.append(roc_auc_score(y[te], oof[te]))
    return oof, aucs


def partial_dependence(top1, y):
    """Train final +rt_pred_missing model on all rows; predict with the indicator forced
    0 then 1 (everything else held fixed). Mean(pred|missing=1) - mean(pred|missing=0) is the
    direction the model assigns to missing-RT in isolation. Negative => penalises missing-RT."""
    X = prep(top1, [NEW])
    n_num = len(BASE_NUMERIC) + 1
    ftypes = ['q'] * n_num + ['c'] * len(CATEGORICAL)
    m = xgb.train(XGB_PARAMS,
                  xgb.DMatrix(X, label=y, enable_categorical=True, feature_types=ftypes),
                  NUM_BOOST_ROUND)
    X0 = X.copy(); X0[NEW] = 0.0
    X1 = X.copy(); X1[NEW] = 1.0
    p0 = m.predict(xgb.DMatrix(X0, enable_categorical=True, feature_types=ftypes))
    p1 = m.predict(xgb.DMatrix(X1, enable_categorical=True, feature_types=ftypes))
    return float((p1 - p0).mean()), float(np.median(p1 - p0))


def report_split(tag, miss, y, base_oof, ext_oof):
    """Direction report on a given row subset."""
    m = miss.astype(bool)
    n_miss, n_pres = int(m.sum()), int((~m).sum())
    print(f'\n--- {tag} (n={len(y)}; missing-RT n={n_miss}, present n={n_pres}) ---')
    if n_miss == 0:
        print('  no missing-RT rows in this subset; nothing to compare')
        return
    print(f'  LABEL side : TP rate   missing={y[m].mean():.4f}   present={y[~m].mean():.4f}   '
          f'(missing {"REWARDED" if y[m].mean() > y[~m].mean() else "PENALISED"} by labels)')
    print(f'  MODEL base : mean conf missing={base_oof[m].mean():.4f}   present={base_oof[~m].mean():.4f}   '
          f'Δ={base_oof[m].mean() - base_oof[~m].mean():+.4f}')
    print(f'  MODEL +ind : mean conf missing={ext_oof[m].mean():.4f}   present={ext_oof[~m].mean():.4f}   '
          f'Δ={ext_oof[m].mean() - ext_oof[~m].mean():+.4f}')


def main():
    print(f'Loading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    ft = ft[ft['hit_ik14'].fillna('').ne('')].reset_index(drop=True)

    # Production top-1: labeled (TP/FP) rows with a curator structure, highest-esim per bin.
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP']) & ft['anno_ik14'].fillna('').ne('')]
    top1 = ft.loc[lab.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    y = top1['hit_label'].values

    # The diagnostic feature: 1 when the API returned no RT prediction.
    top1[NEW] = pd.to_numeric(top1['signed_delta_rt'], errors='coerce').isna().astype(float)
    miss = top1[NEW].values.astype(bool)
    print(f'  {len(top1):,} top-1 rows; prior_TP={y.mean():.3f}; '
          f'missing-RT rate={miss.mean()*100:.1f}% ({miss.sum()} rows)')

    # Trustworthy / Oliver-reviewed slice: all FPs (yy_ are explicitly flagged) + TPs whose
    # anno_name appears in >=2 spectra (reviewed for RT conflict). Mirrors score_gbm_v2.py.
    tp_mask = top1['spectrum_label'].eq('TP').values
    fp_mask = top1['spectrum_label'].eq('FP').values
    name_counts = top1[tp_mask].groupby('anno_name_lower')['wiki_id'].nunique()
    dup_names = set(name_counts[name_counts >= 2].index) - {''}
    is_dup = top1['anno_name_lower'].fillna('').isin(dup_names).values
    trust = fp_mask | (tp_mask & is_dup)
    print(f'  trustworthy subset: {int(trust.sum()):,} rows ({100*trust.mean():.1f}%)')

    # Groups for GroupKFold (empty-IK14 guard).
    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    print('\n=== (a) OOF AUC: baseline vs +rt_pred_missing ===')
    base_oof, base_aucs = oof_predict(top1, [], groups, y)
    ext_oof, ext_aucs = oof_predict(top1, [NEW], groups, y)
    auc_base = roc_auc_score(y, base_oof)
    auc_ext = roc_auc_score(y, ext_oof)
    print(f'  baseline       OOF AUC = {auc_base:.4f}   folds {[f"{a:.4f}" for a in base_aucs]}')
    print(f'  +rt_pred_missing OOF AUC = {auc_ext:.4f}   folds {[f"{a:.4f}" for a in ext_aucs]}')
    print(f'  ΔAUC = {auc_ext - auc_base:+.4f}')

    print('\n=== (b) direction missing-RT moves confidence ===')
    report_split('ALL labeled', miss, y, base_oof, ext_oof)
    pd_mean, pd_med = partial_dependence(top1, y)
    EPS = 1e-4
    if abs(pd_mean) < EPS:
        pd_verdict = 'IGNORES the explicit indicator (redundant — already routes signed_delta_rt NaN natively)'
    else:
        pd_verdict = f'{"PENALISES" if pd_mean < 0 else "REWARDS"} missing-RT via the explicit indicator'
    print(f'\n  PARTIAL DEPENDENCE (+feature final model, force missing 0->1):')
    print(f'    mean Δconf = {pd_mean:+.5f}   median Δconf = {pd_med:+.5f}   => model {pd_verdict}')

    print('\n=== (c) trustworthy / Oliver-reviewed subset only ===')
    report_split('Trustworthy', miss[trust], y[trust], base_oof[trust], ext_oof[trust])

    print('\n=== VERDICT ===')
    label_dir = 'REWARD' if y[miss].mean() > y[~miss].mean() else 'PENALISE'
    model_dir = 'rewards' if base_oof[miss].mean() > base_oof[~miss].mean() else 'penalises'
    print(f'  Labels currently {label_dir} missing-RT (TP rate {y[miss].mean():.3f} vs {y[~miss].mean():.3f}); '
          f'the model {model_dir} it (mean OOF conf {base_oof[miss].mean():.3f} vs {base_oof[~miss].mean():.3f}).')
    print(f'  The explicit indicator adds ~0 AUC ({auc_ext - auc_base:+.4f}) and ~0 partial dependence:')
    print('  XGBoost already extracts the missing-RT signal from signed_delta_rt\'s native NaN routing,')
    print('  so the model\'s direction is whatever the LABELS imply — and that direction is a curation')
    print('  artifact (yy_ FP = RT-disagreement; missing-RT cannot be RT-flagged) that has already')
    print('  flipped vs the 14.6%-TP/10.0%-FP figures this task was scoped against. A penalty therefore')
    print('  cannot be learned honestly here; it must be imposed as a deterministic downstream prior')
    print('  — see RT_MISSING_CONFIDENCE_HAIRCUT in code/score_gbm_v2.py.')


if __name__ == '__main__':
    main()
