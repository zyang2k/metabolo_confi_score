"""bench_harness.py — Shared GBM feature-ablation harness.

Extracted from the bench_*.py scripts that all repeat the same shape:
load feature_table_v2.csv -> top-1 per wiki_id -> 5-fold GroupKFold on
anno_ik14 (with empty-IK14 guard) -> baseline vs +feature arms ->
isotonic-calibrated AUC / Brier / ECE summary + feature importance.

Single source of truth for:
  * BASE_NUMERIC / CATEGORICAL feature lists (mirror score_gbm_v2.py)
  * XGB hyperparameters
  * Group construction with empty-IK14 guard
  * ECE definition

A new bench script becomes ~30 lines: import this module, declare the
extra feature column(s) and (optional) side-file merge, call run_bench().
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')

# Mirrors score_gbm_v2.py production feature set.
BASE_NUMERIC = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEGORICAL = ['hit_adduct_cat', 'db', 'polarity']

XGB_PARAMS = {
    'objective': 'binary:logistic', 'eval_metric': 'auc',
    'tree_method': 'hist', 'max_depth': 5, 'learning_rate': 0.05,
    'subsample': 0.85, 'colsample_bytree': 0.85, 'min_child_weight': 5,
    'reg_alpha': 0.1, 'reg_lambda': 1.0, 'seed': 42, 'verbosity': 0,
}
NUM_BOOST_ROUND = 500


@dataclass
class ArmResult:
    name: str
    oof: np.ndarray
    fold_aucs: list
    auc: float
    importance: dict


@dataclass
class BenchSpec:
    """Declarative bench config — one per `python code/bench_<name>.py`."""
    name: str
    extra_numeric: list = field(default_factory=list)
    extra_table: str | None = None          # side CSV with extra feature(s)
    merge_on: tuple = ('wiki_id', 'library_wiki_id')
    pre_hook: Callable | None = None        # df -> df, applied before top-1
    summary_out: str | None = None          # default: data/bench_<name>_summary.csv


def prep(df: pd.DataFrame, extra_num: list) -> pd.DataFrame:
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


def train(X, y, n_num, n_cat, params=None):
    p = dict(XGB_PARAMS)
    if params:
        p.update(params)
    dtr = xgb.DMatrix(X, label=y, enable_categorical=True,
                      feature_types=['q'] * n_num + ['c'] * n_cat)
    return xgb.train(p, dtr, num_boost_round=NUM_BOOST_ROUND)


def ece(p, y, nbins=10) -> float:
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


def build_top1(ft: pd.DataFrame) -> pd.DataFrame:
    """Labeled subset -> highest-entropy_similarity row per wiki_id."""
    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1_idx = labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    return ft.loc[top1_idx].reset_index(drop=True)


def make_groups(top1: pd.DataFrame) -> np.ndarray:
    """IK14 GroupKFold groups with empty-IK14 guard (each empty -> own group)."""
    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'
    return groups


def run_arm(top1: pd.DataFrame, extra: list, groups: np.ndarray, tag: str,
            params: dict | None = None) -> ArmResult:
    X = prep(top1, extra)
    y = top1['hit_label'].values
    n_num = len(BASE_NUMERIC) + len(extra)
    n_cat = len(CATEGORICAL)
    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        m = train(X.iloc[tr], y[tr], n_num, n_cat, params)
        dte = xgb.DMatrix(X.iloc[te], enable_categorical=True,
                          feature_types=['q'] * n_num + ['c'] * n_cat)
        oof[te] = m.predict(dte)
        a = roc_auc_score(y[te], oof[te])
        fold_aucs.append(a)
        print(f'  [{tag}] Fold {fold}: AUC={a:.4f}')
    auc = roc_auc_score(y, oof)
    print(f'  [{tag}] OOF AUC: {auc:.4f}   per-fold: '
          f'[{", ".join(f"{a:.4f}" for a in fold_aucs)}]')
    m_final = train(X, y, n_num, n_cat, params)
    imp = m_final.get_score(importance_type='gain')
    return ArmResult(tag, oof, fold_aucs, auc, imp)


def summarize(arms: list[ArmResult], y: np.ndarray, out_path: str) -> pd.DataFrame:
    base_auc = arms[0].auc
    rows = []
    for a in arms:
        iso = IsotonicRegression(out_of_bounds='clip')
        iso.fit(a.oof, y)
        oc = iso.transform(a.oof)
        rows.append({
            'arm': a.name, 'auc_oof': a.auc,
            **{f'auc_fold{i}': a.fold_aucs[i] for i in range(5)},
            'brier_raw': brier_score_loss(y, a.oof),
            'brier_cal': brier_score_loss(y, oc),
            'ece_raw': ece(a.oof, y), 'ece_cal': ece(oc, y),
            'delta_vs_baseline': a.auc - base_auc,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out_path, index=False)
    print(f'\n=== Final summary ===')
    cols_show = ['arm', 'auc_oof', 'delta_vs_baseline', 'brier_cal', 'ece_cal']
    print(summary[cols_show].to_string(index=False,
                                        float_format=lambda x: f'{x:.4f}'))
    print(f'\nWrote {out_path}')
    return summary


def print_importance(arms: list[ArmResult], extra: list, topk: int = 15):
    last = arms[-1]
    print(f'\n=== Top-{topk} feature importance (gain) — {last.name} ===')
    for f, g in sorted(last.importance.items(), key=lambda x: -x[1])[:topk]:
        flag = '  *' if f in extra else ''
        print(f'  {f:28s}  {g:>12.2f}{flag}')


def run_bench(spec: BenchSpec, feature_table_path: str = FEATURE_TABLE):
    """End-to-end: load, merge optional side table, top-1, baseline + extended."""
    print(f'Loading {feature_table_path}')
    ft = pd.read_csv(feature_table_path, low_memory=False)
    print(f'  {len(ft):,} rows, {ft["wiki_id"].nunique():,} unique spectra')

    if spec.extra_table:
        print(f'Loading side table {spec.extra_table}')
        side = pd.read_csv(spec.extra_table)
        keep = list(spec.merge_on) + spec.extra_numeric
        side = side[[c for c in keep if c in side.columns]]
        ft = ft.merge(side, on=list(spec.merge_on), how='left')
        for c in spec.extra_numeric:
            cov = ft[c].notna().mean() * 100
            print(f'  merged {c}: {cov:.1f}% coverage')

    if spec.pre_hook:
        ft = spec.pre_hook(ft)

    missing = [c for c in spec.extra_numeric if c not in ft.columns]
    if missing:
        raise SystemExit(f'Missing extra feature columns: {missing}')

    top1 = build_top1(ft)
    print(f'  {len(top1):,} top-1 rows; prior_TP={top1["hit_label"].mean():.3f}')

    for c in spec.extra_numeric:
        n_nan = top1[c].isna().sum()
        print(f'  extra feature {c:28s} NaN in top-1: {n_nan}')

    groups = make_groups(top1)

    print('\n=== Baseline GBM ===')
    base = run_arm(top1, [], groups, 'baseline')

    print(f'\n=== Extended GBM (+ {", ".join(spec.extra_numeric)}) ===')
    ext = run_arm(top1, spec.extra_numeric, groups, f'+{spec.name}')

    y = top1['hit_label'].values
    out = spec.summary_out or os.path.join(
        ROOT, 'data', f'bench_{spec.name}_summary.csv')
    summarize([base, ext], y, out)
    print_importance([base, ext], spec.extra_numeric)

    return {'baseline': base, 'extended': ext, 'top1': top1, 'summary_path': out}
