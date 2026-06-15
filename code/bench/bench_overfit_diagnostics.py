"""
Overfitting / leakage / circularity diagnostics for score_gbm_v2.py.

Reuses the EXACT production model (train_xgb params, prep_features, 5-fold
GroupKFold by anno_ik14) on the documented 15-feature base set (the 2 relational
features are excluded: they introduce label-dependent feature construction that
would muddy the label-shuffle test, and the base set reproduces the documented
pre-rel OOF AUC ~0.9131).

Three tests:
  A. Classic overfitting   — train-vs-OOF gap + per-fold AUC spread.
  B. Leakage               — group-level label-shuffle permutation (expect AUC ~0.5).
  C. Circularity           — RT ablation + adduct ablation (how much OOF AUC leans
                             on the features that proxy Oliver's labeling rule).
  D. Min sanity (pos-only) — OOF confidence on the 109 standard-verified compounds.
                             (No negatives in Min => discrimination AUC undefined;
                             this is a recall/calibration sanity check, not an AUC.)

Run from repo root:
  PYTHONPATH=code .venv_bench/bin/python code/bench_overfit_diagnostics.py
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

import score_gbm_v2 as sg
import xgboost as xgb

ROOT = sg.ROOT
MIN_PATH = os.path.join(ROOT, 'benchmark', '_internal', 'min_verified_entries.csv')

# --- 15-feature base set (production minus the 2 relational features) ---
BASE_NUMERIC = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEG = ['hit_adduct_cat', 'db', 'polarity']
ADDUCT_FEATS = ['compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious',
                'hit_isf_no_ok', 'n_candidate_adducts', 'hit_adduct_cat']


def build_training_set():
    """Replicate score_gbm_v2.main()'s training-set construction verbatim."""
    ft = pd.read_csv(sg.FEATURE_TABLE, low_memory=False)
    ft = ft[ft['hit_ik14'].fillna('').ne('')].reset_index(drop=True)
    overrides_path = os.path.join(ROOT, 'data', 'oliver_extra_labels_20260520.csv')
    ft['oliver_override_source'] = ''
    if os.path.exists(overrides_path):
        overrides = pd.read_csv(overrides_path)
        for _, ov in overrides.iterrows():
            mask = ft['wiki_id'] == ov['wiki_id']
            if not mask.any():
                continue
            top1_pos = ft[mask]['entropy_similarity'].idxmax()
            ft.loc[mask, 'spectrum_label'] = ov['spectrum_label']
            ft.loc[top1_pos, 'hit_label'] = int(ov['hit_label'])
            ft.loc[top1_pos, 'anno_ik14'] = ft.loc[top1_pos, 'hit_ik14']
            ft.loc[mask, 'oliver_override_source'] = ov['source']
    labeled = ft['spectrum_label'].isin(['TP', 'FP']) & ft['anno_ik14'].fillna('').ne('')
    top1_idx = ft[labeled].groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    return top1


def make_groups(top1):
    g = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(g)):
        if g[i] == '':
            g[i] = f'__no_ik14_{i}'
    return g


def oof_run(X, y, groups, feature_cols, return_train=False):
    """5-fold GroupKFold OOF using the production xgb params. Optionally also
    return mean in-sample (train-on-fold) AUC to measure the train-vs-OOF gap."""
    n_num = len(feature_cols) - len(CATEG)
    ftypes = ['q'] * n_num + ['c'] * len(CATEG)
    params = dict(objective='binary:logistic', eval_metric='auc', tree_method='hist',
                  max_depth=5, learning_rate=0.05, subsample=0.85, colsample_bytree=0.85,
                  min_child_weight=5, reg_alpha=0.1, reg_lambda=1.0, seed=42, verbosity=0)
    gkf = GroupKFold(n_splits=5)
    oof = np.full(len(y), np.nan)
    fold_aucs, train_aucs = [], []
    for tr, te in gkf.split(X, y, groups):
        dtr = xgb.DMatrix(X.iloc[tr][feature_cols], label=y[tr],
                          enable_categorical=True, feature_types=ftypes)
        dte = xgb.DMatrix(X.iloc[te][feature_cols],
                          enable_categorical=True, feature_types=ftypes)
        m = xgb.train(params, dtr, num_boost_round=500)
        oof[te] = m.predict(dte)
        fold_aucs.append(roc_auc_score(y[te], oof[te]))
        if return_train:
            train_aucs.append(roc_auc_score(y[tr], m.predict(dtr)))
    auc = roc_auc_score(y, oof)
    return oof, auc, fold_aucs, (np.mean(train_aucs) if return_train else None)


def group_shuffle_labels(top1, y, seed):
    """Permute labels at the anno_ik14 GROUP level (preserves group homogeneity),
    so the shuffle destroys signal without creating intra-group label mixing."""
    rng = np.random.default_rng(seed)
    grp = top1['anno_ik14'].fillna('').values
    # group -> its (majority) label
    df = pd.DataFrame({'g': grp, 'y': y})
    # treat empty ik14 rows as their own singleton groups
    df['g'] = [gi if gi != '' else f'__s{i}' for i, gi in enumerate(grp)]
    glabel = df.groupby('g')['y'].agg(lambda s: int(round(s.mean())))
    perm = rng.permutation(glabel.values)
    gmap = dict(zip(glabel.index, perm))
    return df['g'].map(gmap).astype(int).values


def main():
    # patch the shared encoders to the base feature set
    sg.NUMERIC_FEATURES = BASE_NUMERIC
    sg.CATEGORICAL_FEATURES = CATEG
    sg.ALL_FEATURES = BASE_NUMERIC + CATEG

    print('Building training set (verbatim from score_gbm_v2.main)...')
    top1 = build_training_set()
    y = top1['hit_label'].astype(int).values
    groups = make_groups(top1)
    X, _ = sg.prep_features(top1)
    feat = BASE_NUMERIC + CATEG
    print(f'  n={len(y):,} top-1 bins   prior(TP)={y.mean():.3f}   '
          f'n_groups(anno_ik14)={len(set(groups)):,}')

    # ---- A. classic overfitting: train vs OOF gap + fold spread ----
    print('\n=== A. CLASSIC OVERFITTING (train vs OOF) ===')
    oof, auc, folds, train_auc = oof_run(X, y, groups, feat, return_train=True)
    print(f'  OOF AUC          : {auc:.4f}')
    print(f'  in-sample AUC    : {train_auc:.4f}')
    print(f'  train - OOF gap  : {train_auc - auc:.4f}')
    print(f'  per-fold OOF AUC : {", ".join(f"{a:.4f}" for a in folds)}')
    print(f'  fold spread (max-min): {max(folds) - min(folds):.4f}')

    # ---- B. leakage: group-level label shuffle ----
    print('\n=== B. LEAKAGE (group-level label-shuffle permutation) ===')
    print('  expect AUC ~0.50 if there is no leakage')
    shuf_aucs = []
    for seed in (1, 2, 3):
        ys = group_shuffle_labels(top1, y, seed)
        _, a, _, _ = oof_run(X, ys, groups, feat)
        shuf_aucs.append(a)
        print(f'  shuffle seed {seed}: OOF AUC = {a:.4f}')
    print(f'  mean shuffled AUC: {np.mean(shuf_aucs):.4f}  (real={auc:.4f})')

    # ---- C. circularity: RT and adduct ablations ----
    print('\n=== C. CIRCULARITY (ablate the labeling-rule proxies) ===')
    no_rt = [f for f in feat if f != 'signed_delta_rt']
    _, auc_no_rt, _, _ = oof_run(X, y, groups, no_rt)
    no_add = [f for f in feat if f not in ADDUCT_FEATS]
    _, auc_no_add, _, _ = oof_run(X, y, groups, no_add)
    no_both = [f for f in feat if f != 'signed_delta_rt' and f not in ADDUCT_FEATS]
    _, auc_no_both, _, _ = oof_run(X, y, groups, no_both)
    print(f'  full           : {auc:.4f}')
    print(f'  drop RT        : {auc_no_rt:.4f}   (ΔAUC {auc_no_rt - auc:+.4f})')
    print(f'  drop adduct    : {auc_no_add:.4f}   (ΔAUC {auc_no_add - auc:+.4f})')
    print(f'  drop RT+adduct : {auc_no_both:.4f}   (ΔAUC {auc_no_both - auc:+.4f})')

    # ---- D. Min positive-only sanity (OOF confidence on verified-correct) ----
    print('\n=== D. MIN SANITY (109 standard-verified, positive-only) ===')
    print('  NOTE: Min has no negatives => discrimination AUC is undefined.')
    print('  Reporting OOF confidence on Min bins that are held-out in CV.')
    if os.path.exists(MIN_PATH):
        mn = pd.read_csv(MIN_PATH)
        top1 = top1.assign(_oof=oof)
        m = top1[top1['wiki_id'].isin(set(mn['wiki_id']))]
        found = len(m)
        print(f'  Min entries: {len(mn)}   matched in training top-1: {found}')
        if found:
            lab = m['hit_label'].astype(int)
            print(f'    of matched: {int((lab==1).sum())} labeled TP, {int((lab==0).sum())} labeled FP')
            s = m['_oof'].values
            print(f'    OOF confidence (raw GBM, pre-calibration): '
                  f'median={np.median(s):.3f}  mean={np.mean(s):.3f}  min={np.min(s):.3f}')
            for t in (0.3, 0.5, 0.7):
                print(f'      frac >= {t:.1f}: {(s >= t).mean():.3f}')
            low = m[m['_oof'] < 0.5][['wiki_id', 'name', 'adduct', 'hit_label', '_oof']] \
                if 'name' in m.columns else m[m['_oof'] < 0.5][['wiki_id', 'hit_label', '_oof']]
            if len(low):
                print(f'    verified-correct compounds scored <0.5 (n={len(low)}):')
                print(low.to_string(index=False, max_rows=20))
    else:
        print(f'  Min file not found at {MIN_PATH}')

    print('\nDone.')


if __name__ == '__main__':
    main()
