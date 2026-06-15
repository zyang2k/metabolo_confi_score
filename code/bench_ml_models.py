"""
bench_ml_models.py — Compare RF, Logistic Regression, MLP against GBM on the
same feature table, folds, and calibration pipeline.

Loads the same training set as score_gbm_v2 (top-1 rows, TP+FP only), runs
5-fold GroupKFold OOF, applies isotonic, and reports AUC / Brier / ECE + deck
cases + Min validation per model. Uses the GBM OOF AUC 0.913 as the reference
point.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
MIN_VERIFIED  = os.path.join(ROOT, 'benchmark', '_internal', 'min_verified_entries.csv')

NUMERIC_FEATURES = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEGORICAL_FEATURES = ['hit_adduct_cat', 'db', 'polarity']
ALL_FEATURES = NUMERIC_FEATURES + CATEGORICAL_FEATURES


def encode(df, cat_maps=None):
    """Numeric coerce + label-encode categoricals using a stable mapping."""
    X = df[ALL_FEATURES].copy()
    for c in NUMERIC_FEATURES:
        X[c] = pd.to_numeric(X[c], errors='coerce')
    if cat_maps is None:
        cat_maps = {}
        for c in CATEGORICAL_FEATURES:
            s = X[c].astype(str).fillna('missing')
            cats = sorted(s.unique().tolist())
            cat_maps[c] = {v: i for i, v in enumerate(cats)}
    for c in CATEGORICAL_FEATURES:
        s = X[c].astype(str).fillna('missing')
        max_code = max(cat_maps[c].values()) + 1 if cat_maps[c] else 0
        X[c] = s.map(lambda v: cat_maps[c].get(v, max_code)).astype('int32')
    return X, cat_maps


def bayesian_ece(p, y, nbins=10):
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0: continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / n


MODELS = {
    'GBM (XGBoost)': 'xgb',
    'Random Forest': 'rf',
    'Logistic Regression': 'lr',
    'MLP (2×64)': 'mlp',
}


def train_and_score(model_key, X_train, y_train, X_test):
    """Fit each model and return P(y=1) for X_test."""
    if model_key == 'xgb':
        ft_types = ['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES)
        dtrain = xgb.DMatrix(X_train, label=y_train, enable_categorical=True, feature_types=ft_types)
        dtest = xgb.DMatrix(X_test, enable_categorical=True, feature_types=ft_types)
        params = dict(objective='binary:logistic', eval_metric='auc', tree_method='hist',
                      max_depth=5, learning_rate=0.05, subsample=0.85, colsample_bytree=0.85,
                      min_child_weight=5, reg_alpha=0.1, reg_lambda=1.0, seed=42, verbosity=0)
        model = xgb.train(params, dtrain, num_boost_round=500)
        return model.predict(dtest)
    if model_key == 'rf':
        # RF doesn't need scaling; uses SimpleImputer on NaN
        pipe = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('rf', RandomForestClassifier(
                n_estimators=500, max_depth=None, min_samples_leaf=5,
                max_features='sqrt', n_jobs=-1, random_state=42)),
        ])
        pipe.fit(X_train, y_train)
        return pipe.predict_proba(X_test)[:, 1]
    if model_key == 'lr':
        # LR needs scaling; treat categorical encodings as ordinal (imperfect but consistent with others)
        pipe = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('sc', StandardScaler()),
            ('lr', LogisticRegression(max_iter=2000, C=1.0, random_state=42)),
        ])
        pipe.fit(X_train, y_train)
        return pipe.predict_proba(X_test)[:, 1]
    if model_key == 'mlp':
        pipe = Pipeline([
            ('imp', SimpleImputer(strategy='median')),
            ('sc', StandardScaler()),
            ('mlp', MLPClassifier(hidden_layer_sizes=(64, 64), max_iter=500,
                                  early_stopping=True, random_state=42)),
        ])
        pipe.fit(X_train, y_train)
        return pipe.predict_proba(X_test)[:, 1]
    raise ValueError(model_key)


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1 = ft.loc[labeled.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    labels = top1['hit_label'].values
    print(f'  {len(top1):,} top-1 training rows (TP+FP)   prior={labels.mean():.3f}')

    X_train, cat_maps = encode(top1)

    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    results = {}
    for name, key in MODELS.items():
        print(f'\n=== {name} ===')
        oof = np.full(len(top1), np.nan)
        fold_aucs = []
        for fold, (tr, te) in enumerate(GroupKFold(n_splits=5).split(top1, labels, groups)):
            preds = train_and_score(key, X_train.iloc[tr], labels[tr], X_train.iloc[te])
            oof[te] = preds
            auc = roc_auc_score(labels[te], preds)
            fold_aucs.append(auc)
            print(f'  Fold {fold}: AUC={auc:.4f}')

        valid = ~np.isnan(oof)
        auc_raw = roc_auc_score(labels[valid], oof[valid])
        brier_raw = brier_score_loss(labels[valid], oof[valid])
        ece_raw = bayesian_ece(oof[valid], labels[valid])

        # Isotonic calibration
        iso = IsotonicRegression(out_of_bounds='clip')
        iso.fit(oof[valid], labels[valid])
        oof_cal = iso.transform(oof[valid])
        brier_cal = brier_score_loss(labels[valid], oof_cal)
        ece_cal = bayesian_ece(oof_cal, labels[valid])

        results[name] = dict(
            auc=auc_raw, brier_raw=brier_raw, brier_cal=brier_cal,
            ece_raw=ece_raw, ece_cal=ece_cal,
            fold_aucs=fold_aucs, oof=oof, oof_cal=oof_cal,
        )
        print(f'  AUC: {auc_raw:.4f}  (fold spread {min(fold_aucs):.4f}–{max(fold_aucs):.4f})')
        print(f'  Brier raw → cal:  {brier_raw:.4f} → {brier_cal:.4f}')
        print(f'  ECE   raw → cal:  {ece_raw:.4f} → {ece_cal:.4f}')

    # --- Summary table ---
    print('\n' + '=' * 80)
    print('SUMMARY (OOF metrics, TP+FP labeled spectra)')
    print('=' * 80)
    print(f'  {"Model":<24s}  {"AUC":>7s}  {"Brier":>7s}  {"ECE (cal)":>10s}  {"fold spread":>14s}')
    for name, r in results.items():
        print(f'  {name:<24s}  {r["auc"]:>7.4f}  {r["brier_cal"]:>7.4f}  {r["ece_cal"]:>10.4f}  '
              f'[{min(r["fold_aucs"]):.3f}, {max(r["fold_aucs"]):.3f}]')

    # --- Deck-case + Min checks for the top model ---
    best = max(results.items(), key=lambda kv: kv[1]['auc'])
    print(f'\nTop model: {best[0]} (AUC {best[1]["auc"]:.4f})')

    # Score each bin's row using per-fold OOF predictions (for labeled bins only);
    # use those for a quick deck-case check without re-fitting on full data here.
    for name, r in results.items():
        top1_cal = top1.copy()
        top1_cal['cal_conf'] = r['oof_cal']
        dc = []
        for wid in ['aPUDE1U/2244', 'aEKJ9AS/10123', 'aEKJ9AS/994', 'aEKJ9AS/1674']:
            row = top1_cal[top1_cal['wiki_id'] == wid]
            if len(row):
                dc.append(f'{wid}: {row["cal_conf"].iloc[0]*100:.0f}%')
            else:
                dc.append(f'{wid}: n/a')
        print(f'  {name:<24s}  deck: {"  |  ".join(dc)}')

    try:
        min_df = pd.read_csv(MIN_VERIFIED)[['wiki_id', 'host_ik14']]
        print(f'\nMin validation (top-1 hit IK14 match to Min truth):')
        for name, r in results.items():
            top1_cal = top1.copy()
            top1_cal['cal_conf'] = r['oof_cal']
            top1_cal['hit_ik14'] = top1_cal['hit_ik14'].fillna('')
            m = min_df.merge(top1_cal[['wiki_id', 'hit_ik14', 'cal_conf']], on='wiki_id', how='inner')
            m['match'] = (m['hit_ik14'] == m['host_ik14']).astype(int)
            mc = m[m['match'] == 1]
            high = (mc['cal_conf'] >= 0.9).sum() if len(mc) else 0
            print(f'  {name:<24s}  match {m["match"].sum()}/{len(m)}  '
                  f'median conf on correct {mc["cal_conf"].median()*100:.1f}%  '
                  f'≥90% on correct {high}/{len(mc)}')
    except FileNotFoundError:
        pass


if __name__ == '__main__':
    main()
