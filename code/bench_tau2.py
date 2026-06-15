"""
bench_tau2.py — Ablation: does the channel-disagreement τ² feature help GBM?

Motivation: Huang (Metrology 2026) shows that the LCD / multiplicative pooling
that the Bayesian scorer implements (sum of per-channel logLRs ≡ normalized
product of channel likelihoods) discards heterogeneity information — the
between-channel disagreement variance τ². The weighted-linear-pooling form
preserves both the within-source variance σ²_WA and τ². If our 3-channel
Bayesian leaves channel-disagreement signal on the table, surfacing τ² as an
explicit feature should let the GBM pick up on it.

Per-row τ² is computed exactly as in Huang Eq. 8:
    τ²_row = Σ_k w_k · (μ_k − μ_L)²    where w_k = (1/σ²_k) / Σ(1/σ²_i),
                                       μ_L = Σ w_k μ_k

Concretely:
    μ_k = per-row, per-channel logLR (the channel's TP-vs-FP "estimate")
    σ²_k = empirical variance of channel-k logLR across the training fold
    Channels with non-finite raw input on a row are dropped from that row's τ².

To avoid leakage, σ²_k and the Bayesian channel fits are recomputed within
each GroupKFold fold; τ² for held-out rows uses fold-fit parameters.

Reports baseline GBM (the production feature set) vs baseline + τ², comparing
5-fold OOF AUC, Brier, ECE; also dumps feature importance for τ².

Usage:
    python code/bench_tau2.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct

FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')

NUMERIC_FEATURES = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CATEGORICAL_FEATURES = ['hit_adduct_cat', 'db', 'polarity']

# Channels for τ²: same as production Bayesian scorer (entropy_sim, sim_gap, RT)
TAU2_CHANNELS = [
    ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal',
                transform=logit_upper_half),
    ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),
    ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                higher_means_tp=True, tp_family='student_t', fp_family='student_t'),
]

RT_UNRELIABLE_NORMALIZED_ADDUCTS = {'Cat', 'Anion'}


def mask_rt_for_unreliable_adducts(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    norm = out['adduct'].apply(norm_adduct)
    mask = norm.isin(RT_UNRELIABLE_NORMALIZED_ADDUCTS)
    out.loc[mask, 'signed_delta_rt'] = np.nan
    return out


def per_channel_logLR_table(df: pd.DataFrame, fitted_channels) -> np.ndarray:
    """Return (n_rows, n_channels) array of per-row, per-channel logLR.
    Sets logLR to NaN where the underlying feature value is non-finite, so
    downstream τ² computation can drop missing channels per-row."""
    n = len(df)
    K = len(fitted_channels)
    out = np.full((n, K), np.nan, dtype=float)
    for k, fc in enumerate(fitted_channels):
        raw = df[fc.spec.feature_col].values.astype(float)
        valid = np.isfinite(raw)
        if valid.sum() == 0:
            continue
        lr = fc.logLR(raw)  # already 0 at non-finite, but we'll mask explicitly
        out[valid, k] = lr[valid]
    return out


def compute_tau2(lr_matrix: np.ndarray, sigma2: np.ndarray) -> np.ndarray:
    """Per-row τ² = inverse-variance weighted variance of available channel logLRs.

    lr_matrix: (n, K), NaN where channel is missing for that row.
    sigma2:    (K,), per-channel variance of logLR (estimated on training fold).

    For each row, channels with NaN are dropped; remaining channels get weights
    w_k = (1/σ²_k) / Σ(1/σ²_i), μ_L = Σ w_k μ_k, τ² = Σ w_k (μ_k − μ_L)².
    Rows with <2 valid channels return τ²=0 (no disagreement defined).
    """
    n, K = lr_matrix.shape
    inv_sigma2 = np.where(sigma2 > 0, 1.0 / sigma2, 0.0)  # (K,)

    valid_mask = np.isfinite(lr_matrix)                  # (n, K)
    n_valid = valid_mask.sum(axis=1)                     # (n,)

    # Substitute zero for NaN so contributions are zeroed where masked
    lr_safe = np.where(valid_mask, lr_matrix, 0.0)
    w_unnorm = valid_mask * inv_sigma2[None, :]          # (n, K), 0 where invalid
    W = w_unnorm.sum(axis=1)                             # (n,)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu_L = (w_unnorm * lr_safe).sum(axis=1) / np.where(W > 0, W, 1.0)
        diff2 = (lr_safe - mu_L[:, None]) ** 2
        tau2 = (w_unnorm * diff2 * valid_mask).sum(axis=1) / np.where(W > 0, W, 1.0)

    tau2[n_valid < 2] = 0.0
    tau2[~np.isfinite(tau2)] = 0.0
    return tau2


def prep_features(df: pd.DataFrame, feature_list, cat_maps=None):
    """Coerce numeric + label-encode categoricals."""
    X = df.copy()
    numeric = [f for f in feature_list if f not in CATEGORICAL_FEATURES]
    categorical = [f for f in feature_list if f in CATEGORICAL_FEATURES]

    for c in numeric:
        if c not in X.columns:
            X[c] = np.nan
        X[c] = pd.to_numeric(X[c], errors='coerce')

    if cat_maps is None:
        cat_maps = {}
        for c in categorical:
            if c not in X.columns:
                X[c] = 'missing'
            s = X[c].astype(str).fillna('missing')
            categories = sorted(s.unique().tolist())
            mapping = {v: i for i, v in enumerate(categories)}
            X[c] = s.map(mapping).astype('int32')
            cat_maps[c] = mapping
    else:
        for c in categorical:
            s = X[c].astype(str).fillna('missing')
            mapping = cat_maps[c]
            max_code = max(mapping.values()) + 1 if mapping else 0
            X[c] = s.map(lambda v: mapping.get(v, max_code)).astype('int32')

    return X[numeric + categorical], cat_maps


def train_xgb(X, y, feature_list):
    numeric = [f for f in feature_list if f not in CATEGORICAL_FEATURES]
    categorical = [f for f in feature_list if f in CATEGORICAL_FEATURES]
    types = ['q'] * len(numeric) + ['c'] * len(categorical)
    dtrain = xgb.DMatrix(X[numeric + categorical], label=y,
                         enable_categorical=True, feature_types=types)
    params = {
        'objective': 'binary:logistic', 'eval_metric': 'auc',
        'tree_method': 'hist', 'max_depth': 5, 'learning_rate': 0.05,
        'subsample': 0.85, 'colsample_bytree': 0.85, 'min_child_weight': 5,
        'reg_alpha': 0.1, 'reg_lambda': 1.0, 'seed': 42, 'verbosity': 0,
    }
    return xgb.train(params, dtrain, num_boost_round=500), types


def predict_xgb(model, X, feature_list, types):
    numeric = [f for f in feature_list if f not in CATEGORICAL_FEATURES]
    categorical = [f for f in feature_list if f in CATEGORICAL_FEATURES]
    dmat = xgb.DMatrix(X[numeric + categorical],
                       enable_categorical=True, feature_types=types)
    return model.predict(dmat)


def ece_score(p, y, nbins=10):
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0: continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / n


def run_oof(top1_train, labels, feature_list, tau2_oof=None):
    """Return (oof_predictions, fold_aucs, final_model_importance).

    If feature_list contains 'tau2', tau2_oof must be provided (n,) array
    aligned with top1_train rows."""
    df = top1_train.copy()
    if 'tau2' in feature_list:
        df['tau2'] = tau2_oof

    X_full, cat_maps = prep_features(df, feature_list)

    groups = top1_train['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    gkf = GroupKFold(n_splits=5)
    oof = np.full(len(top1_train), np.nan)
    fold_aucs = []
    final_imp = None
    for fold, (tr, te) in enumerate(gkf.split(top1_train, labels, groups)):
        model, types = train_xgb(X_full.iloc[tr], labels[tr], feature_list)
        oof[te] = predict_xgb(model, X_full.iloc[te], feature_list, types)
        fold_aucs.append(roc_auc_score(labels[te], oof[te]))

    # Train one final model on all training data for importance reporting
    model_full, _ = train_xgb(X_full, labels, feature_list)
    final_imp = model_full.get_score(importance_type='gain')
    return oof, fold_aucs, final_imp


def compute_tau2_oof(top1_train, ft, labels, channels):
    """Per-fold Bayesian fit → per-fold τ² for held-out rows.
    Returns:
      tau2_top1_oof: (n_top1,) τ² for training rows, computed under fold's
                     Bayesian fit on the OTHER folds (no leakage).
      tau2_full:     (n_ft,)   τ² for all rows in `ft`, computed by fitting
                     Bayesian on ALL labeled top-1 rows.
    """
    groups = top1_train['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    tau2_top1_oof = np.full(len(top1_train), np.nan)
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1_train, labels, groups)):
        fitted = [fit_channel(s, top1_train.iloc[tr][s.feature_col].values, labels[tr])
                  for s in channels]
        # Compute σ²_k from the channel logLR distribution on the training fold
        lr_train = per_channel_logLR_table(top1_train.iloc[tr], fitted)
        sigma2 = np.array([np.nanvar(lr_train[:, k]) if np.isfinite(lr_train[:, k]).any() else 0.0
                           for k in range(len(channels))])
        sigma2 = np.where(sigma2 > 0, sigma2, 1e-6)
        # Apply to held-out rows
        lr_te = per_channel_logLR_table(top1_train.iloc[te], fitted)
        tau2_top1_oof[te] = compute_tau2(lr_te, sigma2)

    # For non-training rows, fit on ALL labeled top-1 rows once
    fitted_full = [fit_channel(s, top1_train[s.feature_col].values, labels)
                   for s in channels]
    lr_full = per_channel_logLR_table(ft, fitted_full)
    lr_train_full = per_channel_logLR_table(top1_train, fitted_full)
    sigma2_full = np.array([np.nanvar(lr_train_full[:, k]) if np.isfinite(lr_train_full[:, k]).any() else 0.0
                            for k in range(len(channels))])
    sigma2_full = np.where(sigma2_full > 0, sigma2_full, 1e-6)
    tau2_full = compute_tau2(lr_full, sigma2_full)

    return tau2_top1_oof, tau2_full, sigma2_full


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows × {len(ft.columns)} cols, {ft["wiki_id"].nunique():,} spectra')

    print('\nApplying RT mask for [Cat]+/[Anion]- adducts...')
    ft = mask_rt_for_unreliable_adducts(ft)

    labeled_mask = ft['spectrum_label'].isin(['TP', 'FP'])
    ft_labeled = ft[labeled_mask]
    top1_idx = ft_labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1_train = ft.loc[top1_idx].reset_index(drop=True)
    labels = top1_train['hit_label'].values
    prior_tp = labels.mean()
    print(f'\nTraining set: {len(top1_train):,} top-1 rows (TP+FP)   prior={prior_tp:.3f}')

    # ── τ² computation ──
    print('\nFitting per-fold Bayesian channels and computing τ²...')
    tau2_top1_oof, tau2_full, sigma2_full = compute_tau2_oof(
        top1_train, ft, labels, TAU2_CHANNELS)
    print(f'  Per-channel σ² (full-fit): '
          f'entropy_sim={sigma2_full[0]:.3f}  sim_gap={sigma2_full[1]:.3f}  '
          f'signed_delta_rt={sigma2_full[2]:.3f}')
    print(f'  τ² (training rows): mean={tau2_top1_oof.mean():.3f}  '
          f'median={np.median(tau2_top1_oof):.3f}  '
          f'95p={np.percentile(tau2_top1_oof, 95):.3f}  '
          f'frac_zero={(tau2_top1_oof == 0).mean():.3f}')

    # Sanity: TP vs FP τ² distribution
    tp_tau = tau2_top1_oof[labels == 1]
    fp_tau = tau2_top1_oof[labels == 0]
    print(f'  TP rows: τ² mean={tp_tau.mean():.3f}  median={np.median(tp_tau):.3f}')
    print(f'  FP rows: τ² mean={fp_tau.mean():.3f}  median={np.median(fp_tau):.3f}')
    auc_tau_alone = roc_auc_score(labels, -tau2_top1_oof)  # higher τ² → less TP-like
    print(f'  τ² alone as inverse predictor: AUC={auc_tau_alone:.4f}')

    # ── Baseline GBM ──
    print('\n=== Baseline GBM (production feature set) ===')
    feat_baseline = NUMERIC_FEATURES + CATEGORICAL_FEATURES
    oof_b, folds_b, imp_b = run_oof(top1_train, labels, feat_baseline)
    valid = ~np.isnan(oof_b)
    auc_b = roc_auc_score(labels[valid], oof_b[valid])
    brier_b = brier_score_loss(labels[valid], oof_b[valid])
    iso_b = IsotonicRegression(out_of_bounds='clip').fit(oof_b[valid], labels[valid])
    cal_b = iso_b.transform(oof_b[valid])
    ece_b = ece_score(cal_b, labels[valid])
    print(f'  OOF AUC: {auc_b:.4f}   folds: [{min(folds_b):.4f}, {max(folds_b):.4f}]')
    print(f'  OOF Brier: {brier_b:.4f}   ECE (calibrated): {ece_b:.4f}')

    # ── GBM + τ² ──
    print('\n=== GBM + τ² ===')
    feat_with_tau = feat_baseline + ['tau2']
    oof_t, folds_t, imp_t = run_oof(top1_train, labels, feat_with_tau,
                                    tau2_oof=tau2_top1_oof)
    valid = ~np.isnan(oof_t)
    auc_t = roc_auc_score(labels[valid], oof_t[valid])
    brier_t = brier_score_loss(labels[valid], oof_t[valid])
    iso_t = IsotonicRegression(out_of_bounds='clip').fit(oof_t[valid], labels[valid])
    cal_t = iso_t.transform(oof_t[valid])
    ece_t = ece_score(cal_t, labels[valid])
    print(f'  OOF AUC: {auc_t:.4f}   folds: [{min(folds_t):.4f}, {max(folds_t):.4f}]')
    print(f'  OOF Brier: {brier_t:.4f}   ECE (calibrated): {ece_t:.4f}')

    # ── Comparison ──
    print('\n=== Δ (τ² − baseline) ===')
    print(f'  ΔAUC:   {auc_t - auc_b:+.4f}')
    print(f'  ΔBrier: {brier_t - brier_b:+.4f}  (lower is better)')
    print(f'  ΔECE:   {ece_t - ece_b:+.4f}  (lower is better)')

    # ── Feature importance for τ² ──
    print('\n=== Feature importance (gain) — model with τ² ===')
    imp_sorted = sorted(imp_t.items(), key=lambda x: -x[1])
    for feat, g in imp_sorted[:20]:
        marker = '  ←' if feat == 'tau2' else ''
        print(f'  {feat:30s}  {g:>12.2f}{marker}')

    # tau2 rank
    tau2_rank = next((i for i, (f, _) in enumerate(imp_sorted) if f == 'tau2'), None)
    if tau2_rank is not None:
        print(f'\n  τ² rank: {tau2_rank + 1} of {len(imp_sorted)} features')
    else:
        print('\n  τ² did not split — model deemed it uninformative')


if __name__ == '__main__':
    main()
