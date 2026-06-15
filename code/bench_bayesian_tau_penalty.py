"""
bench_bayesian_tau_penalty.py — Does penalizing logit-posterior by τ² help (b)?

Following Huang's diagnosis (Metrology 2026) that the multiplicative Bayesian
combination (≡ LCD / conflation) discards the between-channel disagreement
variance τ², we test a corrective term on the (b) track:

    logit(P) ← prior_log_odds + Σ_k logLR_k(x_k) − α · τ²(x)

where τ²(x) is the inverse-variance-weighted between-channel variance of
per-row, per-channel logLRs. Higher τ² ⇒ channels disagree ⇒ down-weight the
posterior.

α is tuned on a per-fold validation grid; we report:
  - baseline Bayesian OOF AUC (no penalty)
  - τ²-penalized OOF AUC at best α
  - the α that wins
  - calibrated Brier and ECE for both

Per-fold no-leakage:
  Within each GroupKFold fold, the Bayesian channel densities, σ²_k for τ²
  weighting, AND α are fit on the training fold only and applied to the held-out
  fold. No global tuning.

Usage:
    python code/bench_bayesian_tau_penalty.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct
from bench_tau2 import per_channel_logLR_table, compute_tau2

FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')

CHANNELS = [
    ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal',
                transform=logit_upper_half),
    ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),
    ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                higher_means_tp=True, tp_family='student_t', fp_family='student_t'),
]

RT_UNRELIABLE_NORMALIZED_ADDUCTS = {'Cat', 'Anion'}

# Sweep α on a coarse-then-fine grid
ALPHA_GRID = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0]


def mask_rt_for_unreliable_adducts(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    norm = out['adduct'].apply(norm_adduct)
    mask = norm.isin(RT_UNRELIABLE_NORMALIZED_ADDUCTS)
    out.loc[mask, 'signed_delta_rt'] = np.nan
    return out


def fit_fold(top1_train_tr, labels_tr):
    """Fit the per-channel Bayesian densities on a training fold.
    Returns (fitted_channels, sigma2_k_array, prior_log_odds)."""
    fitted = [fit_channel(s, top1_train_tr[s.feature_col].values, labels_tr)
              for s in CHANNELS]
    lr_tr = per_channel_logLR_table(top1_train_tr, fitted)
    sigma2 = np.array([np.nanvar(lr_tr[:, k]) if np.isfinite(lr_tr[:, k]).any() else 0.0
                       for k in range(len(CHANNELS))])
    sigma2 = np.where(sigma2 > 0, sigma2, 1e-6)
    prior_tp = labels_tr.mean()
    prior_log_odds = np.log(prior_tp / (1 - prior_tp))
    return fitted, sigma2, prior_log_odds


def score_with_penalty(df, fitted, sigma2, prior_log_odds, alpha):
    """Return (posterior, total_lr, tau2) for rows of `df` under penalty α."""
    lr_table = per_channel_logLR_table(df, fitted)
    # logLR contribution: NaN channels → 0 (correct: missing evidence)
    lr_safe = np.where(np.isfinite(lr_table), lr_table, 0.0)
    total_lr = lr_safe.sum(axis=1)
    tau2 = compute_tau2(lr_table, sigma2)
    log_post_odds = prior_log_odds + total_lr - alpha * tau2
    posterior = 1.0 / (1.0 + np.exp(-log_post_odds))
    return posterior, total_lr, tau2


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


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows × {len(ft.columns)} cols, {ft["wiki_id"].nunique():,} spectra')

    print('\nApplying RT mask for [Cat]+/[Anion]- adducts...')
    ft = mask_rt_for_unreliable_adducts(ft)

    labeled_mask = ft['spectrum_label'].isin(['TP', 'FP'])
    ft_labeled = ft[labeled_mask]
    top1_idx = ft_labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    labels = top1['hit_label'].values
    print(f'\nTraining set: {len(top1):,} top-1 rows (TP+FP)   prior={labels.mean():.3f}')

    # GroupKFold by anno_ik14
    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'
    gkf = GroupKFold(n_splits=5)

    # We collect per-α OOF posteriors so we can score each α globally
    n = len(top1)
    oof_baseline = np.full(n, np.nan)
    oof_by_alpha = {a: np.full(n, np.nan) for a in ALPHA_GRID}
    fold_per_alpha = {a: [] for a in ALPHA_GRID}

    print('\nFitting per-fold...')
    for fold, (tr, te) in enumerate(gkf.split(top1, labels, groups)):
        fitted, sigma2, prior_log_odds = fit_fold(top1.iloc[tr], labels[tr])
        # Baseline (α=0 by definition, but we keep an explicit copy for clarity)
        post0, _, tau2_te = score_with_penalty(top1.iloc[te], fitted, sigma2,
                                                prior_log_odds, alpha=0.0)
        oof_baseline[te] = post0

        # All α values
        for a in ALPHA_GRID:
            post_a, _, _ = score_with_penalty(top1.iloc[te], fitted, sigma2,
                                              prior_log_odds, alpha=a)
            oof_by_alpha[a][te] = post_a
            fold_per_alpha[a].append(roc_auc_score(labels[te], post_a))

        print(f'  Fold {fold}: '
              f'baseline AUC={roc_auc_score(labels[te], post0):.4f}  '
              f'τ² range [{tau2_te.min():.2f}, {tau2_te.max():.2f}]  '
              f'mean τ²={tau2_te.mean():.2f}')

    # Per-α global OOF AUC and Brier (raw, then isotonic-calibrated)
    print('\n=== α sweep — global OOF metrics ===')
    print(f'{"α":>5}  {"AUC":>7}  {"Brier_raw":>10}  {"Brier_cal":>10}  {"ECE_cal":>8}')
    results = []
    for a in ALPHA_GRID:
        p = oof_by_alpha[a]
        v = ~np.isnan(p)
        auc = roc_auc_score(labels[v], p[v])
        br_raw = brier_score_loss(labels[v], p[v])
        iso = IsotonicRegression(out_of_bounds='clip').fit(p[v], labels[v])
        p_cal = iso.transform(p[v])
        br_cal = brier_score_loss(labels[v], p_cal)
        ece = ece_score(p_cal, labels[v])
        results.append((a, auc, br_raw, br_cal, ece))
        marker = '  ←' if a == 0.0 else ''
        print(f'{a:>5.2f}  {auc:>7.4f}  {br_raw:>10.4f}  {br_cal:>10.4f}  {ece:>8.4f}{marker}')

    # Identify the best α by AUC (with tie-breaker: lower Brier_cal)
    best = max(results, key=lambda r: (r[1], -r[3]))
    baseline = results[0]
    print(f'\nBaseline (α=0): AUC={baseline[1]:.4f}  Brier_cal={baseline[3]:.4f}')
    print(f'Best α={best[0]:.2f}: AUC={best[1]:.4f}  Brier_cal={best[3]:.4f}')
    print(f'Δ vs baseline: AUC={best[1]-baseline[1]:+.4f}  Brier_cal={best[3]-baseline[3]:+.4f}')

    # Honest read on α stability
    print('\n=== Per-fold AUCs at best α ===')
    fold_aucs_best = fold_per_alpha[best[0]]
    print(f'  α={best[0]:.2f}: folds = {[f"{a:.4f}" for a in fold_aucs_best]}')
    print(f'  range [{min(fold_aucs_best):.4f}, {max(fold_aucs_best):.4f}]')

    # Compare TP / FP τ² distribution explicitly to confirm signal direction
    fitted_full, sigma2_full, _ = fit_fold(top1, labels)
    _, _, tau2_full = score_with_penalty(top1, fitted_full, sigma2_full, 0.0, alpha=0.0)
    tp_tau = tau2_full[labels == 1]
    fp_tau = tau2_full[labels == 0]
    print('\n=== τ² distribution (full-fit) ===')
    print(f'  TP rows: mean={tp_tau.mean():.3f}  median={np.median(tp_tau):.3f}  '
          f'95p={np.percentile(tp_tau, 95):.3f}')
    print(f'  FP rows: mean={fp_tau.mean():.3f}  median={np.median(fp_tau):.3f}  '
          f'95p={np.percentile(fp_tau, 95):.3f}')


if __name__ == '__main__':
    main()
