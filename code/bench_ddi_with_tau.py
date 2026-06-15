"""
bench_ddi_with_tau.py — Add τ² to the DDI / AUC per-channel comparison.

τ² is a derived quantity (between-channel disagreement), not a raw evidence
feature, but for diagnostic purposes we can ask: how does τ² compare to the
three input channels as a TP-vs-FP signal?

τ² is non-negative and right-skewed; fitting a normal is wrong. We use KDE for
the TP and FP τ² densities, then numerically integrate BC = ∫√(p_TP·p_FP) dy.

Usage:
    python code/bench_ddi_with_tau.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct
from bench_tau2 import per_channel_logLR_table, compute_tau2

FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
OUT_FIG = os.path.join(ROOT, 'figures', 'ddi_with_tau.png')

CHANNELS = [
    ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal',
                transform=logit_upper_half),
    ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),
    ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                higher_means_tp=True, tp_family='student_t', fp_family='student_t'),
]


def mask_rt(df):
    out = df.copy()
    norm = out['adduct'].apply(norm_adduct)
    mask = norm.isin({'Cat', 'Anion'})
    out.loc[mask, 'signed_delta_rt'] = np.nan
    return out


def bc_two_normals(mu1, s1, mu2, s2):
    var = s1**2 + s2**2
    return np.sqrt(2 * s1 * s2 / var) * np.exp(-0.25 * (mu1 - mu2)**2 / var)


def bc_kde(tp_vals, fp_vals, lo, hi, n=20000):
    """BC via KDE for non-Gaussian samples (e.g. τ², which is right-skewed)."""
    p_kde = stats.gaussian_kde(tp_vals)
    q_kde = stats.gaussian_kde(fp_vals)
    y = np.linspace(lo, hi, n)
    return float(np.trapz(np.sqrt(p_kde(y) * q_kde(y)), y))


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    ft = mask_rt(ft)
    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1 = ft.loc[labeled.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    labels = top1['hit_label'].values
    print(f'Top-1: {len(top1):,}  prior_TP={labels.mean():.3f}\n')

    fitted = [fit_channel(s, top1[s.feature_col].values, labels) for s in CHANNELS]

    rows = []

    # Channels 1-3 — closed form for normals, numeric for student-t
    for spec, fc in zip(CHANNELS, fitted):
        v = top1[spec.feature_col].values.astype(float)
        if spec.tp_family == 'normal' and spec.fp_family == 'normal':
            bc = bc_two_normals(fc.tp_dist.mean(), fc.tp_dist.std(),
                                fc.fp_dist.mean(), fc.fp_dist.std())
        else:
            tp_v = v[(labels == 1) & np.isfinite(v)]
            fp_v = v[(labels == 0) & np.isfinite(v)]
            all_v = np.concatenate([tp_v, fp_v])
            lo, hi = np.percentile(all_v, [0.1, 99.9])
            span = hi - lo
            bc = bc_kde(tp_v, fp_v, lo - 0.5 * span, hi + 0.5 * span)
        ddi = 1.0 - bc

        # Per-channel AUC of raw value
        valid = np.isfinite(v)
        x = v.copy()
        if spec.transform is not None:
            x = spec.transform(x)
            valid = valid & np.isfinite(x)
        if not spec.higher_means_tp:
            x = -x
        auc = roc_auc_score(labels[valid], x[valid])
        rows.append((spec.name, ddi, auc))
        print(f'  {spec.name:20s}  DDI={ddi:.3f}  AUC={auc:.3f}')

    # Channel 4 — τ² (derived from the three above)
    lr_full = per_channel_logLR_table(top1, fitted)
    sigma2 = np.array([np.nanvar(lr_full[:, k]) for k in range(len(CHANNELS))])
    sigma2 = np.where(sigma2 > 0, sigma2, 1e-6)
    tau2 = compute_tau2(lr_full, sigma2)

    tp_tau = tau2[labels == 1]
    fp_tau = tau2[labels == 0]
    # Grid: τ² is non-negative and right-skewed; use 0 to TP/FP combined 99.5%ile + pad
    hi = np.percentile(np.concatenate([tp_tau, fp_tau]), 99.5) * 1.5
    bc_tau = bc_kde(tp_tau, fp_tau, 0.0, hi)
    ddi_tau = 1.0 - bc_tau
    # τ² is an INVERSE predictor (higher = less TP-like)
    auc_tau = roc_auc_score(labels, -tau2)
    rows.append(('tau2 (derived)', ddi_tau, auc_tau))
    print(f'  {"tau2 (derived)":20s}  DDI={ddi_tau:.3f}  AUC={auc_tau:.3f}  [inverse predictor]')
    print(f'    TP τ²: mean={tp_tau.mean():.3f}  median={np.median(tp_tau):.3f}')
    print(f'    FP τ²: mean={fp_tau.mean():.3f}  median={np.median(fp_tau):.3f}')

    df = pd.DataFrame(rows, columns=['channel', 'DDI', 'AUC'])
    df['DDI_rank'] = df['DDI'].rank(ascending=False).astype(int)
    df['AUC_rank'] = df['AUC'].rank(ascending=False).astype(int)
    print('\n' + df.to_string(index=False))

    # Plot — same two-panel layout, now with 4 bars including τ²
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    x = np.arange(len(df))
    colors = ['#3a6fa5'] * 3 + ['#7a4ca8']  # τ² bar in purple to distinguish
    colors_auc = ['#c46a3a'] * 3 + ['#a07050']

    ax = axes[0]
    ax.bar(x, df['DDI'], color=colors)
    ax.set_xticks(x); ax.set_xticklabels(df['channel'], rotation=15, ha='right')
    ax.set_ylim(0, max(0.25, df['DDI'].max() * 1.4))
    ax.set_ylabel('DDI  (1 − Bhattacharyya coefficient)')
    ax.set_title('DDI — model-side TP/FP separability')
    for i, (ddi, r) in enumerate(zip(df['DDI'], df['DDI_rank'])):
        ax.text(i, ddi + 0.005, f'{ddi:.3f}\n(rank {r})', ha='center', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    ax = axes[1]
    ax.bar(x, df['AUC'], color=colors_auc)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='chance')
    ax.set_xticks(x); ax.set_xticklabels(df['channel'], rotation=15, ha='right')
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel('per-channel ROC AUC')
    ax.set_title('Rank-AUC (test-sample dependent, rank-only)')
    for i, (auc, r) in enumerate(zip(df['AUC'], df['AUC_rank'])):
        ax.text(i, auc + 0.005, f'{auc:.3f}\n(rank {r})', ha='center', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='lower right')

    fig.suptitle('Per-channel separability — original 3 channels + τ² (derived)',
                 fontsize=12, y=1.02)
    plt.tight_layout()
    os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)
    fig.savefig(OUT_FIG, dpi=140, bbox_inches='tight')
    print(f'\nWrote {OUT_FIG}')


if __name__ == '__main__':
    main()
