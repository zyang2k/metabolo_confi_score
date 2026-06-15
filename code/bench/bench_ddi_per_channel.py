"""
bench_ddi_per_channel.py — Per-channel TP-vs-FP separability via Distribution
Discrepancy Index (Huang, Axioms 2026).

DDI is a normalized, symmetric [0, 1] index:
    DDI(p, q) = 1 − BC(p, q)
where BC is the Bhattacharyya coefficient:
    BC(p, q) = ∫ √(p(y) · q(y)) dy

For two normals N(μ₁, σ₁²) and N(μ₂, σ₂²) the closed form is:
    BC = √(2σ₁σ₂ / (σ₁² + σ₂²)) · exp[−¼ · (μ₁−μ₂)² / (σ₁² + σ₂²)]

For Student-t (and any non-Gaussian fit), we integrate numerically on a wide
grid covering both densities' support.

Compared to per-channel AUC:
  - DDI is calibration-invariant in the same sense AUC is (rank-only),
  - BUT it's a property of the FITTED density — independent of how many
    samples land in the test set,
  - AND it satisfies normalization + symmetry (KL and χ² do not).

Output:
  - data/ddi_per_channel.csv  (DDI + per-channel AUC for comparison)
  - figures/ddi_per_channel.png  (bar chart, deck-ready)

Usage:
    python code/bench_ddi_per_channel.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct

FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
OUT_CSV = os.path.join(ROOT, 'data', 'ddi_per_channel.csv')
OUT_FIG = os.path.join(ROOT, 'figures', 'ddi_per_channel.png')

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


def mask_rt_for_unreliable_adducts(df):
    out = df.copy()
    norm = out['adduct'].apply(norm_adduct)
    mask = norm.isin(RT_UNRELIABLE_NORMALIZED_ADDUCTS)
    out.loc[mask, 'signed_delta_rt'] = np.nan
    return out


def bc_two_normals(mu1, s1, mu2, s2):
    """Bhattacharyya coefficient for two univariate normals (closed form)."""
    var = s1**2 + s2**2
    return np.sqrt(2 * s1 * s2 / var) * np.exp(-0.25 * (mu1 - mu2)**2 / var)


def bc_numeric(p_dist, q_dist, lo, hi, n=20000):
    """BC = ∫ √(p·q) dy by trapezoid rule on [lo, hi]."""
    y = np.linspace(lo, hi, n)
    p = np.exp(p_dist.logpdf(y))
    q = np.exp(q_dist.logpdf(y))
    return np.trapz(np.sqrt(p * q), y)


def channel_ddi(fitted_channel, tp_vals, fp_vals):
    """Return (DDI, BC) for a fitted channel.

    Tries the closed-form normal-vs-normal first; falls back to numerical
    integration over a grid spanning both distributions' ±6σ region.
    """
    spec = fitted_channel.spec
    tp = fitted_channel.tp_dist
    fp = fitted_channel.fp_dist

    # Closed form when both are normals
    if spec.tp_family == 'normal' and spec.fp_family == 'normal':
        bc = bc_two_normals(tp.mean(), tp.std(), fp.mean(), fp.std())
    else:
        # Numerical integration; pick grid bounds covering both fits
        # Use the union of TP and FP empirical 0.1–99.9 percentiles, padded
        all_vals = np.concatenate([tp_vals[np.isfinite(tp_vals)],
                                   fp_vals[np.isfinite(fp_vals)]])
        if spec.transform is not None:
            all_vals = spec.transform(all_vals)
            all_vals = all_vals[np.isfinite(all_vals)]
        if len(all_vals) < 10:
            return np.nan, np.nan
        lo = np.percentile(all_vals, 0.1)
        hi = np.percentile(all_vals, 99.9)
        # Pad
        span = hi - lo if hi > lo else 1.0
        lo -= 0.5 * span
        hi += 0.5 * span
        bc = bc_numeric(tp, fp, lo, hi)

    bc = float(np.clip(bc, 0.0, 1.0))
    ddi = 1.0 - bc
    return ddi, bc


def channel_auc(top1, channel_spec, labels):
    """Per-channel rank-AUC of the raw feature value, ignoring missing rows."""
    vals = top1[channel_spec.feature_col].values.astype(float)
    if channel_spec.transform is not None:
        vals = channel_spec.transform(vals)
    valid = np.isfinite(vals)
    if valid.sum() < 10:
        return np.nan
    if not channel_spec.higher_means_tp:
        vals = -vals
    return roc_auc_score(labels[valid], vals[valid])


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    ft = mask_rt_for_unreliable_adducts(ft)

    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1 = ft.loc[labeled.groupby('wiki_id')['entropy_similarity'].idxmax()].reset_index(drop=True)
    labels = top1['hit_label'].values

    print(f'Top-1 training rows: {len(top1):,}   prior_TP={labels.mean():.3f}\n')

    # Fit channels on full top-1 set
    fitted_channels = []
    for spec in CHANNELS:
        fc = fit_channel(spec, top1[spec.feature_col].values, labels)
        fitted_channels.append(fc)

    rows = []
    for spec, fc in zip(CHANNELS, fitted_channels):
        # Get TP/FP raw values for this channel (used for AUC + numeric grid)
        v = top1[spec.feature_col].values.astype(float)
        tp_v = v[labels == 1]
        fp_v = v[labels == 0]

        ddi, bc = channel_ddi(fc, tp_v, fp_v)
        auc = channel_auc(top1, spec, labels)

        # Pull out fitted parameters for the table
        tp_mu = fc.tp_dist.mean() if fc.tp_dist is not None else np.nan
        tp_sd = fc.tp_dist.std() if fc.tp_dist is not None else np.nan
        fp_mu = fc.fp_dist.mean() if fc.fp_dist is not None else np.nan
        fp_sd = fc.fp_dist.std() if fc.fp_dist is not None else np.nan

        rows.append({
            'channel': spec.name,
            'tp_family': spec.tp_family,
            'fp_family': spec.fp_family,
            'tp_n': fc.tp_n,
            'fp_n': fc.fp_n,
            'tp_mu': tp_mu, 'tp_sd': tp_sd,
            'fp_mu': fp_mu, 'fp_sd': fp_sd,
            'BC': bc,
            'DDI': ddi,
            'AUC': auc,
        })
        print(f'  {spec.name:18s}  TP=N({tp_mu:+.2f}, {tp_sd:.2f})  '
              f'FP=N({fp_mu:+.2f}, {fp_sd:.2f})  '
              f'BC={bc:.3f}  DDI={ddi:.3f}  AUC={auc:.3f}')

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f'\nWrote {OUT_CSV}')

    # Plot — two panels because DDI [0, 1] and AUC [0.5, 1] have different scales.
    # Highlight that DDI ranks the channels differently than rank-AUC.
    os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(df))

    # Left: DDI
    ax = axes[0]
    ddi_rank = df['DDI'].rank(ascending=False).astype(int).tolist()
    bars = ax.bar(x, df['DDI'], color='#3a6fa5')
    ax.set_xticks(x); ax.set_xticklabels(df['channel'])
    ax.set_ylim(0, max(0.25, df['DDI'].max() * 1.4))
    ax.set_ylabel('DDI  (1 − Bhattacharyya coefficient)')
    ax.set_title('DDI — Distribution Discrepancy Index\n(calibration-invariant, property of fitted PDFs)')
    for i, (ddi, r) in enumerate(zip(df['DDI'], ddi_rank)):
        ax.text(i, ddi + 0.005, f'{ddi:.3f}\n(rank {r})', ha='center', fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    # Right: per-channel AUC
    ax = axes[1]
    auc_rank = df['AUC'].rank(ascending=False).astype(int).tolist()
    ax.bar(x, df['AUC'], color='#c46a3a')
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, label='chance (0.5)')
    ax.set_xticks(x); ax.set_xticklabels(df['channel'])
    ax.set_ylim(0.45, 1.0)
    ax.set_ylabel('per-channel ROC AUC')
    ax.set_title('Rank-AUC of raw feature\n(rank-based, ignores distribution shape)')
    for i, (auc, r) in enumerate(zip(df['AUC'], auc_rank)):
        ax.text(i, auc + 0.005, f'{auc:.3f}\n(rank {r})', ha='center', fontsize=9)
    ax.grid(axis='y', alpha=0.3)
    ax.legend(loc='lower right')

    fig.suptitle('Per-channel TP-vs-FP separability — DDI vs rank-AUC',
                 fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_FIG, dpi=140, bbox_inches='tight')
    print(f'Wrote {OUT_FIG}')


if __name__ == '__main__':
    main()
