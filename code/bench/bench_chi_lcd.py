"""
bench_chi_lcd.py — T3 sanity check: is the Huang information-loss diagnostic
I² (and χ_LCD) distinct from `ensemble_sd`, or redundant?

Definitions (Huang Metrology 2026):
  σ²_WA(row) = K_valid / Σ_{k∈valid}(1/σ²_k)
  I²(row)    = τ² / (τ² + σ²_WA)               ∈ [0, 1]
  χ_LCD(row) = 1 − √( σ²_WA / (K_valid·(σ²_WA + τ²)) )

τ² is the inverse-variance-weighted variance of per-row, per-channel logLRs
around their pooled mean — exactly what `score_confidence_v2.py` already ships
as the `tau2` column on `deliverable_scores_bayesian_v2.csv`.

This bench:
  1. Fits the 3-channel Bayesian on the labeled top-1 rows (no fold split — this
     is a *diagnostic* on the production-fit row order, not an OOF benchmark).
  2. Computes per-row τ², σ²_WA, I², χ_LCD on the scored picks
     (`pick_source` in {curator_ik14, curator_name_fallback, top1_by_sim}).
  3. Joins to `deliverable_scores_v2.csv` to grab `ensemble_sd` (GBM bootstrap
     epistemic uncertainty).
  4. Reports:
      - Spearman ρ between I² and ensemble_sd on the labeled slice.
        Threshold for "non-redundant" is ρ < 0.5.
      - Distribution of τ²/σ²_WA across rows (median, 90p, fraction > 0.1).
        If 90% of rows have τ²/σ²_WA < 0.1, I² is floor-bound and not useful.

Usage:
    python code/bench_chi_lcd.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct
from score_confidence_v2 import (
    _per_channel_logLR_table, _compute_tau2, _compute_chi_lcd,
    _channel_sigma2, mask_rt_for_unreliable_adducts, pick_scored_rows,
)

FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
GBM_DELIVERABLE = os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv')

CHANNELS = [
    ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal',
                transform=logit_upper_half),
    ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                higher_means_tp=True, tp_family='normal', fp_family='normal'),
    ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                higher_means_tp=True, tp_family='student_t', fp_family='student_t'),
]


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows × {len(ft.columns)} cols')

    print('\nMasking RT for [Cat]+/[Anion]- adducts (matches production scorer)...')
    ft = mask_rt_for_unreliable_adducts(ft)

    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    top1_idx = labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1_train = ft.loc[top1_idx].reset_index(drop=True)
    labels = top1_train['hit_label'].values
    print(f'  Top-1 training rows: {len(top1_train):,}   prior(TP)={labels.mean():.3f}')

    print('\nFitting 3-channel Bayesian on top-1 (full-data fit, same as production)...')
    fitted = [fit_channel(s, top1_train[s.feature_col].values, labels) for s in CHANNELS]
    lr_train = _per_channel_logLR_table(top1_train, fitted)
    sigma2 = _channel_sigma2(lr_train)
    print(f'  σ²_k per channel: '
          f'entropy_sim={sigma2[0]:.3f}, sim_gap={sigma2[1]:.3f}, signed_delta_rt={sigma2[2]:.3f}')

    # Compute diagnostics over all candidates (so we can join to picks afterward).
    print('\nComputing τ², σ²_WA, I², χ_LCD on every candidate row...')
    lr_full = _per_channel_logLR_table(ft, fitted)
    tau2 = _compute_tau2(lr_full, sigma2)
    sigma2_WA, I_sq, chi_lcd = _compute_chi_lcd(lr_full, sigma2, tau2)

    ft = ft.copy()
    ft['tau2'] = tau2
    ft['sigma2_WA'] = sigma2_WA
    ft['I_squared'] = I_sq
    ft['chi_lcd'] = chi_lcd

    # Pick scored rows (one per spectrum) the same way production does.
    picks = pick_scored_rows(ft)
    print(f'  Picks: {len(picks):,} spectra (one row per wiki_id)')

    # --- T3a: distribution of τ²/σ²_WA across rows ---------------------------
    print('\n=== τ²/σ²_WA distribution (across all candidates) ===')
    ratio = ft['tau2'] / ft['sigma2_WA']
    ratio_finite = ratio[np.isfinite(ratio)]
    print(f'  finite rows: {len(ratio_finite):,} / {len(ft):,}')
    print(f'  median = {ratio_finite.median():.4f}')
    print(f'  90p    = {ratio_finite.quantile(0.90):.4f}')
    print(f'  99p    = {ratio_finite.quantile(0.99):.4f}')
    print(f'  fraction > 0.1: {(ratio_finite > 0.1).mean():.3f}')
    print(f'  fraction > 0.5: {(ratio_finite > 0.5).mean():.3f}')
    print(f'  fraction > 1.0: {(ratio_finite > 1.0).mean():.3f}')

    # Same on picks only (the row that actually gets a confidence score).
    pick_ratio = (picks['tau2'] / picks['sigma2_WA'])
    pick_ratio = pick_ratio[np.isfinite(pick_ratio)]
    print(f'\n  on picks: median={pick_ratio.median():.4f}  '
          f'90p={pick_ratio.quantile(0.90):.4f}  '
          f'fraction>0.1: {(pick_ratio > 0.1).mean():.3f}')

    # --- T3b: I² vs ensemble_sd correlation ---------------------------------
    print(f'\nReading GBM deliverable for ensemble_sd: {GBM_DELIVERABLE}')
    gbm = pd.read_csv(GBM_DELIVERABLE, low_memory=False)
    print(f'  {len(gbm):,} rows')

    # Merge on wiki_id (one row per spectrum on both sides).
    merged = picks[['wiki_id', 'spectrum_label', 'hit_label',
                    'tau2', 'sigma2_WA', 'I_squared', 'chi_lcd']].merge(
        gbm[['wiki_id', 'ensemble_sd', 'confidence']].rename(
            columns={'confidence': 'gbm_confidence'}),
        on='wiki_id', how='inner'
    )
    print(f'  After join: {len(merged):,} spectra')

    labeled_slice = merged[merged['spectrum_label'].isin(['TP', 'FP'])].copy()
    print(f'  Labeled (TP+FP) slice: {len(labeled_slice):,}')

    # T3 verdict
    print('\n=== T3: Spearman correlation (labeled slice) ===')
    for var in ['I_squared', 'chi_lcd', 'tau2']:
        valid = np.isfinite(labeled_slice[var]) & np.isfinite(labeled_slice['ensemble_sd'])
        rho, p = spearmanr(labeled_slice.loc[valid, var], labeled_slice.loc[valid, 'ensemble_sd'])
        n = valid.sum()
        verdict = 'NON-REDUNDANT' if abs(rho) < 0.5 else ('weak overlap' if abs(rho) < 0.7 else 'REDUNDANT')
        print(f'  {var:>10s} ↔ ensemble_sd   ρ={rho:+.4f}  p={p:.2e}  n={n:,}   → {verdict}')

    # Also: I² vs gbm_confidence and ensemble_sd vs gbm_confidence as anchors
    print('\n=== Anchor correlations (sanity) ===')
    for left, right in [('I_squared', 'gbm_confidence'),
                        ('ensemble_sd', 'gbm_confidence'),
                        ('tau2', 'gbm_confidence')]:
        valid = np.isfinite(labeled_slice[left]) & np.isfinite(labeled_slice[right])
        rho, _ = spearmanr(labeled_slice.loc[valid, left], labeled_slice.loc[valid, right])
        print(f'  {left:>12s} ↔ {right:<16s}  ρ={rho:+.4f}  n={valid.sum():,}')

    # Distribution of I² in labeled slice
    print('\n=== I² distribution on labeled slice ===')
    for lab in [1, 0]:
        s = labeled_slice.loc[labeled_slice['hit_label'] == lab, 'I_squared']
        print(f'  {"TP" if lab == 1 else "FP"} rows (n={len(s):,}): '
              f'median={s.median():.4f}  90p={s.quantile(0.9):.4f}  '
              f'fraction>0.1: {(s > 0.1).mean():.3f}')

    # Dump for follow-up (T1, T2, T4 can read this)
    out_path = os.path.join(ROOT, 'data', 'bench_chi_lcd_join.csv')
    labeled_slice.to_csv(out_path, index=False)
    print(f'\nWrote {out_path}: {len(labeled_slice):,} rows')


if __name__ == '__main__':
    main()
