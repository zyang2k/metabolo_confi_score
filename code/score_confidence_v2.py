"""
score_confidence_v2.py — Per-bin Bayesian confidence scoring for annotation delivery.

Two output sheets:
  * `data/deliverable_scores_v2.csv` — client-facing, one row per bin.
      Columns: annotation, confidence_pct, n_close_alternatives.
      For curator-annotated bins (spectrum_label='TP' with non-empty anno_ik14),
      the scored row is the one matching anno_ik14. Otherwise, the top-1 by
      entropy_similarity is scored. This way the client sees the curator's
      compound — not a pipeline-picked override — whenever the curator
      expressed a preference.
  * `data/curator_review_v2.csv` — internal review sheet.
      Surfaces pipeline-vs-curator disagreements: pipeline_pick, pipeline_conf,
      oliver_pick, agrees, score_on_oliver_pick, and a simple `flag`.

Scoring fixes applied here (consistent with the top-1 training):
  1. `signed_delta_rt` NaN passthrough — no fillna(0) bias for missing RT.
  2. `signed_delta_rt` masked to NaN for permanent-cation adducts
     (`[Cat]+`, `[Anion]-`). Fanzhou's RT predictor is unreliable for cations.
  3. τ² penalty (Huang Metrology 2026): logit(P) ← prior + Σ_k logLR_k − α·τ²,
     where τ² is the inverse-variance weighted variance of the per-channel logLRs.
     α=0.5 validated in code/bench_bayesian_tau_penalty.py (+0.0093 OOF AUC,
     −0.0038 Brier). Set ALPHA_TAU=0 to disable.

Design decisions locked in `project_deliverable_shape_20260422.md`:
  - Scoring unit = bin = spectrum (each wiki_id is one bin).
  - `n_close_alternatives` uses entropy_similarity, NOT posterior — avoids the
    rank-bias of applying the top-1-trained model to non-top-1 rows.

Usage:
    python code/score_confidence_v2.py
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold
from sklearn.isotonic import IsotonicRegression

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'code'))

from bayesian_score_v2 import ChannelSpec, fit_channel, logit_upper_half
from build_features_v2 import norm_adduct

FEATURE_TABLE    = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
# Bayesian output paths — companion to the primary GBM scorer (see score_gbm_v2.py).
# GBM is the primary production model; Bayesian is kept for interpretability / diagnostic.
DELIVERABLE_OUT  = os.path.join(ROOT, 'data', 'deliverable_scores_bayesian_v2.csv')
CANDIDATES_OUT   = os.path.join(ROOT, 'data', 'candidate_scores_bayesian_v2.csv')
REVIEW_OUT       = os.path.join(ROOT, 'data', 'curator_review_v2.csv')

# Entropy-similarity distance within which another IK14 counts as a "close alternative".
# Chosen to roughly match Oliver's "within 10 percentage points" framing from his email,
# interpreted on the entropy_sim [0, 1] scale. Not posterior-based — see module docstring.
ALT_THRESHOLD_SIM = 0.10

# Adducts whose RT prediction is unreliable (Fanzhou's model is trained on neutrals).
# Expand cautiously — every adduct listed here will drop the RT evidence for that hit.
RT_UNRELIABLE_NORMALIZED_ADDUCTS = {'Cat', 'Anion'}

# Confidence threshold above which a pipeline-vs-curator disagreement is "confident".
CONFIDENT_DISAGREE_THRESHOLD = 0.70

# τ² penalty coefficient (Huang Metrology 2026). Set to 0 to revert to plain LCD posterior.
ALPHA_TAU = 0.5


def _per_channel_logLR_table(df, fitted_channels):
    """(n_rows, n_channels) per-channel logLR; NaN where the input feature is non-finite."""
    n = len(df)
    K = len(fitted_channels)
    out = np.full((n, K), np.nan, dtype=float)
    for k, fc in enumerate(fitted_channels):
        raw = df[fc.spec.feature_col].values.astype(float)
        valid = np.isfinite(raw)
        if valid.sum() == 0:
            continue
        lr = fc.logLR(raw)
        out[valid, k] = lr[valid]
    return out


def _compute_tau2(lr_matrix, sigma2):
    """Per-row inverse-variance weighted variance of per-channel logLRs.
    NaN-channel rows are dropped from that row's τ². Rows with <2 valid channels
    return τ²=0 (no disagreement defined)."""
    n, K = lr_matrix.shape
    inv_sigma2 = np.where(sigma2 > 0, 1.0 / sigma2, 0.0)
    valid_mask = np.isfinite(lr_matrix)
    n_valid = valid_mask.sum(axis=1)
    lr_safe = np.where(valid_mask, lr_matrix, 0.0)
    w_unnorm = valid_mask * inv_sigma2[None, :]
    W = w_unnorm.sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu_L = (w_unnorm * lr_safe).sum(axis=1) / np.where(W > 0, W, 1.0)
        diff2 = (lr_safe - mu_L[:, None]) ** 2
        tau2 = (w_unnorm * diff2 * valid_mask).sum(axis=1) / np.where(W > 0, W, 1.0)
    tau2[n_valid < 2] = 0.0
    tau2[~np.isfinite(tau2)] = 0.0
    return tau2


def _compute_chi_lcd(lr_matrix, sigma2, tau2):
    """Per-row information-loss diagnostics from Huang Metrology 2026.

    Returns (sigma2_WA, I_squared, chi_lcd):
      σ²_WA(row) = K_valid / Σ_{k∈valid}(1/σ²_k)  — weighted-average within-channel variance
      I²(row)    = τ² / (τ² + σ²_WA)              — heterogeneity fraction in [0, 1]
      χ_LCD(row) = 1 − √(σ²_WA / (K_valid · (σ²_WA + τ²)))  — total LCD information loss

    Rows with <2 valid channels: σ²_WA = NaN, I² = 0, χ_LCD = NaN (LCD-vs-WLP not defined)."""
    inv_sigma2 = np.where(sigma2 > 0, 1.0 / sigma2, 0.0)
    valid_mask = np.isfinite(lr_matrix)
    n_valid = valid_mask.sum(axis=1).astype(float)
    inv_sigma2_sum = (valid_mask * inv_sigma2[None, :]).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        sigma2_WA = np.where(inv_sigma2_sum > 0, n_valid / inv_sigma2_sum, np.nan)
        I_sq = np.where(sigma2_WA > 0, tau2 / (tau2 + sigma2_WA), 0.0)
        chi_lcd = 1.0 - np.sqrt(sigma2_WA / (n_valid * (sigma2_WA + tau2)))
    too_few = n_valid < 2
    sigma2_WA[too_few] = np.nan
    I_sq[too_few] = 0.0
    chi_lcd[too_few] = np.nan
    return sigma2_WA, I_sq, chi_lcd


def _channel_sigma2(lr_matrix):
    """Empirical variance of each channel's logLR over a (training) row set.
    Floored at 1e-6 to keep inverse-variance weights finite."""
    K = lr_matrix.shape[1]
    sigma2 = np.array([np.nanvar(lr_matrix[:, k]) if np.isfinite(lr_matrix[:, k]).any() else 0.0
                       for k in range(K)])
    return np.where(sigma2 > 0, sigma2, 1e-6)


def mask_rt_for_unreliable_adducts(df: pd.DataFrame) -> pd.DataFrame:
    """Set `signed_delta_rt` to NaN for rows whose adduct normalizes into
    `RT_UNRELIABLE_NORMALIZED_ADDUCTS`. NaN means "no RT evidence" downstream."""
    out = df.copy()
    norm = out['adduct'].apply(norm_adduct)
    mask = norm.isin(RT_UNRELIABLE_NORMALIZED_ADDUCTS)
    n = mask.sum()
    if n > 0:
        out.loc[mask, 'signed_delta_rt'] = np.nan
    print(f'  Masked signed_delta_rt for {n:,} rows '
          f'(adducts: {sorted(RT_UNRELIABLE_NORMALIZED_ADDUCTS)})')
    return out


def fit_and_score(ft: pd.DataFrame, top1: pd.DataFrame, labels, channels,
                  alpha_tau: float = ALPHA_TAU):
    """Fit channels on top-1 (labeled) rows, then score every row of `ft`.

    Posterior includes the τ² penalty: logit(P) ← prior + Σ_k logLR_k − α·τ²,
    where σ²_k for τ² weighting is estimated from the channel logLR distribution
    on the training top-1 rows.
    """
    fitted = [fit_channel(s, top1[s.feature_col].values, labels) for s in channels]
    prior_tp = labels.mean()
    prior_log_odds = np.log(prior_tp / (1 - prior_tp))

    lr_train = _per_channel_logLR_table(top1, fitted)
    sigma2 = _channel_sigma2(lr_train)

    lr_full = _per_channel_logLR_table(ft, fitted)
    total_lr = np.where(np.isfinite(lr_full), lr_full, 0.0).sum(axis=1)
    tau2 = _compute_tau2(lr_full, sigma2)
    sigma2_WA, I_sq, chi_lcd = _compute_chi_lcd(lr_full, sigma2, tau2)

    log_post_odds = prior_log_odds + total_lr - alpha_tau * tau2
    posterior = 1.0 / (1.0 + np.exp(-log_post_odds))
    return fitted, posterior, prior_log_odds, tau2, I_sq, chi_lcd


def oof_top1(top1: pd.DataFrame, labels, channels, n_splits: int = 5,
             alpha_tau: float = ALPHA_TAU):
    """GroupKFold OOF posterior on the top-1 training rows. For AUC reporting.

    Per-fold no leakage: channels and σ²_k are re-fit on the training fold and
    applied to held-out rows.
    """
    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'
    prior_tp = labels.mean()
    prior_log_odds = np.log(prior_tp / (1 - prior_tp))
    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    for fold, (tr, te) in enumerate(GroupKFold(n_splits=n_splits).split(top1, labels, groups)):
        fitted = [fit_channel(s, top1.iloc[tr][s.feature_col].values, labels[tr])
                  for s in channels]
        lr_tr = _per_channel_logLR_table(top1.iloc[tr], fitted)
        sigma2 = _channel_sigma2(lr_tr)

        lr_te = _per_channel_logLR_table(top1.iloc[te], fitted)
        fold_lr = np.where(np.isfinite(lr_te), lr_te, 0.0).sum(axis=1)
        tau2_te = _compute_tau2(lr_te, sigma2)

        oof[te] = 1.0 / (1.0 + np.exp(-(prior_log_odds + fold_lr - alpha_tau * tau2_te)))
        fold_aucs.append(roc_auc_score(labels[te], oof[te]))
    return oof, fold_aucs


def pick_scored_rows(ft: pd.DataFrame) -> pd.DataFrame:
    """Return one row per wiki_id — the row we score + report as the annotation.

    Selection priority (highest first):
      1. **IK14 match**: curator-annotated bin where hit_ik14 == anno_ik14.
         Highest-sim row among matches.
      2. **Name fallback**: curator-annotated bin where an empty-IK14 row's
         name matches Oliver's annotated compound name. This rescues cases
         where the library deposit has MSMS but no SMILES — typically trusted
         libraries (NIST23) with incomplete metadata. Without this, the
         pipeline would fall back to top-1-by-sim and display a different
         compound than Oliver chose.
      3. **Top-1 fallback**: yy_/FP bins, unannotated bins, or bins where
         neither an IK14 match nor a name match succeeded. Pick top-1 by
         entropy_similarity.

    Rationale: for curator-annotated bins we display THE CURATOR'S compound and
    score it; the client should not see the pipeline overriding a curator
    decision. When pipeline top-1 disagrees with the curator, that disagreement
    is surfaced separately in the review sheet, not on the client-facing sheet.
    """
    ft = ft.copy()
    ft['hit_ik14'] = ft['hit_ik14'].fillna('')
    ft['anno_ik14'] = ft['anno_ik14'].fillna('')

    # Which bins are curator-annotated?
    has_curator = (ft['spectrum_label'] == 'TP') & (ft['anno_ik14'] != '')

    # Pass 1 — IK14 match
    ik14_match = ft[has_curator & (ft['hit_ik14'] == ft['anno_ik14'])]
    ik14_picks_idx = ik14_match.groupby('wiki_id')['entropy_similarity'].idxmax()
    ik14_picks = ft.loc[ik14_picks_idx].copy()
    ik14_picks['pick_source'] = 'curator_ik14'

    # Pass 2 — name fallback (only for bins not already covered by IK14 match).
    # Match hit name to curator's annotated name (both lowercased / stripped).
    # Only considers empty-IK14 rows — IK14-populated rows that don't match anno_ik14
    # are genuinely different compounds and shouldn't rescue anything.
    covered = set(ik14_picks['wiki_id'])
    hit_name_lower = ft['name'].fillna('').str.strip().str.lower()
    anno_name_lower = ft['anno_name_lower'].fillna('').str.strip().str.lower() \
        if 'anno_name_lower' in ft.columns else pd.Series('', index=ft.index)
    name_match = ft[
        has_curator &
        (~ft['wiki_id'].isin(covered)) &
        (ft['hit_ik14'] == '') &
        (anno_name_lower != '') &
        (hit_name_lower == anno_name_lower)
    ]
    name_picks_idx = name_match.groupby('wiki_id')['entropy_similarity'].idxmax()
    name_picks = ft.loc[name_picks_idx].copy()
    name_picks['pick_source'] = 'curator_name_fallback'

    # Pass 3 — top-1 by entropy_similarity for everything else
    covered.update(name_picks['wiki_id'])
    still_remaining = ft[~ft['wiki_id'].isin(covered)]
    top1_idx = still_remaining.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1_picks = ft.loc[top1_idx].copy()
    top1_picks['pick_source'] = 'top1_by_sim'

    picks = pd.concat([ik14_picks, name_picks, top1_picks], ignore_index=True)
    assert picks['wiki_id'].nunique() == len(picks), \
        f'Duplicate wiki_ids in pick_scored_rows: {len(picks) - picks["wiki_id"].nunique()} extra'
    return picks


def count_alternatives_by_sim(ft: pd.DataFrame, picks: pd.DataFrame,
                              threshold_sim: float = ALT_THRESHOLD_SIM) -> pd.DataFrame:
    """For each wiki_id, count candidates with different non-empty hit_ik14
    whose entropy_similarity is within `threshold_sim` of the scored pick.

    Symmetric: a candidate above OR below the pick by ≤ threshold_sim counts.
    Catches the case where a higher-sim candidate loses to Oliver's picked
    lower-sim candidate — it's still a legitimate alternative.
    """
    ft = ft.copy()
    ft['hit_ik14'] = ft['hit_ik14'].fillna('')
    pick_info = picks[['wiki_id', 'entropy_similarity', 'hit_ik14']].rename(
        columns={'entropy_similarity': '_pick_sim', 'hit_ik14': '_pick_ik14'})
    ft = ft.merge(pick_info, on='wiki_id', how='left')

    is_alt = (
        (ft['hit_ik14'] != '') &
        (ft['hit_ik14'] != ft['_pick_ik14']) &
        ((ft['entropy_similarity'] - ft['_pick_sim']).abs() <= threshold_sim)
    )
    alt_counts = (ft[is_alt].groupby('wiki_id').size()
                  .rename('n_close_alternatives').reset_index())
    # Include zero-alt spectra
    all_wids = pd.DataFrame({'wiki_id': ft['wiki_id'].unique()})
    alt_counts = all_wids.merge(alt_counts, on='wiki_id', how='left').fillna({'n_close_alternatives': 0})
    alt_counts['n_close_alternatives'] = alt_counts['n_close_alternatives'].astype(int)
    return alt_counts


def assemble_curator_review(ft: pd.DataFrame, picks: pd.DataFrame) -> pd.DataFrame:
    """Build the internal review sheet that surfaces pipeline-vs-curator disagreements.

    Only rows with curator annotations (spectrum_label='TP' AND anno_ik14 non-empty)
    are included — yy_ bins have no curator pick to compare against.
    """
    ft = ft.copy()
    ft['hit_ik14'] = ft['hit_ik14'].fillna('')
    ft['anno_ik14'] = ft['anno_ik14'].fillna('')
    has_anno = (ft['spectrum_label'] == 'TP') & (ft['anno_ik14'] != '')
    anno_ft = ft[has_anno]

    # Pipeline top-1 (by entropy_sim) per spectrum
    top1_idx = anno_ft.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = anno_ft.loc[top1_idx][['wiki_id', 'name', 'hit_ik14', 'adduct',
                                  'entropy_similarity', 'candidate_posterior']].copy()
    top1 = top1.rename(columns={
        'name': 'pipeline_pick', 'hit_ik14': 'pipeline_ik14',
        'adduct': 'pipeline_adduct', 'entropy_similarity': 'pipeline_sim',
        'candidate_posterior': 'pipeline_conf',
    })

    # Curator's pick row (if present in candidates)
    curator_match = anno_ft[anno_ft['hit_ik14'] == anno_ft['anno_ik14']]
    curator_idx = curator_match.groupby('wiki_id')['entropy_similarity'].idxmax()
    curator = anno_ft.loc[curator_idx][['wiki_id', 'name', 'adduct',
                                        'entropy_similarity', 'candidate_posterior',
                                        'anno_ik14']].copy()
    curator = curator.rename(columns={
        'name': 'oliver_pick', 'adduct': 'oliver_adduct',
        'entropy_similarity': 'oliver_sim',
        'candidate_posterior': 'score_on_oliver_pick',
    })

    review = top1.merge(curator, on='wiki_id', how='left')
    review['agrees'] = (review['pipeline_ik14'] == review['anno_ik14']).astype(int)

    # Flag confidently-disagreeing calls (pipeline very sure, Oliver picked different)
    confident_disagree = (
        (review['agrees'] == 0)
        & (review['pipeline_conf'] >= CONFIDENT_DISAGREE_THRESHOLD)
    )
    review['flag'] = np.where(
        confident_disagree, 'pipeline_disagrees_confidently',
        np.where(review['oliver_pick'].isna(), 'oliver_pick_not_in_hits', '')
    )
    return review


def main():
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    ft['hit_ik14'] = ft['hit_ik14'].fillna('')
    print(f'  {len(ft):,} rows × {len(ft.columns)} cols, '
          f'{ft["wiki_id"].nunique():,} spectra')

    print('\nApplying RT fixes...')
    ft = mask_rt_for_unreliable_adducts(ft)

    # Top-1 training set — labeled spectra only (TP + FP). Blank spectra get
    # scored at inference but excluded from training because they have no
    # ground-truth hit_label (always 0 by construction since spectrum_label != 'TP').
    labeled_mask = ft['spectrum_label'].isin(['TP', 'FP'])
    ft_labeled = ft[labeled_mask]
    top1_train = (ft_labeled.loc[ft_labeled.groupby('wiki_id')['entropy_similarity'].idxmax()]
                    .reset_index(drop=True))
    labels = top1_train['hit_label'].values
    prior_tp = labels.mean()
    n_blank = (ft['spectrum_label'] == 'blank').nunique() if (ft['spectrum_label']=='blank').any() else 0
    n_blank_spectra = ft[ft['spectrum_label']=='blank']['wiki_id'].nunique()
    print(f'\nTop-1 training set (TP+FP): {len(top1_train):,} spectra   prior={prior_tp:.3f}')
    if n_blank_spectra:
        print(f'Blank spectra in table (scored but not trained on): {n_blank_spectra:,}')

    CHANNELS = [
        ChannelSpec('entropy_sim', 'entropy_similarity', 'continuous',
                    higher_means_tp=True, tp_family='normal', fp_family='normal',
                    transform=logit_upper_half),
        ChannelSpec('sim_gap', 'sim_gap', 'continuous',
                    higher_means_tp=True, tp_family='normal', fp_family='normal'),
        ChannelSpec('signed_delta_rt', 'signed_delta_rt', 'continuous',
                    higher_means_tp=True, tp_family='student_t', fp_family='student_t'),
    ]

    # OOF AUC on top-1 (honest, for reporting)
    print('\nFitting per-fold for OOF AUC on top-1...')
    oof, fold_aucs = oof_top1(top1_train, labels, CHANNELS)
    auc = roc_auc_score(labels[~np.isnan(oof)], oof[~np.isnan(oof)])
    print(f'  OOF AUC: {auc:.4f}   folds: [{min(fold_aucs):.4f}, {max(fold_aucs):.4f}]')

    # Full-data fit → score every candidate
    print(f'\nFitting full-data model (α_τ={ALPHA_TAU}); scoring every candidate...')
    _, candidate_posterior, _, candidate_tau2, candidate_I2, candidate_chi = fit_and_score(
        ft, top1_train, labels, CHANNELS)
    ft['candidate_posterior_raw'] = candidate_posterior
    ft['tau2'] = candidate_tau2
    ft['I_squared'] = candidate_I2
    ft['chi_lcd'] = candidate_chi
    print(f'  τ² range: [{candidate_tau2.min():.2f}, {candidate_tau2.max():.2f}]   '
          f'mean={candidate_tau2.mean():.2f}')
    valid_I2 = np.isfinite(candidate_I2)
    print(f'  I² range: [{candidate_I2[valid_I2].min():.3f}, {candidate_I2[valid_I2].max():.3f}]   '
          f'mean={candidate_I2[valid_I2].mean():.3f}   median={np.median(candidate_I2[valid_I2]):.3f}')
    valid_chi = np.isfinite(candidate_chi)
    floor = 1.0 - 1.0 / np.sqrt(3)  # 3-channel constant LCD floor at τ²=0
    print(f'  χ_LCD range: [{candidate_chi[valid_chi].min():.3f}, {candidate_chi[valid_chi].max():.3f}]   '
          f'mean={candidate_chi[valid_chi].mean():.3f}   '
          f'(3-channel floor at τ²=0 is {floor:.3f})')

    # ── Isotonic calibration ────────────────────────────────────────────────
    # Raw Bayesian posteriors are systematically under-confident across the
    # mid-range (~0.1–0.8). See project_min_validation_20260422.md and the
    # 2026-04-22 reliability-diagram run (ECE ~0.15). Fit isotonic on the honest
    # OOF posteriors of top-1 rows (labels = hit_label), apply to every row's
    # raw posterior to produce a client-facing calibrated score.
    print('\nFitting isotonic calibration on OOF top-1 posteriors...')
    oof_valid = ~np.isnan(oof)
    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(oof[oof_valid], labels[oof_valid])
    ft['candidate_posterior'] = iso.transform(ft['candidate_posterior_raw'].values)
    # Report calibration impact on the top-1 OOF values (honest self-check)
    oof_calibrated = iso.transform(oof[oof_valid])
    brier_raw = brier_score_loss(labels[oof_valid], oof[oof_valid])
    brier_cal = brier_score_loss(labels[oof_valid], oof_calibrated)
    print(f'  Brier: raw={brier_raw:.4f} → calibrated={brier_cal:.4f}')
    print(f'  AUC (calibration-invariant sanity): '
          f'{roc_auc_score(labels[oof_valid], oof_calibrated):.4f}')

    # Pick the row to score + report per spectrum
    print('\nSelecting scored row per spectrum...')
    picks = pick_scored_rows(ft)
    for src in ['curator_ik14', 'curator_name_fallback', 'top1_by_sim']:
        n = (picks['pick_source'] == src).sum()
        print(f'  {src:30s}  {n:,}')

    # Count close alternatives using entropy-similarity distance from the pick
    print(f'\nCounting close alternatives (within ±{ALT_THRESHOLD_SIM} entropy_sim of the pick)...')
    alts = count_alternatives_by_sim(ft, picks, threshold_sim=ALT_THRESHOLD_SIM)
    print(f'  n_close_alternatives: median={int(alts["n_close_alternatives"].median())}   '
          f'=0: {(alts["n_close_alternatives"]==0).sum():,}   '
          f'≥1: {(alts["n_close_alternatives"]>=1).sum():,}')

    # Attach alternatives + oof posterior + derive confidence columns
    out = picks.merge(alts, on='wiki_id', how='left')
    out['confidence'] = out['candidate_posterior']            # calibrated (client-facing)
    out['confidence_raw'] = out['candidate_posterior_raw']    # raw Bayesian (diagnostic)
    out['confidence_pct'] = (out['confidence'] * 100).round(1)
    top1_oof = pd.DataFrame({'wiki_id': top1_train['wiki_id'].values,
                             'oof_posterior_top1_model': oof})
    out = out.merge(top1_oof, on='wiki_id', how='left')

    deliverable_cols = [
        'wiki_id', 'spectrum_label', 'hit_label', 'pick_source',
        'name', 'adduct', 'hit_ik14', 'anno_ik14',
        'confidence', 'confidence_pct', 'confidence_raw',
        'oof_posterior_top1_model', 'n_close_alternatives',
        'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
        'tau2', 'I_squared', 'chi_lcd', 'hit_adduct_cat', 'db',
    ]
    deliverable = out[[c for c in deliverable_cols if c in out.columns]].rename(
        columns={'name': 'annotation'})
    deliverable.to_csv(DELIVERABLE_OUT, index=False)
    print(f'\nWrote {DELIVERABLE_OUT}: {len(deliverable):,} rows × {len(deliverable.columns)} cols')

    # Candidate-level dump (one row per (spectrum, candidate))
    cand_cols = ['wiki_id', 'hit_ik14', 'name', 'adduct', 'library_wiki_id',
                 'entropy_similarity', 'sim_gap', 'signed_delta_rt',
                 'candidate_posterior_raw', 'candidate_posterior', 'hit_label']
    ft[[c for c in cand_cols if c in ft.columns]].to_csv(CANDIDATES_OUT, index=False)
    print(f'Wrote {CANDIDATES_OUT}: {len(ft):,} rows (one per candidate)')

    # Curator review sheet
    print('\nAssembling curator review sheet...')
    review = assemble_curator_review(ft, picks)
    n_agree = (review['agrees'] == 1).sum()
    n_disagree_confident = (review['flag'] == 'pipeline_disagrees_confidently').sum()
    n_oliver_missing = (review['flag'] == 'oliver_pick_not_in_hits').sum()
    print(f'  Curator-annotated bins: {len(review):,}')
    print(f'    pipeline agrees: {n_agree:,} ({100*n_agree/max(1,len(review)):.1f}%)')
    print(f'    pipeline disagrees confidently (≥ {CONFIDENT_DISAGREE_THRESHOLD:.2f}): {n_disagree_confident:,}')
    print(f"    Oliver's pick not in candidate list: {n_oliver_missing:,}")
    review.to_csv(REVIEW_OUT, index=False)
    print(f'Wrote {REVIEW_OUT}: {len(review):,} rows × {len(review.columns)} cols')

    # Sanity check on deck cases (shows raw → calibrated)
    print('\n=== Deck cases (conf_raw → conf_calibrated) ===')
    for wid, label in [('aPUDE1U/2244', 'PEP slide 4'),
                       ('aEKJ9AS/10123', 'trig [Cat]+ slide 7'),
                       ('aEKJ9AS/994', 'deoxycarnitine [Cat]+ slide 8'),
                       ('aEKJ9AS/1674', 'trig baseline')]:
        r = deliverable[deliverable['wiki_id'] == wid]
        if len(r):
            r = r.iloc[0]
            print(f'  {wid:22s} ({label:28s}): '
                  f'source={r["pick_source"]:12s}  '
                  f'annotation="{r["annotation"][:25]:25s}"  '
                  f'raw={r["confidence_raw"]*100:4.0f}% → cal={r["confidence_pct"]:4.0f}%  '
                  f'n_alt={int(r["n_close_alternatives"])}')

    return deliverable, review


if __name__ == '__main__':
    main()
