"""add_noise_features.py — Add Oliver's 5 noise/quality features to feature_table_v2.csv.

Source: Fiehn CanMetCon Toronto deck 2026-04-25, slides 8–12.
Oliver's noise framework: normalized entropy is the single best separator;
99.93% of 11,363 annotated spectra had S_norm < 0.987. A one-class SVM with
5 features further refines noise classification.

Features added (per spectrum, joined to all candidate rows by wiki_id):
  - n_ions             : number of MS2 peaks
  - normalized_entropy : spectral_entropy / ln(n_ions)            [0, 1]
  - top5_pct           : sum(top-5 intensities) / sum(all)        [0, 1]
  - ion_density        : n_ions / (max_mz - min_mz)               peaks/Da
  - median_to_base     : median(intensity) / max(intensity)       [0, 1]

base_peak_intensity from Oliver's SVM is not added — raw_intensity column in
curated CSV is all zeros, and the peaks cache stores intensities already
normalized to max=1, so absolute base-peak amplitude is unrecoverable here.

Output: overwrites data/feature_table_v2.csv with the 5 new columns appended.
A backup is written to data/feature_table_v2_pre_noise_feats.csv.bak.
"""

from __future__ import annotations

import json
import os
import shutil

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
BACKUP = os.path.join(ROOT, 'data', 'feature_table_v2_pre_noise_feats.csv.bak')
PEAKS_CACHE = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')

NEW_FEATURES = ['n_ions', 'normalized_entropy', 'top5_pct',
                'ion_density', 'median_to_base']


def compute_per_spectrum_features(peaks: list[list[float]],
                                   spectral_entropy: float | None) -> dict:
    """Compute Oliver's 5 noise features from a peak list.

    Args:
        peaks: list of [mz, intensity] pairs. Intensities expected normalized
               to max=1 in the v2 cache, but we re-normalize defensively.
        spectral_entropy: pre-computed S = -Σ p_i ln p_i (already in feature table).
                          Used to derive normalized_entropy. If None, recompute.
    """
    n = len(peaks)
    if n == 0:
        return {f: np.nan for f in NEW_FEATURES}

    arr = np.asarray(peaks, dtype=np.float64)
    mz = arr[:, 0]
    inten = arr[:, 1]
    inten_max = inten.max()
    if inten_max <= 0:
        return {f: np.nan for f in NEW_FEATURES}

    # Renormalize intensity-to-base for top5 / median_to_base; entropy uses sum-norm
    inten_rel = inten / inten_max
    inten_sum = inten.sum()

    top5_pct = (np.sort(inten)[-5:].sum() / inten_sum) if inten_sum > 0 else np.nan

    mz_range = float(mz.max() - mz.min())
    ion_density = n / mz_range if mz_range > 0 else np.nan

    median_to_base = float(np.median(inten_rel))

    if spectral_entropy is None or pd.isna(spectral_entropy) or n < 2:
        # Recompute entropy if not provided or only one peak
        if inten_sum > 0:
            p = inten / inten_sum
            p = p[p > 0]
            S = float(-(p * np.log(p)).sum())
        else:
            S = np.nan
    else:
        S = float(spectral_entropy)

    if n >= 2 and not np.isnan(S):
        norm_ent = S / np.log(n)
    elif n < 2:
        # Single peak → S=0, normalized entropy undefined; treat as 0 (no diversity)
        norm_ent = 0.0
    else:
        norm_ent = np.nan

    return {
        'n_ions': int(n),
        'normalized_entropy': float(norm_ent) if not np.isnan(norm_ent) else np.nan,
        'top5_pct': float(top5_pct) if not np.isnan(top5_pct) else np.nan,
        'ion_density': float(ion_density) if not np.isnan(ion_density) else np.nan,
        'median_to_base': median_to_base,
    }


def main():
    print(f'Loading peaks cache: {PEAKS_CACHE}')
    with open(PEAKS_CACHE) as f:
        peaks_cache = json.load(f)
    print(f'  {len(peaks_cache):,} cached spectra')

    print(f'Loading feature table: {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows, {ft["wiki_id"].nunique():,} unique spectra')

    # Per-spectrum feature computation (one row per wiki_id)
    spec_entropy_lookup = (
        ft.drop_duplicates('wiki_id')[['wiki_id', 'spectral_entropy']]
        .set_index('wiki_id')['spectral_entropy'].to_dict()
    )

    print('Computing per-spectrum noise features...')
    rows = []
    n_missing = 0
    for wid in ft['wiki_id'].unique():
        peaks = peaks_cache.get(wid)
        if peaks is None:
            n_missing += 1
            rows.append({'wiki_id': wid, **{f: np.nan for f in NEW_FEATURES}})
            continue
        feats = compute_per_spectrum_features(peaks, spec_entropy_lookup.get(wid))
        rows.append({'wiki_id': wid, **feats})
    spec_feats = pd.DataFrame(rows)
    print(f'  computed for {len(spec_feats) - n_missing:,} spectra; '
          f'{n_missing:,} missing peaks (will be NaN)')

    # Sanity: print summary statistics + Oliver's threshold on labeled subset
    labeled_mask = ft['spectrum_label'].isin(['TP', 'FP'])
    labeled_spec = ft[labeled_mask].drop_duplicates('wiki_id')[['wiki_id']].merge(
        spec_feats, on='wiki_id', how='left')
    print('\nSummary on labeled (TP+FP) subset:')
    for c in NEW_FEATURES:
        v = labeled_spec[c].dropna()
        if len(v):
            print(f'  {c:22s}  median={v.median():.4f}   p05={v.quantile(0.05):.4f}'
                  f'   p95={v.quantile(0.95):.4f}')

    # Oliver's threshold check: 99.93% of annotated had S_norm < 0.987
    over_threshold = (labeled_spec['normalized_entropy'] >= 0.987).sum()
    n_total = labeled_spec['normalized_entropy'].notna().sum()
    print(f'\n  S_norm ≥ 0.987 on labeled (TP+FP): '
          f'{over_threshold} / {n_total} ({100*over_threshold/max(n_total,1):.2f}%)')
    print('  (Oliver: 99.93% of annotated <0.987, i.e. only 0.07% should be ≥0.987)')

    # Drop existing columns if any (re-runnable)
    for c in NEW_FEATURES:
        if c in ft.columns:
            ft = ft.drop(columns=[c])
            print(f'  dropping existing column "{c}" before re-merge')

    # Merge: every candidate row picks up its spectrum's noise features
    ft_out = ft.merge(spec_feats, on='wiki_id', how='left')
    assert len(ft_out) == len(ft), 'Row count changed after merge — duplicate wiki_id in spec_feats?'

    # Backup + overwrite
    if not os.path.exists(BACKUP):
        print(f'\nBacking up to {BACKUP}')
        shutil.copy(FEATURE_TABLE, BACKUP)
    print(f'Writing extended feature table to {FEATURE_TABLE}')
    ft_out.to_csv(FEATURE_TABLE, index=False)
    print(f'Done. Added 5 columns; {len(ft_out):,} rows × {len(ft_out.columns)} cols')


if __name__ == '__main__':
    main()
