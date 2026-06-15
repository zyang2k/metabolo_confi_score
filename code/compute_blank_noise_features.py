"""compute_blank_noise_features.py — Apply Oliver's noise features to fetched blanks.

Reads:
  data/blanks_query_peaks_cache_pos.json    (9,714 spectra)
  data/blanks_query_peaks_cache_neg.json    (1,702 spectra)
  data/blank_spectra_projected_scores.csv   (identity_score, projected_posterior)

For each blank spectrum, computes the same 5 features added to feature_table_v2.csv
in code/add_noise_features.py:
  - n_ions, normalized_entropy, top5_pct, ion_density, median_to_base

Reports:
  - Distribution stats per polarity
  - Snorm > 0.987 pass rate (Oliver's noise gate from CanMetCon deck slide 11)
  - Estimated usable corpus size after gate
  - Comparison with annotated subset

Output:
  data/blanks_noise_features.csv  — one row per blank wiki_id
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POS_CACHE = os.path.join(ROOT, 'data', 'blanks_query_peaks_cache_pos.json')
NEG_CACHE = os.path.join(ROOT, 'data', 'blanks_query_peaks_cache_neg.json')
PROJ_SCORES = os.path.join(ROOT, 'data', 'blank_spectra_projected_scores.csv')
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
OUT = os.path.join(ROOT, 'data', 'blanks_noise_features.csv')

NOISE_GATE = 0.987   # Oliver's Snorm threshold for Orbitrap


def compute_features(peaks: list[list[float]]) -> dict:
    """Same 5-feature spec as code/add_noise_features.py."""
    if not peaks:
        return {f: np.nan for f in
                ['n_ions', 'normalized_entropy', 'top5_pct',
                 'ion_density', 'median_to_base']}
    arr = np.asarray(peaks, dtype=np.float64)
    mz, inten = arr[:, 0], arr[:, 1]
    n = len(peaks)
    inten_max = inten.max()
    inten_sum = inten.sum()
    if inten_max <= 0 or inten_sum <= 0:
        return {f: np.nan for f in
                ['n_ions', 'normalized_entropy', 'top5_pct',
                 'ion_density', 'median_to_base']}

    # Compute spectral entropy (sum-normalized)
    p = inten / inten_sum
    p_pos = p[p > 0]
    S = float(-(p_pos * np.log(p_pos)).sum())
    norm_ent = (S / np.log(n)) if n >= 2 else 0.0

    top5_pct = float(np.sort(inten)[-5:].sum() / inten_sum)
    mz_range = float(mz.max() - mz.min())
    ion_density = (n / mz_range) if mz_range > 0 else np.nan
    median_to_base = float(np.median(inten / inten_max))

    return {
        'n_ions': int(n),
        'normalized_entropy': float(norm_ent),
        'top5_pct': top5_pct,
        'ion_density': ion_density,
        'median_to_base': median_to_base,
    }


def load_and_compute(path: str, polarity: str) -> pd.DataFrame:
    print(f'\nLoading {path}')
    with open(path) as f:
        cache = json.load(f)
    print(f'  {len(cache):,} spectra in cache')
    rows = []
    n_empty = 0
    for wid, peaks in cache.items():
        if peaks is None or len(peaks) == 0:
            n_empty += 1
            rows.append({'wiki_id': wid, 'polarity': polarity,
                         **{k: np.nan for k in ['n_ions', 'normalized_entropy',
                                                'top5_pct', 'ion_density',
                                                'median_to_base']}})
            continue
        rows.append({'wiki_id': wid, 'polarity': polarity, **compute_features(peaks)})
    print(f'  computed for {len(rows) - n_empty:,} spectra; {n_empty:,} empty')
    return pd.DataFrame(rows)


def main():
    pos = load_and_compute(POS_CACHE, 'pos')
    neg = load_and_compute(NEG_CACHE, 'neg')
    blanks = pd.concat([pos, neg], ignore_index=True)

    # Merge with projected_posterior + identity_score
    print(f'\nMerging projected_posterior from {PROJ_SCORES}')
    proj = pd.read_csv(PROJ_SCORES)
    blanks = blanks.merge(proj, on='wiki_id', how='left')
    print(f'  rows after merge: {len(blanks):,}; with projected_posterior: '
          f'{blanks["projected_posterior"].notna().sum():,}')

    blanks.to_csv(OUT, index=False)
    print(f'\nWrote {OUT}: {len(blanks):,} rows')

    # ── Distribution summary ──
    def summary(df: pd.DataFrame, label: str) -> None:
        print(f'\n=== {label} (n={len(df):,}) ===')
        for c in ['n_ions', 'normalized_entropy', 'top5_pct',
                  'ion_density', 'median_to_base']:
            v = df[c].dropna()
            if len(v):
                print(f'  {c:22s}  median={v.median():.4f}   p05={v.quantile(0.05):.4f}'
                      f'   p95={v.quantile(0.95):.4f}')

    summary(blanks[blanks['polarity'] == 'pos'], 'POS blanks')
    summary(blanks[blanks['polarity'] == 'neg'], 'NEG blanks')

    # ── Oliver's Snorm 0.987 gate ──
    print(f'\n=== Oliver Snorm > {NOISE_GATE} noise gate ===')
    for pol in ['pos', 'neg', 'all']:
        sub = blanks if pol == 'all' else blanks[blanks['polarity'] == pol]
        sn = sub['normalized_entropy'].dropna()
        n_total = len(sn)
        n_noise = (sn > NOISE_GATE).sum()
        print(f'  {pol:5s}: {n_noise:,} / {n_total:,} '
              f'noise ({100*n_noise/max(n_total,1):.1f}%) | '
              f'usable: {n_total - n_noise:,} ({100*(n_total-n_noise)/max(n_total,1):.1f}%)')

    # Compare to annotated baseline (from feature_table_v2 which already has noise features)
    print(f'\n=== Annotated baseline (from feature_table_v2.csv) ===')
    if os.path.exists(FEATURE_TABLE):
        ft = pd.read_csv(FEATURE_TABLE, low_memory=False,
                         usecols=['wiki_id', 'spectrum_label', 'normalized_entropy', 'polarity'])
        ann = ft[ft['spectrum_label'].isin(['TP', 'FP'])].drop_duplicates('wiki_id')
        sn_ann = ann['normalized_entropy'].dropna()
        n_noise_ann = (sn_ann > NOISE_GATE).sum()
        print(f'  annotated Snorm > {NOISE_GATE}: {n_noise_ann}/{len(sn_ann):,} '
              f'({100*n_noise_ann/max(len(sn_ann),1):.2f}%)')
        print(f'  annotated median Snorm: {sn_ann.median():.4f}')
        print(f'  blank     median Snorm: {blanks["normalized_entropy"].median():.4f}')

    # ── Stratified deliverable view ──
    print(f'\n=== Stratified usability (HILIC Orbi blanks) ===')
    print(f'  Total fetched              : {len(blanks):,}')
    pass_gate = blanks[blanks['normalized_entropy'] <= NOISE_GATE]
    print(f'  Pass Snorm gate            : {len(pass_gate):,}'
          f' ({100*len(pass_gate)/len(blanks):.1f}%)')

    # Stratify by projected_posterior bucket
    pp = blanks['projected_posterior']
    for thr_lo, thr_hi, label in [(0.0, 0.3, 'low conf <0.3'),
                                   (0.3, 0.5, 'mid 0.3-0.5'),
                                   (0.5, 0.7, 'high 0.5-0.7'),
                                   (0.7, 1.01, 'top ≥0.7')]:
        m = (pp >= thr_lo) & (pp < thr_hi)
        n_seg = m.sum()
        n_seg_pass = (m & (blanks['normalized_entropy'] <= NOISE_GATE)).sum()
        print(f'    {label:20s}: {n_seg:,} ({n_seg_pass:,} pass gate)')


if __name__ == '__main__':
    main()
