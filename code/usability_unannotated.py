"""usability_unannotated.py — Sanity-check unannotated (blank) spectra.

Question: are the 11,416 in-house unannotated spectra distributionally similar
enough to the 5,638 annotated spectra to be useful as self-supervised
pretraining data, or are they degenerate (all noise / all at one m/z / etc.)?

Tests use only what's already on disk (no fetch). If blanks pass these checks,
the next step is fetching their raw peaks to apply Oliver's full noise filter
(Snorm > 0.987 gate + 5-feature SVM). If blanks fail, save the fetch cost.

Tests:
  T1. Spectral entropy distribution overlap (annotated vs blanks).
      Blanks should NOT be all entropy=0 (single-peak / dead spectra).
  T2. Precursor m/z range coverage (annotated vs blanks).
      Blanks shouldn't pile up at one m/z (system contaminant).
  T3. Retention time coverage.
      Blanks shouldn't all elute at t=0 (void volume / matrix).
  T4. Hits-per-spectrum from blanks_hits files.
      If blanks have very few candidate hits each, the library search itself
      is failing — bad sign for self-supervised value.
  T5. identity_score and projected_posterior distributions.
      How "blank" are these blanks? identity_score < 0.7 by definition, but
      projected_posterior tells us how the calibrated GBM scores them.
  T6. Estimated effective corpus size after Oliver-style threshold.
      We can't compute Snorm without n_ions, but we can estimate using a
      surrogate: spectra with entropy > 0.3 (low entropy proxy for
      "non-trivial information content").
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEG_CUR = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_negESI_curated_042126.csv')
POS_CUR = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_posESI_curated_042126.csv')
NEG_BL_HITS = os.path.join(ROOT, 'data', 'neg_blanks_hits.csv')
POS_BL_HITS = os.path.join(ROOT, 'data', 'pos_blanks_hits.csv')
PROJ_SCORES = os.path.join(ROOT, 'data', 'blank_spectra_projected_scores.csv')


def load_curated_split() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (annotated, blank) DataFrames pooling neg + pos curated CSVs."""
    rows_anno = []
    rows_blank = []
    for path, polarity in [(NEG_CUR, 'neg'), (POS_CUR, 'pos')]:
        if not os.path.exists(path):
            print(f'  skipping (missing): {path}')
            continue
        df = pd.read_csv(path, low_memory=False)
        df['polarity'] = polarity
        is_blank = (
            (df['identity_score'] < 0.7) &
            (df['name'].isna() | (df['name'].astype(str).str.strip() == ''))
        )
        rows_anno.append(df[~is_blank])
        rows_blank.append(df[is_blank])
        print(f'  {polarity}: total={len(df):,}  annotated={(~is_blank).sum():,}  '
              f'blank={is_blank.sum():,}')
    anno = pd.concat(rows_anno, ignore_index=True)
    blank = pd.concat(rows_blank, ignore_index=True)
    return anno, blank


def quantile_table(s: pd.Series, name: str) -> pd.Series:
    return pd.Series({
        f'{name}_n':     int(s.notna().sum()),
        f'{name}_p05':   s.quantile(0.05),
        f'{name}_p25':   s.quantile(0.25),
        f'{name}_p50':   s.quantile(0.50),
        f'{name}_mean':  s.mean(),
        f'{name}_p75':   s.quantile(0.75),
        f'{name}_p95':   s.quantile(0.95),
    })


def t1_entropy(anno: pd.DataFrame, blank: pd.DataFrame) -> None:
    print('\n=== T1: Spectral entropy distribution ===')
    a = pd.to_numeric(anno['entropy'], errors='coerce').dropna()
    b = pd.to_numeric(blank['entropy'], errors='coerce').dropna()
    print(f'  annotated  n={len(a):,}  median={a.median():.3f}  '
          f'mean={a.mean():.3f}  p05={a.quantile(0.05):.3f}  p95={a.quantile(0.95):.3f}')
    print(f'  blank      n={len(b):,}  median={b.median():.3f}  '
          f'mean={b.mean():.3f}  p05={b.quantile(0.05):.3f}  p95={b.quantile(0.95):.3f}')
    n_zero_a = (a == 0).sum()
    n_zero_b = (b == 0).sum()
    print(f'  entropy=0 (single-peak / dead): annotated={n_zero_a} ({100*n_zero_a/max(len(a),1):.1f}%), '
          f'blank={n_zero_b} ({100*n_zero_b/max(len(b),1):.1f}%)')
    n_tiny_b = (b < 0.3).sum()
    print(f'  blank entropy < 0.3 (likely noise): {n_tiny_b}/{len(b)} ({100*n_tiny_b/max(len(b),1):.1f}%)')


def t2_precursor(anno: pd.DataFrame, blank: pd.DataFrame) -> None:
    print('\n=== T2: Precursor m/z coverage ===')
    a = pd.to_numeric(anno['precursor_mz'], errors='coerce').dropna()
    b = pd.to_numeric(blank['precursor_mz'], errors='coerce').dropna()
    print(f'  annotated  range [{a.min():.1f}, {a.max():.1f}]  '
          f'median={a.median():.1f}  IQR=[{a.quantile(0.25):.1f}, {a.quantile(0.75):.1f}]')
    print(f'  blank      range [{b.min():.1f}, {b.max():.1f}]  '
          f'median={b.median():.1f}  IQR=[{b.quantile(0.25):.1f}, {b.quantile(0.75):.1f}]')
    # Look for pile-up at single m/z (contamination tell-tale)
    bins_b = pd.cut(b, bins=20)
    top_bin_frac = bins_b.value_counts(normalize=True).iloc[0]
    print(f'  blank pile-up: top 5%-mass-bin holds {100*top_bin_frac:.1f}% of blanks')


def t3_rt(anno: pd.DataFrame, blank: pd.DataFrame) -> None:
    print('\n=== T3: Retention time coverage ===')
    a = pd.to_numeric(anno['rt'], errors='coerce').dropna()
    b = pd.to_numeric(blank['rt'], errors='coerce').dropna()
    print(f'  annotated  range [{a.min():.2f}, {a.max():.2f}]  '
          f'median={a.median():.2f}  IQR=[{a.quantile(0.25):.2f}, {a.quantile(0.75):.2f}]')
    print(f'  blank      range [{b.min():.2f}, {b.max():.2f}]  '
          f'median={b.median():.2f}  IQR=[{b.quantile(0.25):.2f}, {b.quantile(0.75):.2f}]')
    # Void / matrix peak indicator: large fraction at small RT
    n_void_a = (a < 0.5).sum()
    n_void_b = (b < 0.5).sum()
    print(f'  rt < 0.5 min (void-volume): '
          f'annotated={n_void_a} ({100*n_void_a/max(len(a),1):.1f}%), '
          f'blank={n_void_b} ({100*n_void_b/max(len(b),1):.1f}%)')


def t4_hits(blank: pd.DataFrame) -> None:
    print('\n=== T4: Hits per blank spectrum (library-search productivity) ===')
    blank_wids = set(blank['wiki_id'].astype(str))
    for path, label in [(NEG_BL_HITS, 'neg'), (POS_BL_HITS, 'pos')]:
        if not os.path.exists(path):
            continue
        hits = pd.read_csv(path, usecols=['wiki_id'], low_memory=False)
        hits['wiki_id'] = hits['wiki_id'].astype(str)
        cnt = hits.groupby('wiki_id').size()
        # Filter to blank wiki_ids that match this polarity's curated set
        cnt = cnt[cnt.index.isin(blank_wids)]
        if len(cnt) == 0:
            print(f'  {label}: 0 blank wiki_ids overlap with hit file')
            continue
        zero_hits = len(blank_wids) - len(cnt)
        print(f'  {label}: {len(cnt):,} blanks have hits   |   '
              f'mean={cnt.mean():.1f}  median={int(cnt.median())}  '
              f'p95={int(cnt.quantile(0.95))}  max={int(cnt.max())}   |   '
              f'~{zero_hits:,} blanks may have 0 hits')


def t5_projected_score(blank: pd.DataFrame) -> None:
    print('\n=== T5: Identity score & projected posterior on blanks ===')
    a = pd.to_numeric(blank['identity_score'], errors='coerce').dropna()
    print(f'  identity_score (in curated CSV, def < 0.7 for blanks)  '
          f'n={len(a):,}  median={a.median():.3f}  '
          f'p05={a.quantile(0.05):.3f}  p95={a.quantile(0.95):.3f}')
    if not os.path.exists(PROJ_SCORES):
        print(f'  {PROJ_SCORES} not found, skipping projected_posterior')
        return
    proj = pd.read_csv(PROJ_SCORES)
    p = pd.to_numeric(proj['projected_posterior'], errors='coerce').dropna()
    print(f'  projected_posterior (GBM on blanks)  n={len(p):,}  '
          f'median={p.median():.3f}  p05={p.quantile(0.05):.3f}  '
          f'p95={p.quantile(0.95):.3f}')
    for thr in [0.3, 0.5, 0.7, 0.9]:
        n_above = (p >= thr).sum()
        print(f'    projected_posterior ≥ {thr}: {n_above} ({100*n_above/max(len(p),1):.2f}%)')


def t6_corpus_estimate(anno: pd.DataFrame, blank: pd.DataFrame) -> None:
    print('\n=== T6: Estimated usable corpus size for self-supervised pretraining ===')
    b_ent = pd.to_numeric(blank['entropy'], errors='coerce')

    # Surrogate filter without n_ions: entropy > 0.3
    surrogate_passing = (b_ent > 0.3).sum()
    print(f'  surrogate filter (entropy > 0.3): {surrogate_passing:,} of {len(blank):,} '
          f'blanks ({100*surrogate_passing/max(len(blank),1):.1f}%)')

    # Compare to known annotated entropy distribution
    a_ent = pd.to_numeric(anno['entropy'], errors='coerce')
    a_passing = (a_ent > 0.3).sum()
    print(f'  same surrogate on annotated: {a_passing:,}/{len(anno):,} '
          f'({100*a_passing/max(len(anno),1):.1f}%)')

    # Reference points
    print('\n  reference scales for self-supervised LC-MS:')
    print('    DreaMS              : 24,000,000 spectra')
    print('    MS2DeepScore        :    109,734 spectra (15,062 IK14)')
    print('    Spec2Vec            :     95,320 spectra')
    print('    MIST (NPLIB1)       :      8,030 spectra (7,131 IK14)')
    print('    --- our scale ---')
    print(f'    OUR labeled         :      5,638 spectra')
    print(f'    OUR annotated total :      {a_passing:,} (entropy>0.3 filter)')
    print(f'    OUR blanks usable   :      {surrogate_passing:,} (entropy>0.3 filter)')
    print(f'    OUR combined potl.  :      {a_passing + surrogate_passing:,}')


def main():
    print('Loading curated CSVs...')
    anno, blank = load_curated_split()
    print(f'\nTotal: annotated={len(anno):,}  blank={len(blank):,}')

    t1_entropy(anno, blank)
    t2_precursor(anno, blank)
    t3_rt(anno, blank)
    t4_hits(blank)
    t5_projected_score(blank)
    t6_corpus_estimate(anno, blank)

    print('\n=== Verdict logic ===')
    print('  Pass criteria for "worth fetching peaks":')
    print('    - blank entropy median should be > 0.5 (not all dead spectra)')
    print('    - blank precursor m/z range should overlap annotated >= 80%')
    print('    - blank rt distribution should not pile up at void volume')
    print('    - >= 50% of blanks should have ≥ 1 library hit')
    print('    - corpus size after surrogate filter should be > 3,000')
    print('  See above for actual numbers — judge by inspection.')


if __name__ == '__main__':
    main()
