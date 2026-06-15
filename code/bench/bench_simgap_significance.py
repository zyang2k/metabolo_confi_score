"""bench_simgap_significance.py — ion-counting significance of the similarity gap.

Builds per-spectrum gap-significance columns (see simgap_significance.py) and evaluates
them as GBM features against the production baseline, plus the diagnostics that justify
the metric independent of any AUC bump:

  1. API-vs-local entropy-similarity agreement (tolerance pick).
  2. Methylcytidine-style near-tie reproduction: a small gap with low z and ~coin-flip
     direction stability — the bug Oliver flagged.
  3. N_eff sensitivity (the one free parameter): does the tied/not-tied verdict survive
     a 10x swing in the assumed counting budget?
  4. GBM ablation: baseline vs +sim_gap_z (+ gap_sigma, p_top1_stable), full slice and
     restricted to the confusable band (the only place the feature can act).
  5. Deterministic-gate preview: how many top-1 calls are MS²-tied (z < Z_TIED) and would
     defer to RT — the lever that fixes the methylcytidine call.

Side table cached at data/sim_gap_significance.csv (slow to build; ~minutes). Pass
--rebuild to regenerate. Diagnostic only — writes no model, changes no production feature.
"""
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

import bench_harness as bh
from bench_harness import BenchSpec, run_bench, build_top1, make_groups, run_arm, summarize
from simgap_significance import (
    build_significance_table, gap_significance_for_spectrum, local_entropy_sim,
    DEFAULT_N_EFF, DEFAULT_B, ENTROPY_TOL_DA,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIDE = os.path.join(ROOT, 'data', 'sim_gap_significance.csv')
QC = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')
LC = os.path.join(ROOT, 'data', 'library_peaks_cache.json')
Z_TIED = 2.0          # gap < 2σ of its own counting noise → MS²-tied → defer to RT
EXTRA = ['sim_gap_z', 'gap_sigma', 'p_top1_stable']


def load_caches():
    return json.load(open(QC)), json.load(open(LC))


def diag_tolerance(ft, qc, lc):
    print('\n=== (1) local-vs-API entropy_similarity agreement ===')
    m = ft[ft['library_wiki_id'].astype(str).isin(lc.keys())
           & ft['wiki_id'].isin(qc.keys())
           & pd.to_numeric(ft['entropy_similarity'], errors='coerce').notna()]
    samp = m.sample(min(1500, len(m)), random_state=0)
    for tol in [0.01, 0.02, 0.05]:
        loc, api = [], []
        for _, r in samp.iterrows():
            s = local_entropy_sim(qc[r['wiki_id']], lc[str(r['library_wiki_id'])], tol)
            if np.isfinite(s):
                loc.append(s); api.append(float(r['entropy_similarity']))
        loc, api = np.array(loc), np.array(api)
        print(f'  tol={tol:.2f}: corr={np.corrcoef(loc, api)[0, 1]:.4f}  '
              f'median|Δ|={np.median(np.abs(loc - api)):.4f}  (n={len(loc)})')


def _competitors(ft, wid, lc):
    g = ft[(ft['wiki_id'] == wid) & ft['hit_ik14'].fillna('').ne('')]
    return [{'library_wiki_id': l, 'lib_peaks': lc.get(str(l))}
            for l in g['library_wiki_id'].tolist()]


def diag_methylcytidine(ft, qc, lc):
    print('\n=== (2) methylcytidine-style near-ties (low z, ~coin-flip direction) ===')
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    sub = ft[ft['wiki_id'].isin(set(lab['wiki_id'])) & ft['hit_ik14'].fillna('').ne('')].copy()
    sub = sub[sub['library_wiki_id'].astype(str).isin(lc.keys()) & sub['wiki_id'].isin(qc.keys())]
    nc = sub.groupby('wiki_id').size()
    cand = nc[nc >= 2].index
    rows = []
    for wid in cand:
        r = gap_significance_for_spectrum(qc[wid], _competitors(ft, wid, lc),
                                          B=120, seed=abs(hash(wid)) % 2**32)
        if np.isfinite(r['sim_gap_z']):
            rows.append((wid, r['sim_top1_local'], r['sim_top2_local'],
                         r['gap_point'], r['gap_sigma'], r['sim_gap_z'], r['p_top1_stable']))
    d = pd.DataFrame(rows, columns=['wiki_id', 'sim1', 'sim2', 'gap', 'sigma', 'z', 'p_stable'])
    ties = d[d['z'] < Z_TIED].sort_values('gap')
    print(f'  {len(d)} confusable spectra scored; {len(ties)} are MS²-tied (z<{Z_TIED}, '
          f'{100*len(ties)/max(len(d),1):.1f}%)')
    print('  Sample near-ties (small gap, low z, ~0.5 direction stability):')
    print(ties.head(8).to_string(index=False, float_format=lambda x: f'{x:.4f}'))
    print('  Sample clear calls (large z):')
    print(d.sort_values('z', ascending=False).head(4).to_string(
        index=False, float_format=lambda x: f'{x:.4f}'))


def diag_neff(ft, qc, lc):
    print('\n=== (3) N_eff sensitivity (tied verdict should survive a 10x swing) ===')
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    sub = ft[ft['wiki_id'].isin(set(lab['wiki_id'])) & ft['hit_ik14'].fillna('').ne('')]
    sub = sub[sub['library_wiki_id'].astype(str).isin(lc.keys()) & sub['wiki_id'].isin(qc.keys())]
    nc = sub.groupby('wiki_id').size()
    wids = nc[nc >= 2].index.to_series().sample(min(300, (nc >= 2).sum()), random_state=3).tolist()
    res = {}
    for neff in [300, 1000, 3000]:
        zs = []
        for wid in wids:
            r = gap_significance_for_spectrum(qc[wid], _competitors(ft, wid, lc),
                                              n_eff=neff, B=100, seed=abs(hash(wid)) % 2**32)
            zs.append(r['sim_gap_z'])
        res[neff] = np.array(zs)
    base = res[1000]
    tied_1000 = base < Z_TIED
    print(f'  N_eff=1000: {tied_1000.sum()}/{len(base)} tied (z<{Z_TIED})')
    for neff in [300, 3000]:
        tied = res[neff] < Z_TIED
        agree = np.mean(tied == tied_1000)
        print(f'  N_eff={neff:5d}: {tied.sum():3d} tied;  verdict agreement vs N_eff=1000 = {100*agree:.1f}%')


def build_side_table(ft, qc, lc, n_eff, B):
    lab = ft[ft['spectrum_label'].isin(['TP', 'FP'])]
    wiki_ids = lab['wiki_id'].unique()
    print(f'\nBuilding significance table for {len(wiki_ids):,} labeled spectra '
          f'(N_eff={n_eff}, B={B})...')
    tab = build_significance_table(ft, qc, lc, wiki_ids=wiki_ids, n_eff=n_eff, B=B)
    tab.to_csv(SIDE, index=False)
    print(f'  wrote {SIDE} ({len(tab):,} rows, '
          f'{tab["sim_gap_z"].notna().mean()*100:.1f}% with finite z)')
    return tab


def bench_gbm(ft):
    print('\n=== (4) GBM ablation: baseline vs +significance features ===')
    spec = BenchSpec(name='simgap_significance', extra_numeric=EXTRA,
                     extra_table=SIDE, merge_on=('wiki_id',))
    out = run_bench(spec)

    # confusable-band-restricted AUC: where the feature can actually act.
    top1 = out['top1']
    band = top1['sim_gap_z'].notna() & (top1['n_candidates'].fillna(1) > 1)
    print(f'\n=== Confusable-band AUC (n={band.sum():,}, the slice the feature touches) ===')
    from sklearn.metrics import roc_auc_score
    y = top1['hit_label'].values
    for arm in [out['baseline'], out['extended']]:
        sub_auc = roc_auc_score(y[band.values], arm.oof[band.values])
        print(f'  {arm.name:22s} band AUC = {sub_auc:.4f}')

    # (5) deterministic-gate preview
    z = top1['sim_gap_z']
    tied = (z < Z_TIED) & z.notna()
    print(f'\n=== (5) Gate preview: MS²-tied top-1 calls (z<{Z_TIED}) ===')
    print(f'  {tied.sum():,} of {z.notna().sum():,} scored top-1 calls are MS²-tied '
          f'({100*tied.sum()/max(z.notna().sum(),1):.1f}%)')
    print('  → on these, MS² cannot separate the candidates; the gate would defer to RT.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rebuild', action='store_true', help='regenerate the side table')
    ap.add_argument('--n_eff', type=int, default=DEFAULT_N_EFF)
    ap.add_argument('--B', type=int, default=DEFAULT_B)
    ap.add_argument('--skip-diag', action='store_true')
    args = ap.parse_args()

    ft = pd.read_csv(bh.FEATURE_TABLE, low_memory=False)
    qc, lc = load_caches()

    if not args.skip_diag:
        diag_tolerance(ft, qc, lc)
        diag_methylcytidine(ft, qc, lc)
        diag_neff(ft, qc, lc)

    if args.rebuild or not os.path.exists(SIDE):
        build_side_table(ft, qc, lc, args.n_eff, args.B)

    bench_gbm(ft)


if __name__ == '__main__':
    main()
