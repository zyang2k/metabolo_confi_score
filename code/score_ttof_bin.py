"""score_ttof_bin.py — Reproduce the cross-platform TTOF confidence for one bin.

Same path as code/bench_ttof_golden.py: build the TTOF feature table from library_hits +
spectrum metadata, train the Orbitrap-trained GBM (BASE_NUMERIC + struct_logit), and predict
on the requested bin's top-1-by-entropy_similarity candidate. Unlike the bench, this scores a
SINGLE bin regardless of whether it carries a golden/yy_ label, and fills the target bin's MS2
cosines on demand so its feature row is complete.

Caveats (memory: project_ttof_cross_platform_20260526):
  * TTOF is ADVISORY tier — the Orbitrap-trained model transfers at AUC ~0.85; the number is a
    raw GBM probability, NOT isotonic-calibrated and NOT comparable 1:1 to the Orbitrap
    deliverable's confidence_pct.
  * Mirrors the bench's struct_logit feature set, which PREDATES the gate-2 relational features
    now in score_gbm_v2.py. It is the repo's canonical TTOF scorer, not current Orbitrap prod.
  * The missing-RT haircut (RT_MISSING_CONFIDENCE_HAIRCUT) lives only in the Orbitrap deliverable
    step. The TTOF path drops candidates without a predicted RT upstream, so it never applies
    here; we report whether the target bin's top-1 candidate has an RT prediction.

Usage: python code/score_ttof_bin.py <wiki_id> [pos|neg]
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import xgboost as xgb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bench_harness import BASE_NUMERIC, CATEGORICAL, FEATURE_TABLE, build_top1, prep, train
from build_features_v2 import compute_ms2_scores
import bench_ttof_golden as btg


def main():
    wiki_id = sys.argv[1] if len(sys.argv) > 1 else 'aJWMOO8/7191'
    polarity = (sys.argv[2] if len(sys.argv) > 2 else 'pos').lower()
    pol_val = 1 if polarity == 'pos' else 0
    hits_p = btg.TTOF_POS_HITS if pol_val else btg.TTOF_NEG_HITS
    spec_p = btg.TTOF_POS_SPEC if pol_val else btg.TTOF_NEG_SPEC
    peaks_p = btg.TTOF_POS_PEAKS if pol_val else btg.TTOF_NEG_PEAKS

    print(f'Target: {wiki_id}  (TTOF {polarity})')
    adduct_tax = btg.load_adduct_taxonomy()
    print('Loading library peaks cache (~335 MB) ...')
    with open(btg.LIB_PEAKS) as f:
        lib_peaks = json.load(f)

    ttof = btg.build_ttof_table(hits_p, spec_p, peaks_p, lib_peaks, pol_val, adduct_tax)

    # ---- the target bin's top-1 candidate by entropy_similarity ----
    rows = ttof[ttof.wiki_id == wiki_id]
    if rows.empty:
        # Bin may have been dropped (no candidate with delta_predicted_rt, or empty hit_ik14).
        print(f'\n{wiki_id} not present in the built TTOF table — likely all candidates lack a '
              f'predicted RT (dropped upstream) or have empty hit_ik14. No score.')
        return
    target = rows.loc[rows.entropy_similarity.idxmax()].copy()

    # Fill MS2 cosines for the target if not already (bench only computes them on labeled bins).
    if not np.isfinite(target.get('forward_cosine', np.nan)):
        with open(peaks_p) as f:
            q_peaks = json.load(f)
        q = q_peaks.get(str(wiki_id))
        lid = target.get('library_wiki_id')
        l = lib_peaks.get(str(lid)) if isinstance(lid, str) else None
        if q is not None and l is not None:
            fwd, rev, cc, ci = compute_ms2_scores(q, l)
            target['forward_cosine'], target['reverse_cosine'] = fwd, rev
            target['cov_count'], target['cov_int'] = cc, ci
            print(f'  filled target MS2 cosines: fwd={fwd:.3f} rev={rev:.3f} cov_count={cc} cov_int={ci:.3f}')
        else:
            print('  target MS2 peaks unavailable (query or library cache miss) — cosines left NaN')

    # ---- train the Orbitrap-trained model exactly as the bench does ----
    print('\nTraining Orbitrap GBM (BASE_NUMERIC + struct_logit) ...')
    orb = pd.read_csv(FEATURE_TABLE, low_memory=False)
    side = pd.read_csv(btg.MOLREX_FEATURES)
    side['wiki_id'] = side.wiki_id.astype(str)
    orb['wiki_id'] = orb.wiki_id.astype(str)
    orb = orb.merge(side[['wiki_id', 'hit_ik14', 'struct_logit']], on=['wiki_id', 'hit_ik14'], how='left')
    orb_top1 = build_top1(orb).reset_index(drop=True)
    X_train = prep(orb_top1, ['struct_logit'])
    y_train = orb_top1['hit_label'].values
    n_num, n_cat = len(BASE_NUMERIC) + 1, len(CATEGORICAL)
    model = train(X_train, y_train, n_num, n_cat)

    # ---- predict on the single target row ----
    tgt_df = pd.DataFrame([target])
    X_tgt = prep(tgt_df, ['struct_logit'])
    dt = xgb.DMatrix(X_tgt, enable_categorical=True, feature_types=['q'] * n_num + ['c'] * n_cat)
    pred = float(model.predict(dt)[0])

    rt_missing = not np.isfinite(pd.to_numeric(pd.Series([target.get('signed_delta_rt')]), errors='coerce')[0])
    print('\n' + '=' * 64)
    print(f'TTOF confidence for {wiki_id}')
    print('=' * 64)
    print(f'  annotation (spec)   : {target.get("name")}   adduct (hit) {target.get("adduct")}')
    print(f'  hit_adduct_cat      : {target.get("hit_adduct_cat")}   db {target.get("db")}   label {target.get("spectrum_label")}')
    print(f'  entropy_similarity  : {target.get("entropy_similarity"):.4f}')
    print(f'  sim_gap (vs 2nd)    : {target.get("sim_gap"):.4f}')
    print(f'  signed_delta_rt     : {target.get("signed_delta_rt")}   (RT prediction {"MISSING" if rt_missing else "present"})')
    print(f'  delta_mda           : {target.get("delta_mda"):.3f}')
    print(f'  forward/reverse cos : {target.get("forward_cosine")}/{target.get("reverse_cosine")}')
    print(f'  struct_logit        : {target.get("struct_logit")}')
    print(f'\n  RAW GBM confidence  : {pred:.4f}  ({pred*100:.1f}%)   [uncalibrated, advisory tier]')
    if rt_missing:
        print(f'  (RT missing → Orbitrap deliverable would apply the {0.15:.0%} haircut → '
              f'{pred*0.85:.4f}; not applied in the TTOF path)')


if __name__ == '__main__':
    main()
