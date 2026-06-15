"""score_oliver_ttof_v2.py — Action 2: honest before/after on Oliver's actual TTOF bins.

Scores his two flagged TTOF-pos bins through OUR pipeline, with the fixes computed WITHIN TTOF
(not against the Orbitrap graph) so we can see whether they fire cross-platform:
  - methylcytidine isomer tie:  aJWMOO8/7191  (2'-O-methylcytidine, esim 1.0)  -> sim_gap fix
  - aminoadipic ISF:            search 2-aminoadipic [M+H-H2O]+ bins           -> RT-confirmation

Base model = Orbitrap-trained (the cross-platform weak spot — report honestly).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.isotonic import IsotonicRegression
import bench_ttof_golden as G
from bench_harness import BASE_NUMERIC, CATEGORICAL, XGB_PARAMS, NUM_BOOST_ROUND

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RT_WIN, ESIM_MIN, SIMGAP_EPS = 10.0, 0.50, 0.05


def add_within_ttof_fixes(ttof):
    """Recompute compound_has_ok_adduct WITHIN TTOF (within-bin OR same-RT co-eluting clean parent),
    recompute hit_isf_no_ok, and add a soft-thresholded sim_gap — the two Oliver fixes, computed
    on TTOF itself rather than the frozen Orbitrap graph."""
    t = ttof.copy()
    t['entropy_similarity'] = pd.to_numeric(t['entropy_similarity'], errors='coerce')
    t['measured_rt'] = pd.to_numeric(t['measured_rt'], errors='coerce')
    # FIX: bench_ttof_golden leaves hit_adduct_cat='unknown' on TTOF (taxonomy keys are
    # bracket-free, TTOF strings have brackets). Reclassify with the working Orbitrap
    # classifier so the ISF/ok-adduct logic actually applies on TTOF.
    import build_features_v2 as B
    lookup = B._load_adduct_taxonomy()
    t['hit_adduct_cat'] = t['adduct'].apply(lambda a: B.classify_adduct(a, lookup))
    t['hit_is_isf'] = (t['hit_adduct_cat'] == 'isf').astype(int)
    print('  reclassified TTOF adducts:', t['hit_adduct_cat'].value_counts().to_dict())
    t['_is_ok'] = (t['hit_adduct_cat'] == 'ok').astype(int)
    # within-bin compound_has_ok_adduct (the deployed/baseline definition)
    wb = t.groupby(['wiki_id', 'hit_ik14'])['_is_ok'].max().rename('choa_within').reset_index()
    t = t.merge(wb, on=['wiki_id', 'hit_ik14'], how='left')
    # RT-confirmation within TTOF: same ik14, clean ok adduct, esim>=thr, co-eluting in another bin
    conf = t[(t['_is_ok'] == 1) & (t['entropy_similarity'] >= ESIM_MIN) & t['measured_rt'].notna()]
    cb = conf.groupby(['hit_ik14', 'wiki_id'])['measured_rt'].first().reset_index()
    by_ik = {ik: g[['wiki_id', 'measured_rt']].values for ik, g in cb.groupby('hit_ik14')}
    rtmap = t.drop_duplicates('wiki_id').set_index('wiki_id')['measured_rt']
    pairs = t[['wiki_id', 'hit_ik14']].drop_duplicates()
    def confd(w, ik):
        arr = by_ik.get(ik); rt = rtmap.get(w)
        if arr is None or rt is None or not np.isfinite(rt): return 0
        return int(any(w2 != w and abs(rt - rt2) <= RT_WIN for w2, rt2 in arr))
    pairs['rtc'] = [confd(w, ik) for w, ik in zip(pairs.wiki_id, pairs.hit_ik14)]
    t = t.merge(pairs, on=['wiki_id', 'hit_ik14'], how='left')
    t['choa_fixed'] = ((t['choa_within'] == 1) | (t['rtc'] == 1)).astype(int)
    # soft-threshold sim_gap
    g = pd.to_numeric(t['sim_gap'], errors='coerce')
    t['simgap_fixed'] = np.sign(g) * np.maximum(0, np.abs(g) - SIMGAP_EPS)
    return t


def prep(df, cm=None):
    X = pd.DataFrame(index=df.index); fit = cm is None; cm = cm or {}
    for c in BASE_NUMERIC: X[c] = pd.to_numeric(df[c], errors='coerce') if c in df else np.nan
    for c in CATEGORICAL:
        s = df[c].astype(str).fillna('missing') if c in df else pd.Series('missing', index=df.index)
        if fit: cm[c] = {v: i for i, v in enumerate(sorted(s.unique()))}
        mx = max(cm[c].values())+1 if cm[c] else 0
        X[c] = s.map(lambda v: cm[c].get(v, mx)).astype('int32')
    return X[BASE_NUMERIC+CATEGORICAL], cm


def main():
    tax = G.load_adduct_taxonomy()
    print('Loading library peaks cache ...')
    lib = json.load(open(G.LIB_PEAKS))
    ttof = G.build_ttof_table(G.TTOF_POS_HITS, G.TTOF_POS_SPEC, G.TTOF_POS_PEAKS, lib, 1, tax)
    ttof = add_within_ttof_fixes(ttof)
    print(f'TTOF pos: {len(ttof):,} rows, {ttof.wiki_id.nunique():,} bins')

    # train Orbitrap GBM (production BASE features) + trustworthy isotonic
    orb = pd.read_csv(G.FEATURE_TABLE, low_memory=False)
    lab = orb[orb.spectrum_label.isin(['TP', 'FP']) & orb.anno_ik14.fillna('').ne('')]
    top1 = orb.loc[lab.groupby('wiki_id').entropy_similarity.idxmax()].reset_index(drop=True)
    Xtr, cm = prep(top1); y = top1.hit_label.values
    ftypes = ['q']*len(BASE_NUMERIC)+['c']*len(CATEGORICAL)
    m = xgb.train(XGB_PARAMS, xgb.DMatrix(Xtr, label=y, enable_categorical=True, feature_types=ftypes), NUM_BOOST_ROUND)
    raw_tr = m.predict(xgb.DMatrix(Xtr, enable_categorical=True, feature_types=ftypes))
    tp = top1.spectrum_label.eq('TP'); fp = top1.spectrum_label.eq('FP')
    nmc = top1[tp].groupby('anno_name_lower')['wiki_id'].nunique(); dup = set(nmc[nmc >= 2].index) - {''}
    trust = (fp | (tp & top1.anno_name_lower.fillna('').isin(dup))).values
    iso = IsotonicRegression(out_of_bounds='clip'); iso.fit(raw_tr[trust], y[trust])

    def score(df, choa_col, simgap_col):
        d = df.copy(); d['compound_has_ok_adduct'] = df[choa_col]
        d['hit_isf_no_ok'] = ((d['hit_is_isf'] == 1) & (d[choa_col] == 0)).astype(int)
        d['sim_gap'] = df[simgap_col]
        X, _ = prep(d, cm)
        r = m.predict(xgb.DMatrix(X, enable_categorical=True, feature_types=ftypes))
        return r, iso.transform(r)

    def show(df, title):
        rb, cb = score(df, 'choa_within', 'sim_gap')          # baseline
        rf, cf = score(df, 'choa_fixed', 'simgap_fixed')      # both fixes
        out = df[['lib_name', 'adduct', 'hit_adduct_cat', 'entropy_similarity', 'sim_gap', 'rtc']].copy()
        out['conf_base'] = (cb*100).round(1); out['conf_fixed'] = (cf*100).round(1)
        print(f'\n=== {title} ===')
        print(out.sort_values('entropy_similarity', ascending=False).head(8).to_string(index=False))

    # Oliver bin 1: methylcytidine isomer tie
    b = ttof[ttof.wiki_id == 'aJWMOO8/7191']
    if len(b): show(b, "methylcytidine  aJWMOO8/7191  (sim_gap / isomer tie)")
    # Oliver bin 2: aminoadipic [M+H-H2O]+ ISF bins
    am = ttof[ttof['lib_name'].fillna('').str.lower().str.contains('aminoadipic') &
              ttof['adduct'].astype(str).str.replace(' ', '').str.contains('M\\+H-H2O', regex=True)]
    for w in am.wiki_id.unique()[:2]:
        show(ttof[ttof.wiki_id == w], f"aminoadipic ISF  {w}  (RT-confirmation)")

    print('\nNOTE: base model is Orbitrap-trained (cross-platform weak spot). '
          'rtc=1 means RT-confirmation fired WITHIN TTOF.')


if __name__ == '__main__':
    main()
