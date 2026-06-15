"""score_oliver_ttof.py — score Oliver's exact TTOF bin (2-aminoadipic acid [M+H-H2O]+)
through the updated pipeline, before vs after the RT-confirmation fix.

His example is TTOF HILIC pos (deployed MassWiki showed ~32%). Our pipeline is Orbitrap-native,
so we: build the TTOF-pos feature table (bench_ttof_golden machinery), locate the aminoadipic
[M+H-H2O]+ bin, check whether its clean [M+H] parent actually co-elutes in TTOF (the thing that
should cancel the ISF penalty), train the production GBM on Orbitrap, and score the bin with the
ISF penalty ON (orphan) vs OFF (confirmed). Reports raw + approx-calibrated.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.isotonic import IsotonicRegression
# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

import bench_ttof_golden as G
from bench_harness import BASE_NUMERIC, CATEGORICAL, XGB_PARAMS, NUM_BOOST_ROUND

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RT_WIN, ESIM_MIN = 10.0, 0.50


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
    lib_peaks = json.load(open(G.LIB_PEAKS))
    ttof = G.build_ttof_table(G.TTOF_POS_HITS, G.TTOF_POS_SPEC, G.TTOF_POS_PEAKS, lib_peaks, 1, tax)
    print(f'TTOF pos: {len(ttof):,} candidate rows, {ttof.wiki_id.nunique():,} bins')

    # locate aminoadipic [M+H-H2O]+ candidate(s)
    nm = ttof['lib_name'].fillna('').str.lower() if 'lib_name' in ttof else ttof['name'].fillna('').str.lower()
    am = ttof[nm.str.contains('aminoadipic')]
    isf = am[am['adduct'].astype(str).str.replace(' ', '').str.contains('M\\+H-H2O', case=False, regex=True)]
    print(f"\naminoadipic candidates: {len(am)};  [M+H-H2O]+ rows: {len(isf)}")
    if not len(isf):
        print('adducts seen for aminoadipic:', am['adduct'].value_counts().to_dict()); return
    # pick the highest-esim [M+H-H2O]+ (Oliver's hit, esim ~0.85)
    hit = isf.sort_values('entropy_similarity', ascending=False).iloc[0]
    ik = hit['hit_ik14']; wid = hit['wiki_id']; rt = hit['measured_rt']
    print(f"  bin {wid}  esim {hit['entropy_similarity']:.3f}  rt {rt:.1f}  adduct {hit['adduct']}  ik {ik}")

    # does a clean [M+H] of this compound co-elute in ANOTHER TTOF bin?  (the confirmation)
    same = ttof[(ttof.hit_ik14 == ik) & (ttof.hit_adduct_cat == 'ok') &
                (pd.to_numeric(ttof.entropy_similarity, errors='coerce') >= ESIM_MIN)]
    same = same[same.wiki_id != wid]
    same = same.assign(drt=(pd.to_numeric(same.measured_rt, errors='coerce') - rt).abs())
    coelut = same[same.drt <= RT_WIN]
    print(f"\nRT-confirmation check: clean-adduct same-compound bins co-eluting (±{RT_WIN}s): {len(coelut)}")
    if len(coelut):
        print(coelut[['wiki_id', 'adduct', 'entropy_similarity', 'measured_rt', 'drt']].head().to_string())
    confirmed = len(coelut) > 0

    # train production GBM on current Orbitrap table (RT-confirm-aware) + approx isotonic on trustworthy
    orb = pd.read_csv(G.FEATURE_TABLE, low_memory=False)
    lab = orb[orb.spectrum_label.isin(['TP', 'FP']) & orb.anno_ik14.fillna('').ne('')]
    top1 = orb.loc[lab.groupby('wiki_id').entropy_similarity.idxmax()].reset_index(drop=True)
    Xtr, cm = prep(top1); ytr = top1.hit_label.values
    ftypes = ['q']*len(BASE_NUMERIC)+['c']*len(CATEGORICAL)
    model = xgb.train(XGB_PARAMS, xgb.DMatrix(Xtr, label=ytr, enable_categorical=True, feature_types=ftypes), NUM_BOOST_ROUND)
    raw_tr = model.predict(xgb.DMatrix(Xtr, enable_categorical=True, feature_types=ftypes))
    # trustworthy slice for isotonic (yy_ FP + duplicate-name TP) — approx, in-sample
    tp = top1.spectrum_label.eq('TP'); fp = top1.spectrum_label.eq('FP')
    nmc = top1[tp].groupby('anno_name_lower')['wiki_id'].nunique(); dup = set(nmc[nmc >= 2].index) - {''}
    trust = (fp | (tp & top1.anno_name_lower.fillna('').isin(dup))).values
    iso = IsotonicRegression(out_of_bounds='clip'); iso.fit(raw_tr[trust], ytr[trust])

    # score Oliver's bin with ISF penalty ON (orphan) vs OFF (confirmed)
    base = pd.DataFrame([hit])
    def score_variant(has_ok, isf_no_ok):
        d = base.copy(); d['compound_has_ok_adduct'] = has_ok; d['hit_isf_no_ok'] = isf_no_ok
        X, _ = prep(d, cm)
        r = float(model.predict(xgb.DMatrix(X, enable_categorical=True, feature_types=ftypes))[0])
        return r, float(iso.transform([r])[0])
    r_on, c_on = score_variant(0, 1)     # orphan ISF (the old/penalized state)
    r_off, c_off = score_variant(1, 0)   # confirmed (penalty wiped)
    print('\n=== Oliver TTOF bin — 2-aminoadipic acid [M+H-H2O]+ ===')
    print(f'  ISF penalty ON  (orphan)    : raw {r_on*100:5.1f}%   calibrated {c_on*100:5.1f}%')
    print(f'  ISF penalty OFF (confirmed) : raw {r_off*100:5.1f}%   calibrated {c_off*100:5.1f}%')
    print(f'  RT-confirmation fires here? : {"YES -> confirmed score applies" if confirmed else "NO -> stays penalized"}')


if __name__ == '__main__':
    main()
