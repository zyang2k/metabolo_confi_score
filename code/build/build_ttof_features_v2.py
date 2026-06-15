"""build_ttof_features_v2.py — run build_features_v2's REAL per-candidate logic on TTOF.

The shortcut builder (bench_ttof_golden) used a per-BIN sim_gap, broken adduct classification,
and dropped no-predicted-RT candidates. This replicates the Orbitrap builder faithfully on TTOF
inputs (reusing build_features_v2 helpers), so Oliver's TTOF bins can be scored the way Orbitrap
bins are. MS2 cosines are skipped (left NaN — minor, correlated with esim); everything else matches.

Scores with the Orbitrap-trained production model (cross-platform — absolute numbers carry that caveat).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, xgboost as xgb
from sklearn.isotonic import IsotonicRegression
# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

import build_features_v2 as B
from bench_harness import BASE_NUMERIC, CATEGORICAL, XGB_PARAMS, NUM_BOOST_ROUND

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RT_WIN, ESIM_MIN = 10.0, 0.50
HITS = {'pos': ROOT+'/data/library_hits/ttof_hilic_pos_masswiki_hits.csv',
        'neg': ROOT+'/data/library_hits/ttof_hilic_neg_masswiki_hits.csv'}
SPEC = {'pos': ROOT+'/data/TTOF_HILIC_posESI_uncurated_041326.csv',
        'neg': ROOT+'/data/TTOF_HILIC_negESI_uncurated_041326.csv'}


def compute_sim_gaps(group):
    comp = group.loc[group['hit_ik14'] != '', 'entropy_similarity'].values
    sims = group['entropy_similarity'].values
    group = group.copy()
    if len(comp) <= 1:
        group['sim_gap'] = sims if len(comp) == 1 else 0.0  # single → distance from 0
        if len(comp) == 1:
            group['sim_gap'] = [s if (ik != '' and s == comp[0]) else 0.0
                                for s, ik in zip(sims, group['hit_ik14'])]
        return group
    sc = np.sort(comp)[::-1]; top, run = sc[0], sc[1]
    group['sim_gap'] = [(s - run) if (ik != '' and s == top) else (s - top)
                        for s, ik in zip(sims, group['hit_ik14'])]
    return group


def build_ttof(pol):
    hits = pd.read_csv(HITS[pol], low_memory=False)
    spec = pd.read_csv(SPEC[pol], low_memory=False)
    sm = spec[['wiki_id', 'precursor_mz', 'rt', 'entropy']].drop_duplicates('wiki_id').rename(
        columns={'precursor_mz': 'obs_mz', 'rt': 'rt_obs', 'entropy': 'spectral_entropy'})
    df = hits.merge(sm, on='wiki_id', how='left')
    df['polarity'] = 0 if pol == 'neg' else 1
    # IK14 + neutral mass via unique-SMILES cache (rdkit per-row is too slow over ~300k)
    uniq = pd.Series(df['smiles'].dropna().unique())
    ik = {s: B.get_ik14(s) for s in uniq}; nm = {s: B.get_neutral_mass(s) for s in uniq}
    df['hit_ik14'] = df['smiles'].map(ik).fillna('')
    df = df[df['hit_ik14'] != ''].copy()
    lookup = B._load_adduct_taxonomy()
    df['hit_adduct_cat'] = df['adduct'].apply(lambda a: B.classify_adduct(a, lookup))
    df['hit_is_isf'] = (df['hit_adduct_cat'] == 'isf').astype(int)
    df['hit_is_dubious'] = (df['hit_adduct_cat'] == 'dubious').astype(int)
    df['hit_neutral_mass'] = df['smiles'].map(nm)
    df['hit_theoretical_mz'] = [B.compute_theoretical_mz(m, a) if pd.notna(m) else np.nan
                                for m, a in zip(df['hit_neutral_mass'], df['adduct'])]
    df['delta_mda'] = (pd.to_numeric(df['obs_mz'], errors='coerce') - df['hit_theoretical_mz']).abs() * 1000
    df['signed_delta_rt'] = pd.to_numeric(df['delta_predicted_rt'], errors='coerce')
    df['_is_ok'] = (df['hit_adduct_cat'] == 'ok').astype(int)
    # compound_has_ok_adduct (within bin) + n_candidate_adducts — BEFORE dedup (over all adduct rows)
    choa = df.groupby(['wiki_id', 'hit_ik14'])['_is_ok'].max().rename('compound_has_ok_adduct').reset_index()
    nadd = df.groupby(['wiki_id', 'hit_ik14'])['adduct'].nunique().rename('n_candidate_adducts').reset_index()
    df = df.sort_values('entropy_similarity', ascending=False).drop_duplicates(['wiki_id', 'hit_ik14'], keep='first')
    df = df.merge(choa, on=['wiki_id', 'hit_ik14'], how='left').merge(nadd, on=['wiki_id', 'hit_ik14'], how='left')
    # RT-confirmation within TTOF
    conf = df[(df['_is_ok'] == 1) & (pd.to_numeric(df['entropy_similarity'], errors='coerce') >= ESIM_MIN)
              & df['rt_obs'].notna()]
    by_ik = {k: g[['wiki_id', 'rt_obs']].values for k, g in conf.groupby('hit_ik14')}
    def rtc(w, ik, rt):
        arr = by_ik.get(ik)
        return int(arr is not None and np.isfinite(rt) and any(w2 != w and abs(rt - r2) <= RT_WIN for w2, r2 in arr))
    df['compound_ok_rt_confirmed'] = [rtc(w, ik, rt) for w, ik, rt in zip(df.wiki_id, df.hit_ik14, df.rt_obs)]
    df['compound_has_ok_adduct'] = ((df['compound_has_ok_adduct'] == 1) | (df['compound_ok_rt_confirmed'] == 1)).astype(int)
    df['hit_isf_no_ok'] = ((df['hit_is_isf'] == 1) & (df['compound_has_ok_adduct'] == 0)).astype(int)
    df = df.groupby('wiki_id', group_keys=False).apply(compute_sim_gaps)
    df['n_candidates'] = df.groupby('wiki_id')['hit_ik14'].transform('count')
    for c in ['forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int']:
        df[c] = np.nan   # MS2 cosines skipped
    df['db'] = df.get('db', 'ttof'); df['name'] = df.get('lib_name', '')
    return df


def main():
    print('Building faithful TTOF pos features (real per-candidate logic)...')
    t = build_ttof('pos')
    print(f'  {len(t):,} candidate rows, {t.wiki_id.nunique():,} bins   '
          f'adduct cats: {t.hit_adduct_cat.value_counts().to_dict()}')

    def prep(df, cm=None):
        X = pd.DataFrame(index=df.index); fit = cm is None; cm = cm or {}
        for c in BASE_NUMERIC: X[c] = pd.to_numeric(df[c], errors='coerce') if c in df else np.nan
        for c in CATEGORICAL:
            s = df[c].astype(str).fillna('missing') if c in df else pd.Series('missing', index=df.index)
            if fit: cm[c] = {v: i for i, v in enumerate(sorted(s.unique()))}
            mx = max(cm[c].values())+1 if cm[c] else 0
            X[c] = s.map(lambda v: cm[c].get(v, mx)).astype('int32')
        return X[BASE_NUMERIC+CATEGORICAL], cm
    orb = pd.read_csv(ROOT+'/data/feature_table_v2.csv', low_memory=False)
    lab = orb[orb.spectrum_label.isin(['TP', 'FP']) & orb.anno_ik14.fillna('').ne('')]
    top1 = orb.loc[lab.groupby('wiki_id').entropy_similarity.idxmax()].reset_index(drop=True)
    Xtr, cm = prep(top1); y = top1.hit_label.values
    fty = ['q']*len(BASE_NUMERIC)+['c']*len(CATEGORICAL)
    m = xgb.train(XGB_PARAMS, xgb.DMatrix(Xtr, label=y, enable_categorical=True, feature_types=fty), NUM_BOOST_ROUND)
    rtr = m.predict(xgb.DMatrix(Xtr, enable_categorical=True, feature_types=fty))
    tp = top1.spectrum_label.eq('TP'); fp = top1.spectrum_label.eq('FP')
    nmc = top1[tp].groupby('anno_name_lower')['wiki_id'].nunique(); dup = set(nmc[nmc >= 2].index) - {''}
    trust = (fp | (tp & top1.anno_name_lower.fillna('').isin(dup))).values
    iso = IsotonicRegression(out_of_bounds='clip'); iso.fit(rtr[trust], y[trust])
    X, _ = prep(t, cm)
    t['conf'] = (iso.transform(m.predict(xgb.DMatrix(X, enable_categorical=True, feature_types=fty)))*100).round(1)

    # Oliver's methylcytidine bin
    b = t[t.wiki_id == 'aJWMOO8/7191'].sort_values('entropy_similarity', ascending=False)
    print('\n=== Oliver methylcytidine bin aJWMOO8/7191 (FAITHFUL per-candidate) ===')
    print(b[['name', 'adduct', 'hit_adduct_cat', 'entropy_similarity', 'sim_gap', 'signed_delta_rt',
             'compound_ok_rt_confirmed', 'conf']].head(8).to_string(index=False))


if __name__ == '__main__':
    main()
