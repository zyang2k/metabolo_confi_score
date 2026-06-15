"""check_rt_consistency.py — does the live rt_prediction endpoint agree with the
delta_predicted_rt MassWiki has ALREADY cached in its get_data analysis blob?

Within-MassWiki consistency only (ignores our local April data). For N spectra:
  1. GET /analysis/get_data  -> spectrum.rt + open_search candidates' (smiles, delta_predicted_rt)
     [open_search is used because identity_search currently returns empty]
  2. POST /reference_library/rt_prediction(smiles, method) -> fresh predicted RT
  3. cached_pred = spectrum.rt + cached_delta   (delta = predicted - observed, verified)
     compare fresh_pred vs cached_pred
Reports agreement distribution. If they disagree, the standalone rt_prediction service and
the model behind MassWiki's stored analysis are out of sync.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np, pandas as pd, requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = 'https://masswiki.us-west-2.elasticbeanstalk.com'
METHOD = {'neg': '5min_hilic_neg', 'pos': '5min_hilic_pos'}
N_SPECTRA = int(sys.argv[1]) if len(sys.argv) > 1 else 60


def main():
    tok = open('/tmp/.mw_token').read().strip()
    s = requests.Session()
    s.headers.update({'Authorization': f'Bearer {tok}', 'Accept': 'application/json',
                      'User-Agent': 'Mozilla/5.0', 'Content-Type': 'application/json'})
    # sample labeled spectra, mixed polarity
    ft = pd.read_csv(os.path.join(ROOT, 'data', 'feature_table_v2.csv'), low_memory=False,
                     usecols=['wiki_id', 'spectrum_label'])
    # polarity is NOT in the get_data JSON — map it from the hits file (neg/pos per wiki_id)
    oh = pd.read_csv(os.path.join(ROOT, 'data', 'orbitrap_hits_v2.csv.bak_apr23'),
                     low_memory=False, usecols=['wiki_id', 'polarity'])
    wid_pol = dict(oh.drop_duplicates('wiki_id').itertuples(index=False))
    wids = ft[ft['spectrum_label'].isin(['TP', 'FP'])]['wiki_id'].drop_duplicates()
    wids = [w for w in wids if wid_pol.get(w) in METHOD]
    wids = pd.Series(wids).sample(min(N_SPECTRA, len(wids)), random_state=0).tolist()

    rows = []          # (wiki_id, polarity, smiles, cached_delta, spec_rt)
    for i, wid in enumerate(wids):
        try:
            r = s.get(f'{BASE}/analysis/get_data',
                      params={'wiki_id': wid, 'source': 'binbase', 'isPublic': 'false'}, timeout=30)
            if r.status_code != 200:
                continue
            d = r.json()
        except Exception:
            continue
        spec_rt = (d.get('spectrum') or {}).get('rt')
        pol = wid_pol.get(wid)
        cands = (d.get('analysis', {}).get('reference_library', {}) or {}).get('open_search') or []
        if spec_rt is None or pol not in METHOD:
            continue
        for h in cands:
            sm = h.get('smiles'); dl = h.get('delta_predicted_rt')
            if sm and dl is not None:
                rows.append((wid, pol, sm, float(dl), float(spec_rt)))
        if i % 20 == 0:
            print(f'  get_data {i}/{len(wids)} ... collected {len(rows)} candidates', flush=True)

    df = pd.DataFrame(rows, columns=['wiki_id', 'polarity', 'smiles', 'cached_delta', 'spec_rt'])
    df = df.drop_duplicates(['smiles', 'polarity'])
    df.to_csv(os.path.join(ROOT, 'data', 'rt_consistency_candidates.csv'), index=False)
    print(f'\ncollected {len(df)} unique (smiles,polarity) with cached delta')

    # fresh predictions TODAY (separate recheck cache so we don't reuse stale values).
    cf = os.path.join(ROOT, 'data', 'rt_prediction_recheck_cache.json')
    pred_cache = json.load(open(cf)) if os.path.exists(cf) else {}
    for pol, g in df.groupby('polarity'):
        smis = [s for s in g['smiles'].tolist() if f'{METHOD[pol]}|{s}' not in pred_cache]
        for k in range(0, len(smis), 50):
            chunk = smis[k:k + 50]
            r = None
            for attempt in range(5):
                try:
                    r = s.post(f'{BASE}/reference_library/rt_prediction',
                               data=json.dumps({'smiles_list': chunk, 'method_id': METHOD[pol]}), timeout=180)
                except Exception:
                    time.sleep(3 * (attempt + 1)); continue
                if r.status_code == 401:
                    json.dump(pred_cache, open(cf, 'w'))
                    print('TOKEN EXPIRED mid-run — refresh and rerun (recheck cache saved)'); return
                if r.status_code in (502, 503, 504):
                    time.sleep(3 * (attempt + 1)); continue
                break
            if r is not None and r.status_code == 200:
                for sm, pr in zip(chunk, r.json().get('predictions', [])):
                    pred_cache[f'{METHOD[pol]}|{sm}'] = pr
            json.dump(pred_cache, open(cf, 'w'))
        print(f'  [{pol}] fresh predictions cached', flush=True)
    df['fresh_pred'] = [pred_cache.get(f'{METHOD[p]}|{s}') for p, s in zip(df['polarity'], df['smiles'])]
    print(f'fresh predictions obtained: {df["fresh_pred"].notna().sum()}/{len(df)}')
    df = df[df['fresh_pred'].notna()].copy()
    df['cached_pred'] = df['spec_rt'] + df['cached_delta']     # delta = pred - obs
    df['diff'] = df['fresh_pred'] - df['cached_pred']
    d = df['diff']
    print(f'\n=== fresh rt_prediction vs MassWiki cached prediction (n={len(df)}) ===')
    print(f'  median diff      : {d.median():.2f} s')
    print(f'  med |diff|       : {d.abs().median():.2f} s')
    print(f'  IQR              : [{d.quantile(.25):.1f}, {d.quantile(.75):.1f}] s')
    print(f'  within +/-2 s    : {100*(d.abs()<=2).mean():.1f}%')
    print(f'  within +/-5 s    : {100*(d.abs()<=5).mean():.1f}%')
    print(f'  within +/-10 s   : {100*(d.abs()<=10).mean():.1f}%')
    print(f'  max |diff|       : {d.abs().max():.1f} s')
    corr = np.corrcoef(df['fresh_pred'], df['cached_pred'])[0, 1]
    print(f'  correlation      : {corr:.4f}')
    df.to_csv(os.path.join(ROOT, 'data', 'rt_consistency_check.csv'), index=False)
    print('  wrote data/rt_consistency_check.csv')


if __name__ == '__main__':
    main()
