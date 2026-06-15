"""refresh_rt_predictions.py — Path B RT refresh via MassWiki rt_prediction endpoint.

Background (memory: project_masswiki_identity_search_broken_20260611). MassWiki's
get_data returns identity_search empty (RT-prediction backend was down; Fanzhou's service
restarted on a new IP). Rather than reanalyze (a shared-prod write), we refresh RT
predictions directly: POST /reference_library/rt_prediction {smiles_list, method_id} for
every candidate SMILES, keeping the intact Apr-23 candidates + labels.

Validated 2026-06-12:
  * sign/units: signed_delta_rt = predicted_rt - observed_rt, both in SECONDS
    (confirmed vs the API's own delta_predicted_rt: pred 19.89 - obs 99.78 = -79.89 ≈ -77.8).
  * method ids: 5min_hilic_neg / 5min_hilic_pos (GET /method/get_all_list).
  * CAVEAT: the restarted rt_prediction service may serve a DIFFERENT model version than
    the one behind the Apr-23 stored delta_predicted_rt (disagreements up to ~14 s) — a
    full refresh is internally consistent but shifts the RT distribution; confirm with
    Fanzhou which model is live.

Resume-safe: predictions cached in data/rt_prediction_cache.json keyed by
'{method_id}|{smiles}'. Re-run after a token refresh to continue. On 401 it stops cleanly
and reports remaining count. This script ONLY builds the prediction cache; applying it to
the feature table is a separate step (apply_rt_refresh, below / --apply).
"""
from __future__ import annotations
import os, sys, json, time, argparse
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HITS_BAK = os.path.join(ROOT, 'data', 'orbitrap_hits_v2.csv.bak_apr23')
FEATURE_TABLE = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
CACHE = os.path.join(ROOT, 'data', 'rt_prediction_cache.json')
TOKEN_FILE = '/tmp/.mw_token'
BASE = 'https://masswiki.us-west-2.elasticbeanstalk.com'
EP = f'{BASE}/reference_library/rt_prediction'
METHOD = {'neg': '5min_hilic_neg', 'pos': '5min_hilic_pos'}
BATCH = 50
SAVE_EVERY = 5  # batches
RETRIES = 4      # per-batch retries on 502/504 (service flaky since restart)
BACKOFF = 3.0    # seconds, exponential


def load_cache() -> dict:
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            return json.load(f)
    return {}


def save_cache(c: dict):
    tmp = CACHE + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(c, f)
    os.replace(tmp, CACHE)


def candidate_smiles() -> pd.DataFrame:
    """Unique (smiles, polarity) pairs that appear as candidates in the feature table."""
    oh = pd.read_csv(HITS_BAK, low_memory=False,
                     usecols=['wiki_id', 'library_wiki_id', 'smiles', 'polarity'])
    oh = oh.dropna(subset=['smiles']).drop_duplicates(['wiki_id', 'library_wiki_id'])
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False, usecols=['wiki_id', 'library_wiki_id'])
    cand = ft.merge(oh, on=['wiki_id', 'library_wiki_id'], how='inner')
    return cand[['smiles', 'polarity']].drop_duplicates().reset_index(drop=True)


def predict_all():
    pairs = candidate_smiles()
    cache = load_cache()
    sess = requests.Session()
    tok = open(TOKEN_FILE).read().strip()
    sess.headers.update({'Authorization': f'Bearer {tok}', 'Accept': 'application/json',
                         'User-Agent': 'Mozilla/5.0', 'Content-Type': 'application/json'})

    todo = [(s, p) for s, p in pairs.itertuples(index=False)
            if f'{METHOD[p]}|{s}' not in cache]
    print(f'{len(pairs):,} unique (smiles,polarity); {len(cache):,} cached; {len(todo):,} to predict')

    # group remaining by polarity so each POST uses one method_id
    by_pol = {'neg': [s for s, p in todo if p == 'neg'],
              'pos': [s for s, p in todo if p == 'pos']}
    n_done = 0
    batch_i = 0
    for pol, smis in by_pol.items():
        mid = METHOD[pol]
        for k in range(0, len(smis), BATCH):
            chunk = smis[k:k + BATCH]
            r = None
            for attempt in range(RETRIES):
                try:
                    r = sess.post(EP, data=json.dumps({'smiles_list': chunk, 'method_id': mid}),
                                  timeout=180)
                except Exception as e:
                    print(f'  [{pol}] batch {batch_i} req-exc {e} (try {attempt})')
                    time.sleep(BACKOFF * (attempt + 1))
                    continue
                if r.status_code == 401:
                    save_cache(cache)
                    print(f'\nTOKEN EXPIRED (401). Saved {len(cache):,} cached. '
                          f'{len(todo)-n_done:,} remaining — refresh /tmp/.mw_token and re-run.')
                    return False
                if r.status_code in (502, 503, 504):
                    time.sleep(BACKOFF * (attempt + 1))   # flaky backend; back off + retry
                    continue
                break
            if r is None or r.status_code != 200:
                code = r.status_code if r is not None else 'exc'
                print(f'  [{pol}] batch {batch_i} HTTP {code} after {RETRIES} tries — leaving uncached')
                save_cache(cache)
                batch_i += 1
                continue
            preds = r.json().get('predictions', [])
            if len(preds) != len(chunk):
                print(f'  [{pol}] length mismatch {len(preds)}!={len(chunk)} — skipping batch')
                batch_i += 1
                continue
            for s, pr in zip(chunk, preds):
                cache[f'{mid}|{s}'] = pr
            n_done += len(chunk)
            batch_i += 1
            if batch_i % SAVE_EVERY == 0:
                save_cache(cache)
                print(f'  [{pol}] {n_done:,}/{len(todo):,} predicted', flush=True)
    save_cache(cache)
    print(f'\nDONE. cache now {len(cache):,} entries. {len(todo)-n_done:,} still missing.')
    return True


def apply_refresh():
    """Recompute signed_delta_rt = predicted - observed for every candidate row, write
    a refreshed feature table (backs up the current one first)."""
    cache = load_cache()
    oh = pd.read_csv(HITS_BAK, low_memory=False,
                     usecols=['wiki_id', 'library_wiki_id', 'smiles', 'polarity'])
    oh = oh.dropna(subset=['smiles']).drop_duplicates(['wiki_id', 'library_wiki_id'])
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    m = ft.merge(oh[['wiki_id', 'library_wiki_id', 'smiles', 'polarity']],
                 on=['wiki_id', 'library_wiki_id'], how='left')

    def pred_for(row):
        s, pol = row['smiles'], row['polarity_y'] if 'polarity_y' in row else row.get('polarity')
        if not isinstance(s, str) or pol not in METHOD:
            return None
        return cache.get(f'{METHOD[pol]}|{s}')

    pol_col = 'polarity_y' if 'polarity_y' in m.columns else 'polarity'
    preds = [cache.get(f'{METHOD[p]}|{s}') if (isinstance(s, str) and p in METHOD) else None
             for s, p in zip(m['smiles'], m[pol_col])]
    m['new_pred_rt'] = preds
    new_sdrt = m['new_pred_rt'] - m['rt_obs']
    cov = m['new_pred_rt'].notna().mean() * 100
    old = pd.to_numeric(ft['signed_delta_rt'], errors='coerce')
    ft_out = ft.copy()
    ft_out['signed_delta_rt'] = new_sdrt.where(m['new_pred_rt'].notna(), other=pd.NA).values
    newcov = pd.to_numeric(ft_out['signed_delta_rt'], errors='coerce').notna().mean() * 100
    bak = FEATURE_TABLE + '.bak_pre_rt_refresh'
    if not os.path.exists(bak):
        pd.read_csv(FEATURE_TABLE, low_memory=False).to_csv(bak, index=False)
    ft_out.to_csv(FEATURE_TABLE, index=False)
    print(f'Applied RT refresh: new_pred coverage {cov:.1f}% of candidate rows; '
          f'signed_delta_rt non-null {newcov:.1f}% (was {old.notna().mean()*100:.1f}%). '
          f'Backup → {bak}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='apply cached predictions to feature table')
    args = ap.parse_args()
    if args.apply:
        apply_refresh()
    else:
        predict_all()
