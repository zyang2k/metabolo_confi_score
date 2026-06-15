"""fetch_blank_query_peaks.py — Fetch raw query peaks for blank (unannotated) spectra.

Blanks = spectra in the curated CSVs with identity_score < 0.7 and no curator name.
We have library hits for these (data/{neg,pos}_blanks_hits.csv) but no raw peaks,
which blocks computing Oliver's noise features (n_ions, normalized_entropy, etc.)
on the unannotated pool.

This script reuses the same MassWiki API endpoint and binbase→zyang2k fallback
pattern as code/fetch_pos_query_peaks.py.

Token is read from MASSWIKI_TOKEN env var (preferred) or argv[2]. Never written
to disk.

Usage:
    MASSWIKI_TOKEN=<token> python code/fetch_blank_query_peaks.py --polarity pos
    MASSWIKI_TOKEN=<token> python code/fetch_blank_query_peaks.py --polarity neg
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CUR_TEMPLATE = os.path.join(ROOT, 'data', 'Orbitrap_HILIC_{POL}ESI_curated_042126.csv')
OUT_TEMPLATE = os.path.join(ROOT, 'data', 'blanks_query_peaks_cache_{pol}.json')
ERR_TEMPLATE = os.path.join(ROOT, 'data', 'blanks_query_peaks_errors_{pol}.csv')

BASE_URL = 'https://masswiki.us-west-2.elasticbeanstalk.com'
ENDPOINT = f'{BASE_URL}/analysis/get_data'
MAX_WORKERS = 8
RPS = 6.0
TIMEOUT = 30
BATCH = 500


def extract_query_peaks(payload):
    """Match the extractor pattern used in fetch_pos_query_peaks.py."""
    top_spec = payload.get('spectrum')
    if isinstance(top_spec, dict):
        for pk in ['peaks', 'peaks_clean', 'msms', 'fragments']:
            v = top_spec.get(pk)
            if v and isinstance(v, list):
                return v
    analysis = payload.get('analysis', {})
    for key in ['spectrum', 'query_spectrum', 'msms']:
        spec = analysis.get(key)
        if isinstance(spec, dict):
            pk = spec.get('peaks') or spec.get('msms') or spec.get('fragments')
            if pk:
                return pk
        elif isinstance(spec, list) and spec:
            return spec
    for key in ['peaks', 'msms']:
        v = payload.get(key)
        if v and isinstance(v, list):
            return v
    return None


def make_fetch_one(session: requests.Session):
    def fetch_one(wiki_id: str):
        for source, is_public in [('binbase', 'false'), ('zyang2k', 'true')]:
            try:
                r = session.get(
                    ENDPOINT,
                    params={'wiki_id': wiki_id, 'source': source, 'isPublic': is_public},
                    timeout=TIMEOUT,
                )
                if r.status_code == 200:
                    peaks = extract_query_peaks(r.json())
                    if peaks:
                        return wiki_id, peaks, None
                    # Fallthrough to next source on empty payload
                    return wiki_id, None, f'{source}: no peaks in payload'
                elif r.status_code != 400:
                    return wiki_id, None, f'{source}: HTTP {r.status_code}'
            except Exception as e:
                return wiki_id, None, f'{source}: {e}'
        return wiki_id, None, 'both sources failed'
    return fetch_one


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--polarity', required=True, choices=['pos', 'neg'])
    p.add_argument('--token-positional', nargs='?', default=None,
                   help='Optional CLI fallback for token (prefer MASSWIKI_TOKEN env)')
    args = p.parse_args()

    token = os.environ.get('MASSWIKI_TOKEN', '') or (args.token_positional or '')
    if not token:
        raise SystemExit('Set MASSWIKI_TOKEN env var (or pass token positionally).')

    cur_path = CUR_TEMPLATE.format(POL=args.polarity)
    out_path = OUT_TEMPLATE.format(pol=args.polarity)
    err_path = ERR_TEMPLATE.format(pol=args.polarity)

    print(f'Loading curated {args.polarity} CSV: {cur_path}')
    df = pd.read_csv(cur_path, low_memory=False)
    blanks = df[
        (df['identity_score'] < 0.7) &
        (df['name'].isna() | (df['name'].astype(str).str.strip() == ''))
    ]
    wids = (blanks['wiki_id'].dropna().astype(str).str.strip()
            .unique().tolist())
    wids = [w for w in wids if w]
    print(f'  blanks to fetch peaks for: {len(wids):,}')

    # Resume from existing cache
    cache: dict = {}
    if os.path.exists(out_path):
        with open(out_path) as f:
            cache = json.load(f)
        print(f'  existing cache: {len(cache):,}')
    todo = [w for w in wids if w not in cache]
    print(f'  remaining to fetch: {len(todo):,}')
    if not todo:
        print('All peaks cached. Nothing to do.')
        return

    # Session with retry + auth header (token only used in this process)
    session = requests.Session()
    retries = Retry(total=4, backoff_factor=0.4,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=['GET'], respect_retry_after_header=True)
    session.mount('https://', HTTPAdapter(max_retries=retries,
                                          pool_connections=64, pool_maxsize=64))
    session.headers.update({
        'Accept': 'application/json',
        'Authorization': f'Bearer {token}',
    })

    fetch_one = make_fetch_one(session)
    err_rows = []
    errors = 0
    t0 = time.time()

    for bs in range(0, len(todo), BATCH):
        batch = todo[bs:bs + BATCH]
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(fetch_one, w): w for w in batch}
            for fut in as_completed(futs):
                wid, peaks, err = fut.result()
                if peaks:
                    cache[wid] = peaks
                else:
                    cache[wid] = None
                    errors += 1
                    err_rows.append({'wiki_id': wid, 'error': err})
                time.sleep(1.0 / RPS)
        # Periodic flush so we don't lose progress
        with open(out_path, 'w') as f:
            json.dump(cache, f)
        done = min(bs + BATCH, len(todo))
        el = time.time() - t0
        rate = done / el if el > 0 else 0
        eta = (len(todo) - done) / rate if rate > 0 else 0
        print(f'  [{args.polarity}] {done}/{len(todo)}  errors={errors}  '
              f'{rate:.1f}req/s  ETA {eta/60:.1f}min', flush=True)

    if err_rows:
        pd.DataFrame(err_rows).to_csv(err_path, index=False)
        print(f'  wrote errors to {err_path}')

    n_with = sum(1 for v in cache.values() if v)
    n_none = sum(1 for v in cache.values() if not v)
    print(f'\nDone. cache size {len(cache):,}  with peaks {n_with:,}  empty {n_none:,}')


if __name__ == '__main__':
    main()
