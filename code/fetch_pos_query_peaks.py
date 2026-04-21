"""
fetch_pos_query_peaks.py
-------------------------
Fetch query peaks for pos annotated spectra (3,738) from MassWiki.
Saves to data/pos_query_peaks_cache.json keyed by wiki_id.

Usage:
    python fetch_pos_query_peaks.py <MASSWIKI_TOKEN>
"""
import json, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = '/Users/ellayoung/Desktop/metabolo_confi_score'
POS_CSV   = f'{ROOT}/data/Orbitrap_HILIC_posESI_curated_041326.csv'
OUT_PEAKS = f'{ROOT}/data/pos_query_peaks_cache.json'

BASE_URL  = "https://masswiki.us-west-2.elasticbeanstalk.com"
ENDPOINT  = f"{BASE_URL}/analysis/get_data"
MAX_WORKERS = 8
RPS = 6.0
TIMEOUT = 30

TOKEN = sys.argv[1] if len(sys.argv) > 1 else os.getenv("MASSWIKI_TOKEN", "")
if not TOKEN:
    raise SystemExit("Need token")

# ── Load pos annotated wiki_ids (TP + FP only) ────────────────────────
print('Loading pos CSV...')
pos = pd.read_csv(POS_CSV, low_memory=False)
names = pos['name'].fillna('').astype(str)
yy = names.str.startswith('yy_')
zz = names.str.startswith('zz_')
annotated = pos[pos['is_manual_annotated'].astype(bool) &
                 ((~yy & ~zz & (names.str.strip() != '')) | yy)].copy()
wiki_ids = annotated['wiki_id'].dropna().astype(str).str.strip().unique().tolist()
print(f'Pos annotated (TP+FP): {len(wiki_ids):,}')

# ── Load existing peak cache ──────────────────────────────────────────
peaks_cache = {}
if os.path.exists(OUT_PEAKS):
    with open(OUT_PEAKS) as f:
        peaks_cache = json.load(f)
    print(f'Peak cache loaded: {len(peaks_cache):,}')

todo = [w for w in wiki_ids if w not in peaks_cache]
print(f'Need to fetch peaks: {len(todo):,}')
if not todo:
    print('All peaks cached.')
    sys.exit(0)

# ── Session ──────────────────────────────────────────────────────────
session = requests.Session()
retries = Retry(total=4, backoff_factor=0.4,
                status_forcelist=(429, 500, 502, 503, 504),
                allowed_methods=["GET"], respect_retry_after_header=True)
session.mount("https://", HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64))
session.headers.update({"Accept": "application/json", "Authorization": f"Bearer {TOKEN}"})

def extract_query_peaks(payload):
    """Match refetch_orbitrap_hits.py extract_query_peaks."""
    top_spec = payload.get("spectrum")
    if isinstance(top_spec, dict):
        for pk in ["peaks", "peaks_clean", "msms", "fragments"]:
            v = top_spec.get(pk)
            if v and isinstance(v, list):
                return v
    analysis = payload.get("analysis", {})
    for key in ["spectrum", "query_spectrum", "msms"]:
        spec = analysis.get(key)
        if isinstance(spec, dict):
            pk = spec.get("peaks") or spec.get("msms") or spec.get("fragments")
            if pk: return pk
        elif isinstance(spec, list) and spec:
            return spec
    for key in ["peaks", "msms"]:
        v = payload.get(key)
        if v and isinstance(v, list):
            return v
    return None

def fetch_one(wiki_id):
    for source, is_public in [("binbase", "false"), ("zyang2k", "true")]:
        try:
            r = session.get(ENDPOINT,
                            params={"wiki_id": wiki_id, "source": source, "isPublic": is_public},
                            timeout=TIMEOUT)
            if r.status_code == 200:
                peaks = extract_query_peaks(r.json())
                if peaks:
                    return wiki_id, peaks, None
                return wiki_id, None, f"{source}: no peaks in payload"
            elif r.status_code != 400:
                return wiki_id, None, f"{source}: HTTP {r.status_code}"
        except Exception as e:
            return wiki_id, None, f"{source}: {e}"
    return wiki_id, None, "both sources failed"

# ── Run in batches ───────────────────────────────────────────────────
BATCH = 500
errors = 0
t0 = time.time()
for bs in range(0, len(todo), BATCH):
    batch = todo[bs:bs+BATCH]
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(fetch_one, w): w for w in batch}
        for fut in as_completed(futs):
            wid, peaks, err = fut.result()
            if peaks:
                peaks_cache[wid] = peaks
            else:
                peaks_cache[wid] = None
                errors += 1
            time.sleep(1.0/RPS)
    with open(OUT_PEAKS, 'w') as f:
        json.dump(peaks_cache, f)
    done = min(bs+BATCH, len(todo))
    el = time.time()-t0
    rate = done/el if el > 0 else 0
    eta = (len(todo)-done)/rate if rate > 0 else 0
    print(f'  {done}/{len(todo)}  errors={errors}  {rate:.1f}req/s  ETA {eta/60:.0f}min')

n_with = sum(1 for v in peaks_cache.values() if v)
n_none = sum(1 for v in peaks_cache.values() if not v)
print(f'\nDone. Cache: {len(peaks_cache):,}. With peaks: {n_with:,}. None: {n_none:,}')
