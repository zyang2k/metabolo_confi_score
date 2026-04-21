"""
fetch_ttof_hits_and_peaks.py
-----------------------------
Single-pass fetch: library hits + query peaks for TTOF HILIC annotated spectra.
Reads both pos and neg CSVs, fetches annotated (TP+FP) wiki_ids only.

Output:
  data/library_hits/ttof_hilic_neg_masswiki_hits.csv
  data/library_hits/ttof_hilic_pos_masswiki_hits.csv
  data/ttof_neg_query_peaks_cache.json
  data/ttof_pos_query_peaks_cache.json
  data/library_hits/ttof_fetch_cache.json  (resume support)

Usage:
    python fetch_ttof_hits_and_peaks.py <MASSWIKI_TOKEN>
"""
import json, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = '/Users/ellayoung/Desktop/metabolo_confi_score'
NEG_CSV = f'{ROOT}/data/TTOF_HILIC_negESI_uncurated_041326.csv'
POS_CSV = f'{ROOT}/data/TTOF_HILIC_posESI_uncurated_041326.csv'
CACHE_FILE = f'{ROOT}/data/library_hits/ttof_fetch_cache.json'

BASE_URL = "https://masswiki.us-west-2.elasticbeanstalk.com"
ENDPOINT = f"{BASE_URL}/analysis/get_data"
MAX_WORKERS = 8
RPS = 6.0
TIMEOUT = 30

TOKEN = sys.argv[1] if len(sys.argv) > 1 else os.getenv("MASSWIKI_TOKEN", "")
if not TOKEN:
    raise SystemExit("Need token: python fetch_ttof_hits_and_peaks.py <TOKEN>")


def extract_query_peaks(payload):
    """Extract raw query peaks from API response."""
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


def extract_hits(payload, wiki_id):
    """Extract ref + anno hits from API response."""
    rows = []
    analysis = payload.get("analysis", {})
    for hit_source in ("reference_library", "annotation_library"):
        lib = analysis.get(hit_source) or {}
        hits = lib.get("identity_search") or []
        if not isinstance(hits, list):
            continue
        source_label = "reference" if hit_source == "reference_library" else "annotation"
        for i, h in enumerate(hits, 1):
            rows.append({
                "wiki_id": wiki_id,
                "hit_source": source_label,
                "db": h.get("db") or h.get("source"),
                "id": h.get("id") or h.get("identifier"),
                "lib_name": h.get("name"),
                "adduct": h.get("adduct"),
                "lib_precursor_mz": h.get("precursor_mz") or h.get("precursor"),
                "entropy_similarity": h.get("entropy_similarity") or h.get("score"),
                "library_type": h.get("library_type"),
                "lib_rt": h.get("rt") or h.get("retention_time"),
                "ri": h.get("ri"),
                "rank": h.get("rank") or i,
                "smiles": h.get("smiles"),
                "predicted_rt_hilic": h.get("predicted_rt_hilic"),
                "predicted_rt_rp": h.get("predicted_rt_rp"),
                "anno_delta_rt": h.get("anno_delta_rt"),
                "delta_predicted_rt": h.get("delta_predicted_rt"),
                "library_wiki_id": h.get("library_wiki_id") or h.get("wiki_id"),
            })
    return rows


# ── 1. Collect annotated wiki_ids from both files ─────────────────────
print('Loading TTOF CSVs...')
wid_to_mode = {}  # wiki_id -> 'neg' or 'pos'

for csv_path, mode in [(NEG_CSV, 'neg'), (POS_CSV, 'pos')]:
    df = pd.read_csv(csv_path, low_memory=False)
    names = df['name'].fillna('').astype(str)
    yy = names.str.startswith('yy_')
    zz = names.str.startswith('zz_')
    annotated = df[df['is_manual_annotated'].astype(bool) &
                    ((~yy & ~zz & (names.str.strip() != '')) | yy)].copy()
    for wid in annotated['wiki_id'].dropna().astype(str).str.strip().unique():
        wid_to_mode[wid] = mode
    print(f'  {mode}: {len(annotated)} annotated spectra')

all_wids = list(wid_to_mode.keys())
print(f'Total annotated wiki_ids: {len(all_wids):,}')

# ── 2. Load cache ────────────────────────────────────────────────────
cache = {}
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE) as f:
        cache = json.load(f)
    print(f'Cache loaded: {len(cache):,}')

todo = [w for w in all_wids if w not in cache]
print(f'Need to fetch: {len(todo):,}')

if todo:
    # ── 3. Session ───────────────────────────────────────────────────
    session = requests.Session()
    retries = Retry(total=4, backoff_factor=0.4,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=["GET"], respect_retry_after_header=True)
    session.mount("https://", HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64))
    session.headers.update({"Accept": "application/json", "Authorization": f"Bearer {TOKEN}"})

    def fetch_one(wiki_id):
        for source, is_public in [("binbase", "false"), ("zyang2k", "true")]:
            try:
                r = session.get(ENDPOINT,
                                params={"wiki_id": wiki_id, "source": source, "isPublic": is_public},
                                timeout=TIMEOUT)
                if r.status_code == 200:
                    payload = r.json()
                    hits = extract_hits(payload, wiki_id)
                    peaks = extract_query_peaks(payload)
                    return wiki_id, {"hits": hits, "peaks": peaks}, None
                elif r.status_code != 400:
                    return wiki_id, None, f"{source}: HTTP {r.status_code}"
            except Exception as e:
                return wiki_id, None, f"{source}: {e}"
        return wiki_id, None, "both sources failed"

    # ── 4. Batch fetch ───────────────────────────────────────────────
    BATCH = 500
    errors = 0
    t0 = time.time()
    for bs in range(0, len(todo), BATCH):
        batch = todo[bs:bs+BATCH]
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(fetch_one, w): w for w in batch}
            for fut in as_completed(futs):
                wid, result, err = fut.result()
                if result is not None:
                    cache[wid] = result
                else:
                    cache[wid] = None
                    errors += 1
                time.sleep(1.0/RPS)
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache, f)
        done = min(bs+BATCH, len(todo))
        el = time.time()-t0
        rate = done/el if el > 0 else 0
        eta = (len(todo)-done)/rate if rate > 0 else 0
        print(f'  {done}/{len(todo)}  errors={errors}  {rate:.1f}req/s  ETA {eta/60:.0f}min')

    print(f'\nFetch complete. Errors: {errors}')

# ── 5. Flatten to per-mode CSVs + peak caches ────────────────────────
print('Flattening...')
for mode in ['neg', 'pos']:
    mode_wids = [w for w, m in wid_to_mode.items() if m == mode]
    hit_rows = []
    peaks_dict = {}
    n_with, n_without = 0, 0

    for wid in mode_wids:
        entry = cache.get(wid)
        if not entry:
            n_without += 1
            continue
        hits = entry.get('hits') or []
        peaks = entry.get('peaks')
        if hits:
            hit_rows.extend(hits)
            n_with += 1
        else:
            n_without += 1
        if peaks:
            peaks_dict[wid] = peaks

    hits_df = pd.DataFrame(hit_rows)
    hits_path = f'{ROOT}/data/library_hits/ttof_hilic_{mode}_masswiki_hits.csv'
    hits_df.to_csv(hits_path, index=False)

    peaks_path = f'{ROOT}/data/ttof_{mode}_query_peaks_cache.json'
    with open(peaks_path, 'w') as f:
        json.dump(peaks_dict, f)

    print(f'  {mode}: {len(hits_df):,} hits ({n_with} spectra), '
          f'{len(peaks_dict):,} with peaks, {n_without} no hits → {hits_path}')

print('\nDone.')
