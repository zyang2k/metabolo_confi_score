"""
fetch_orbitrap_hits.py
---------------------
Fetch ALL reference + annotation library hits from MassWiki for the
Orbitrap HILIC negESI dataset (5,971 spectra).

No filtering by annotation status — we want hits for every spectrum.

Output: data/library_hits/orbitrap_hilic_neg_masswiki_hits.csv
Errors: data/library_hits/orbitrap_hilic_neg_masswiki_errors.csv

Uses ThreadPoolExecutor (8 workers), ~6 RPS throttle, auto-retry on
429/500-504, fallback from binbase→zyang2k on 400.
"""
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = '/Users/ellayoung/Desktop/metabolo_confi_score'
SPREADSHEET = f'{ROOT}/data/masswiki_Orbitrap HILIC negESI_2026-03-19.xlsx'
OUT_HITS = f'{ROOT}/data/library_hits/orbitrap_hilic_neg_masswiki_hits.csv'
OUT_ERRORS = f'{ROOT}/data/library_hits/orbitrap_hilic_neg_masswiki_errors.csv'
CACHE_FILE = f'{ROOT}/data/library_hits/orbitrap_fetch_cache.json'

BASE_URL = "https://masswiki.us-west-2.elasticbeanstalk.com"
ENDPOINT = f"{BASE_URL}/analysis/get_data"
MAX_WORKERS = 8
RPS = 6.0
TIMEOUT = 20

TOKEN = sys.argv[1] if len(sys.argv) > 1 else os.getenv("MASSWIKI_TOKEN", "")

# ── 1. Load wiki_ids ─────────────────────────────────────────────────────
print('Loading Orbitrap spreadsheet...')
df = pd.read_excel(SPREADSHEET, sheet_name='masswiki_result_2026-03-19', header=4)
all_wids = df['wiki_id'].dropna().astype(str).str.strip().unique().tolist()
print(f'Total unique wiki_ids: {len(all_wids)}')

# ── 2. Load cache (resume support) ──────────────────────────────────────
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE) as f:
        cache = json.load(f)
    print(f'Cache loaded: {len(cache)} entries')
else:
    cache = {}

todo = [w for w in all_wids if w not in cache]
print(f'Remaining to fetch: {len(todo)}')

if not todo:
    print('Nothing to fetch — all cached.')
else:
    # ── 3. Build session ─────────────────────────────────────────────────
    session = requests.Session()
    retries = Retry(total=4, backoff_factor=0.4,
                    status_forcelist=(429, 500, 502, 503, 504),
                    allowed_methods=["GET"], respect_retry_after_header=True)
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    session.mount("https://", adapter)
    session.headers.update({
        "Accept": "application/json",
        "Authorization": f"Bearer {TOKEN}",
    })

    # ── 4. Fetch function ────────────────────────────────────────────────
    def fetch_one(wiki_id):
        for source, is_public in [("binbase", "false"), ("zyang2k", "true")]:
            try:
                r = session.get(ENDPOINT,
                                params={"wiki_id": wiki_id, "source": source, "isPublic": is_public},
                                timeout=TIMEOUT)
                if r.status_code == 200:
                    payload = r.json()
                    analysis = payload.get("analysis", {}) if isinstance(payload, dict) else {}
                    ref = (analysis.get("reference_library") or {}).get("identity_search")
                    anno = (analysis.get("annotation_library") or {}).get("identity_search")
                    if isinstance(ref, list) or isinstance(anno, list):
                        return wiki_id, {"ref": ref or [], "anno": anno or []}, None
                elif r.status_code != 400:
                    return wiki_id, None, f"{source}/{is_public}: HTTP {r.status_code}"
                # 400 = try fallback
            except Exception as e:
                return wiki_id, None, f"{source}/{is_public}: {e}"
        return wiki_id, None, "both sources failed"

    # ── 5. Run in batches (save cache every 500) ─────────────────────────
    BATCH = 500
    errors = []
    t0 = time.time()

    for batch_start in range(0, len(todo), BATCH):
        batch = todo[batch_start:batch_start + BATCH]

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(fetch_one, w): w for w in batch}
            for fut in as_completed(futs):
                wid, result, err = fut.result()
                if result is not None:
                    cache[wid] = result
                else:
                    cache[wid] = None
                    errors.append({"wiki_id": wid, "error": err})
                time.sleep(1.0 / RPS)

        # Save cache
        with open(CACHE_FILE, 'w') as f:
            json.dump(cache, f)

        done = min(batch_start + BATCH, len(todo))
        elapsed = time.time() - t0
        rate = done / elapsed if elapsed > 0 else 0
        eta = (len(todo) - done) / rate if rate > 0 else 0
        print(f'  {done}/{len(todo)} done  ({len(errors)} errors)  '
              f'{rate:.1f} req/s  ETA {eta/60:.0f}min')

    print(f'\nFetch complete. Errors: {len(errors)}')

# ── 6. Flatten to CSV ────────────────────────────────────────────────────
print('Flattening hits to CSV...')
rows = []
for wid, hit_bundle in cache.items():
    if not hit_bundle:
        continue
    for source_key, api_key in [("reference", "ref"), ("annotation", "anno")]:
        hits = hit_bundle.get(api_key) or []
        for i, h in enumerate(hits, 1):
            rows.append({
                "wiki_id": wid,
                "hit_source": source_key,
                "db": h.get("db") or h.get("source"),
                "id": h.get("id") or h.get("identifier") or h.get("accession"),
                "lib_name": h.get("name"),
                "adduct": h.get("adduct"),
                "lib_precursor_mz": h.get("precursor_mz") or h.get("precursor"),
                "entropy_similarity": h.get("entropy_similarity") or h.get("score") or h.get("similarity"),
                "library_type": h.get("library_type"),
                "lib_rt": h.get("rt") or h.get("retention_time"),
                "ri": h.get("ri") or h.get("retention_index"),
                "rank": h.get("rank") or i,
                "smiles": h.get("smiles"),
                "predicted_rt_hilic": h.get("predicted_rt_hilic"),
                "predicted_rt_rp": h.get("predicted_rt_rp"),
                "anno_delta_rt": h.get("anno_delta_rt"),
                "delta_predicted_rt": h.get("delta_predicted_rt"),
            })

df_hits = pd.DataFrame(rows)
df_hits.to_csv(OUT_HITS, index=False)
print(f'Saved {len(df_hits)} hit rows ({df_hits["wiki_id"].nunique()} spectra) → {OUT_HITS}')

if errors:
    pd.DataFrame(errors).to_csv(OUT_ERRORS, index=False)
    print(f'Saved {len(errors)} errors → {OUT_ERRORS}')

# Summary
n_with_hits = sum(1 for v in cache.values() if v and (v.get("ref") or v.get("anno")))
n_no_hits = sum(1 for v in cache.values() if not v or (not v.get("ref") and not v.get("anno")))
print(f'\nSummary: {n_with_hits} spectra with hits, {n_no_hits} with no hits')
