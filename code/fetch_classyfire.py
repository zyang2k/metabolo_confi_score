"""
fetch_classyfire.py
-------------------
Classify confirmed lib_names using NPClassifier (UCSD/GNPS).
ClassyFire server is down; NPClassifier works directly from SMILES.

NPClassifier hierarchy:
  pathway   (broadest)  ~= ClassyFire kingdom/superclass
  superclass             ~= ClassyFire class
  class      (finest)   ~= ClassyFire subclass

Output: results/current/lib_classyfire.csv
  lib_name, smiles, npc_pathway, npc_superclass, npc_class
"""
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote

import pandas as pd
import requests

ROOT  = '/Users/ellayoung/Desktop/metabolo_confi_score'
DATA  = f'{ROOT}/data'
OUT   = f'{ROOT}/results/current'
CACHE = f'{OUT}/classyfire_cache.json'

NPCLASSIFIER_URL = 'https://npclassifier.ucsd.edu/classify?smiles={}'
MAX_WORKERS      = 5
SLEEP_PER_REQ    = 0.3   # ~3 RPS per worker → ~15 RPS total

# ── 1. Load data ──────────────────────────────────────────────────────────
print('Loading data...')
assertions = pd.read_csv(f'{OUT}/assertions.csv')
hits       = pd.read_csv(f'{DATA}/hilic_ttof_neg_masswiki_hits_new.csv', low_memory=False)

confirmed   = assertions[assertions['confirmed']].copy()
all_names   = [n for n in confirmed['lib_name'].dropna().unique()
               if not str(n).startswith('yy')]

name_to_smiles = (hits.dropna(subset=['smiles'])
                      .groupby('lib_name')['smiles']
                      .first()
                      .to_dict())

# Deduplicate by SMILES — query each unique SMILES once
smiles_to_names = {}
for name in all_names:
    smi = name_to_smiles.get(name)
    if smi:
        smiles_to_names.setdefault(smi, []).append(name)

unique_smiles = list(smiles_to_names.keys())
print(f'lib_names: {len(all_names)}  |  unique SMILES to classify: {len(unique_smiles)}')

# ── 2. Load cache ─────────────────────────────────────────────────────────
if os.path.exists(CACHE):
    with open(CACHE) as f:
        cache = json.load(f)
    print(f'Cache loaded: {len(cache)} entries')
else:
    cache = {}

todo = [smi for smi in unique_smiles if smi not in cache]
print(f'Remaining: {len(todo)}')

# ── 3. Fetch NPClassifier ─────────────────────────────────────────────────
session = requests.Session()

def fetch_one(smi):
    time.sleep(SLEEP_PER_REQ)
    url = NPCLASSIFIER_URL.format(quote(smi, safe=''))
    try:
        r = session.get(url, timeout=20)
        if r.status_code == 200:
            d = r.json()
            return smi, {
                'npc_pathway':    ', '.join(d.get('pathway_results', [])),
                'npc_superclass': ', '.join(d.get('superclass_results', [])),
                'npc_class':      ', '.join(d.get('class_results', [])),
                'isglycoside':    d.get('isglycoside', False),
            }
        else:
            return smi, None
    except Exception:
        return smi, None

batch_size = 200
failed = []

for batch_start in range(0, len(todo), batch_size):
    batch = todo[batch_start: batch_start + batch_size]
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(fetch_one, smi): smi for smi in batch}
        for fut in as_completed(futures):
            smi, result = fut.result()
            if result is not None:
                cache[smi] = result
            else:
                failed.append(smi)

    with open(CACHE, 'w') as f:
        json.dump(cache, f)

    done = min(batch_start + batch_size, len(todo))
    print(f'  {done}/{len(todo)} done  ({len(failed)} failed)', end='\r')

print(f'\nFetch complete. Cache: {len(cache)}  Failed: {len(failed)}')

# Retry failed
if failed:
    print(f'Retrying {len(failed)}...')
    with ThreadPoolExecutor(max_workers=2) as ex:
        futures = {ex.submit(fetch_one, smi): smi for smi in failed}
        for fut in as_completed(futures):
            smi, result = fut.result()
            if result is not None:
                cache[smi] = result
    with open(CACHE, 'w') as f:
        json.dump(cache, f)

# ── 4. Build output ───────────────────────────────────────────────────────
rows = []
for name in all_names:
    smi = name_to_smiles.get(name, '')
    cf  = cache.get(smi, {}) if smi else {}
    rows.append({
        'lib_name':       name,
        'smiles':         smi,
        'npc_pathway':    cf.get('npc_pathway', ''),
        'npc_superclass': cf.get('npc_superclass', ''),
        'npc_class':      cf.get('npc_class', ''),
        'isglycoside':    cf.get('isglycoside', ''),
    })

cf_df = pd.DataFrame(rows).drop_duplicates('lib_name')
cf_df.to_csv(f'{OUT}/lib_classyfire.csv', index=False)

print(f'\nSaved {len(cf_df)} rows → {OUT}/lib_classyfire.csv')
print()
print('NPC pathway distribution (top 15):')
print(cf_df['npc_pathway'].replace('', 'Unknown').value_counts().head(15))
