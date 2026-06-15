"""Fetch library peaks for all hits passing sim >= 0.75 filter."""
import pandas as pd, json, time, requests, warnings
warnings.filterwarnings('ignore')

HITS_PATH    = 'data/orbitrap_hits_refetched.csv'
LIB_CACHE    = 'data/library_peaks_cache.json'
LIB_URL      = 'https://masswiki.us-west-2.elasticbeanstalk.com/reference_library/get_spectra_data'
BATCH_SIZE   = 50
EXCLUDED_DBS = {'5min_hilic_neg', '5min_lipid_neg', 'MB-EU'}
MS2_MIN      = 0.75

spectra = pd.read_excel('data/masswiki_Orbitrap HILIC negESI_2026-03-19.xlsx', header=4)
spectra['label'] = 'unlabeled'
spectra.loc[:1297, 'label'] = 'TP'
spectra.loc[spectra['name'].str.startswith('yy_', na=False), 'label'] = 'FP'
ann = spectra[spectra['label'].isin(['TP', 'FP'])][['wiki_id']].copy()

hits = pd.read_csv(HITS_PATH)
joint = hits.merge(ann, on='wiki_id', how='inner')
joint_f = joint[
    ~joint['db'].isin(EXCLUDED_DBS) &
    (joint['entropy_similarity'] >= MS2_MIN) &
    joint['library_wiki_id'].notna()
]

lib_cache = json.load(open(LIB_CACHE))
unique_lwids = joint_f['library_wiki_id'].unique().tolist()
to_fetch = [lid for lid in unique_lwids if lid not in lib_cache]
print(f'Need to fetch: {len(to_fetch):,}  (already cached: {len(unique_lwids)-len(to_fetch):,})')

errors = 0
for i in range(0, len(to_fetch), BATCH_SIZE):
    batch = to_fetch[i:i + BATCH_SIZE]
    try:
        resp = requests.post(
            LIB_URL,
            json={'id_list': batch, 'get_details': False, 'include_fields': ['peaks']},
            headers={'accept': 'application/json', 'Content-Type': 'application/json'},
            timeout=30
        )
        if resp.status_code == 200:
            for entry in resp.json():
                wid   = entry.get('wiki_id')
                peaks = entry.get('peaks', [])
                if wid and peaks:
                    lib_cache[wid] = peaks
                elif wid:
                    errors += 1
    except Exception as e:
        errors += 1
    if (i // BATCH_SIZE) % 50 == 0 and i > 0:
        print(f'  {i:,}/{len(to_fetch):,}  cache={len(lib_cache):,}')
        json.dump(lib_cache, open(LIB_CACHE, 'w'))
    time.sleep(0.1)

json.dump(lib_cache, open(LIB_CACHE, 'w'))
covered = sum(1 for lid in unique_lwids if lid in lib_cache)
print(f'Done. Cache: {len(lib_cache):,}. Coverage: {covered:,}/{len(unique_lwids):,}. Errors: {errors}')
