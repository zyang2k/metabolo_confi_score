import os, json, re, time
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import ms_entropy
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
RDLogger.DisableLog('rdApp.*')
import matplotlib.pyplot as plt

ROOT = '/Users/ellayoung/Desktop/metabolo_confi_score'

# Neg Orbitrap (training data)
NEG_XLSX       = f'{ROOT}/data/masswiki_Orbitrap HILIC negESI_2026-03-19.xlsx'
NEG_HITS_CSV   = f'{ROOT}/data/orbitrap_hits_refetched.csv'
QUERY_PEAKS    = f'{ROOT}/data/query_peaks_cache.json'

# TTOF
TTOF_NEG_CSV   = f'{ROOT}/data/TTOF_HILIC_negESI_uncurated_041326.csv'
TTOF_POS_CSV   = f'{ROOT}/data/TTOF_HILIC_posESI_uncurated_041326.csv'
TTOF_NEG_HITS  = f'{ROOT}/data/library_hits/ttof_hilic_neg_masswiki_hits.csv'
TTOF_POS_HITS  = f'{ROOT}/data/library_hits/ttof_hilic_pos_masswiki_hits.csv'
TTOF_NEG_PEAKS = f'{ROOT}/data/ttof_neg_query_peaks_cache.json'
TTOF_POS_PEAKS = f'{ROOT}/data/ttof_pos_query_peaks_cache.json'

# Shared
LIB_PEAKS      = f'{ROOT}/data/library_peaks_cache.json'
ADDUCT_TAX     = f'{ROOT}/data/adduct_taxonomy_oliver.csv'
SOLID_TP_PATH  = f'{ROOT}/data/solid_tp.csv'
INCHIKEY_CACHE = f'{ROOT}/data/inchikey_cache.json'
OUT_DIR        = f'{ROOT}/results/ttof_gbm'
os.makedirs(OUT_DIR, exist_ok=True)

PPM_TOL = 10.0

ik_cache = {}
if os.path.exists(INCHIKEY_CACHE):
    with open(INCHIKEY_CACHE) as f:
        ik_cache = json.load(f)

def get_ik14(smiles):
    if not isinstance(smiles, str) or not smiles.strip(): return ''
    if smiles in ik_cache: return ik_cache[smiles].get('ik14', '')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return ''
    ik = InchiToInchiKey(MolToInchi(mol))
    if ik:
        ik_cache[smiles] = {'ik14': ik[:14], 'inchikey': ik}
        return ik[:14]
    return ''

print('Setup complete')

adduct_tax = pd.read_csv(ADDUCT_TAX)
adduct_lookup = dict(zip(adduct_tax['adduct'].str.strip(), adduct_tax['category']))

# Overrides: M+H/M-H alphabetical artifact + pos-mode ok adducts
for a in ['M+H', 'M-H', 'M+Na', 'M+NH4', 'M+K', '2M+H', '2M+Na', '2M+K', '2M+NH4']:
    adduct_lookup[a] = 'ok'
print(f'Adduct overrides applied')

def norm_adduct(s):
    if not isinstance(s, str): return ''
    s = s.strip()
    s = re.sub(r'^\[', '', s)
    s = re.sub(r'\][\+\-]?\d*[\+\-]?$', '', s)
    s = re.sub(r'[\+\-]$', '', s)
    return s.strip()

def classify_adduct(a):
    norm = norm_adduct(a)
    cat = adduct_lookup.get(norm)
    if cat: return cat
    # Rule-based ISF: M+H-X (pos) or M-H-X (neg) = neutral loss
    if re.match(r'M\+H-', norm) or re.match(r'M-H-', norm):
        return 'isf'
    return 'unknown'

def compute_scores(q_peaks, l_peaks, ppm_tol=PPM_TOL):
    if not q_peaks or not l_peaks: return None
    try:
        q_c = ms_entropy.clean_spectrum(q_peaks)
        l_c = ms_entropy.clean_spectrum(l_peaks)
        q_w = ms_entropy.apply_weight_to_intensity(q_c)
        l_w = ms_entropy.apply_weight_to_intensity(l_c)
        q_a = np.array(q_w, dtype=float)
        l_a = np.array(l_w, dtype=float)
        if len(q_a)==0 or len(l_a)==0: return None
        q_mz,q_i = q_a[:,0], q_a[:,1]
        l_mz,l_i = l_a[:,0], l_a[:,1]
        q_n = q_i/q_i.sum(); l_n = l_i/l_i.sum()
        mL=mQ=0.0; pairs=[]; used_q = np.zeros(len(q_mz), bool)
        for j,(lm,li) in enumerate(zip(l_mz,l_n)):
            tol = lm * ppm_tol / 1e6
            d = np.abs(q_mz - lm)
            cand = np.where((d<=tol) & ~used_q)[0]
            if len(cand):
                b = cand[np.argmin(d[cand])]
                mL += li; mQ += q_n[b]
                pairs.append((q_i[b], l_i[j]))
                used_q[b] = True
        res = {'reverse_score':mL, 'forward_score':mQ, 'n_matched':len(pairs)}
        if len(pairs) >= 2:
            qa = np.array([p[0] for p in pairs])
            la = np.array([p[1] for p in pairs])
            res['max_deviation'] = np.max(np.abs(np.log2((qa+1e-6)/(la+1e-6))))
        else:
            res['max_deviation'] = np.nan
        return res
    except Exception:
        return None

def apply_adduct_features(df):
    df = df.copy()
    df['adduct_cat'] = df['adduct'].apply(classify_adduct)
    df['is_isf_adduct'] = (df['adduct_cat']=='isf').astype(int)
    df['is_dubious_adduct'] = (df['adduct_cat']=='dubious').astype(int)
    df['name_lower'] = df['name'].fillna('').str.strip().str.lower()
    has_ok = (df[df['adduct_cat']=='ok'].groupby('name_lower')['wiki_id'].count().rename('_has_ok'))
    df = df.merge(has_ok.reset_index(), on='name_lower', how='left')
    df['has_ok_adduct'] = df['_has_ok'].fillna(0).clip(upper=1).astype(int)
    df = df.drop(columns='_has_ok')
    df['isf_no_mh'] = ((df['is_isf_adduct']==1) & (df['has_ok_adduct']==0)).astype(int)
    n_add = df.groupby('name_lower')['adduct'].nunique().rename('n_compound_adducts')
    df = df.merge(n_add.reset_index(), on='name_lower', how='left')
    df['is_nacetyl'] = df['name'].fillna('').str.contains(
        r'(?i)^N-?acetyl|^N\d-acetyl|^acetyl-.*(?:amine|alanine|valine|leucine|'
        r'isoleucine|glycine|serine|threonine|cysteine|methionine|phenylalanine|'
        r'tyrosine|tryptophan|aspart|glutam|histid|lysine|arginine|proline|ornithine|'
        r'citrulline|carnosine)'
    ).astype(int)
    return df

def build_anno_features(spectra_df, hits_df, query_peaks_cache, lib_peaks_cache):
    hits_df = hits_df.copy()
    hits_df['hit_ik14'] = hits_df['smiles'].fillna('').apply(get_ik14)
    spectra_df['anno_ik14'] = spectra_df['annotation-smiles'].fillna('').apply(get_ik14)
    named_info = spectra_df[['wiki_id','anno_ik14','label','tier']].copy()
    joint = hits_df[hits_df['hit_source']=='reference'].merge(named_info, on='wiki_id', how='inner')
    joint['entropy_similarity'] = pd.to_numeric(joint['entropy_similarity'], errors='coerce')
    joint['is_anno_hit'] = (joint['hit_ik14']!='') & (joint['anno_ik14']!='') & (joint['hit_ik14']==joint['anno_ik14'])
    joint['dedup_key'] = joint.apply(
        lambda r: (r['wiki_id'], r['hit_ik14']) if r['hit_ik14'] else (r['wiki_id'], r.get('lib_name', r.get('name', ''))), axis=1)
    joint = joint.sort_values('entropy_similarity', ascending=False).drop_duplicates('dedup_key').reset_index(drop=True)
    rows = []
    for wid, g in joint.groupby('wiki_id'):
        anno = g[g['is_anno_hit']]
        others = g[~g['is_anno_hit']]
        if len(anno) > 0:
            best = anno.loc[anno['entropy_similarity'].idxmax()]
            next_sim = others['entropy_similarity'].max() if len(others) > 0 else 0.0
            anno_rank = int((g['entropy_similarity'] >= best['entropy_similarity']).sum())
            lib_pk = lib_peaks_cache.get(best.get('library_wiki_id',''))
            qry_pk = query_peaks_cache.get(wid)
            sc = compute_scores(qry_pk, lib_pk)
            lib_mz = pd.to_numeric(best.get('lib_precursor_mz'), errors='coerce')
            obs_mz = spectra_df.loc[spectra_df['wiki_id']==wid,'precursor_mz'].iloc[0]
            row = {'wiki_id': wid,
                   'anno_entropy_sim': best['entropy_similarity'],
                   'anno_delta_mda': abs(obs_mz-lib_mz)*1000 if pd.notna(lib_mz) and lib_mz>0 else np.nan,
                   'sim_gap': max(0.0, best['entropy_similarity']-next_sim) if pd.notna(next_sim) else best['entropy_similarity'],
                   'anno_rank': anno_rank, 'has_anno_hit': True}
            if sc:
                row['anno_reverse']=sc['reverse_score']; row['anno_forward']=sc['forward_score']; row['max_deviation']=sc['max_deviation']
            else:
                row['anno_reverse']=np.nan; row['anno_forward']=np.nan; row['max_deviation']=np.nan
        else:
            row = dict.fromkeys(['anno_entropy_sim','anno_delta_mda','sim_gap','anno_rank',
                                  'anno_reverse','anno_forward','max_deviation'], np.nan)
            row['wiki_id']=wid; row['has_anno_hit']=False
        rows.append(row)
    return pd.DataFrame(rows)

def build_evidence(spectra_df, anno_df):
    ev = spectra_df[['wiki_id','name','label','tier','precursor_mz',
                      'is_isf_adduct','is_dubious_adduct','has_ok_adduct','isf_no_mh',
                      'n_compound_adducts','is_nacetyl','adduct']].copy()
    ev['spectral_entropy'] = pd.to_numeric(spectra_df['entropy'], errors='coerce')
    ev['delta_rt_abs'] = pd.to_numeric(spectra_df['anno_delta_rt'], errors='coerce').abs()
    ev = ev.merge(anno_df, on='wiki_id', how='left')
    return ev

print('All helpers ready')


# Load + label neg
neg_raw = pd.read_excel(NEG_XLSX, header=4)
neg_raw['label'] = 'unlabeled'
neg_raw.loc[:1297, 'label'] = 'TP'
neg_raw.loc[neg_raw['name'].str.startswith('yy_', na=False), 'label'] = 'FP'
neg_raw.loc[neg_raw['name'].str.startswith('zz_', na=False), 'label'] = 'TN'
solid = pd.read_csv(SOLID_TP_PATH)
neg_raw['tier'] = 'first_pass'
neg_raw.loc[:1297, 'tier'] = 'regular_tp'
neg_raw.loc[neg_raw['wiki_id'].isin(set(solid['wiki_id'])), 'tier'] = 'golden_tp'
neg_raw.loc[neg_raw['label']=='FP', 'tier'] = 'fp'
neg = neg_raw[neg_raw['label'].isin(['TP','FP'])].copy().reset_index(drop=True)
neg = apply_adduct_features(neg)

# Load neg hits + peaks
neg_hits = pd.read_csv(NEG_HITS_CSV, low_memory=False)
neg_hits = neg_hits.rename(columns={'id':'library_id','lib_name':'name'})
with open(QUERY_PEAKS) as f:
    neg_query_peaks = json.load(f)
with open(LIB_PEAKS) as f:
    lib_peaks_cache = json.load(f)
print(f'Neg: {len(neg)} spectra, {len(neg_hits):,} hits, {len(neg_query_peaks)} query peaks, {len(lib_peaks_cache):,} lib peaks')

# Compute neg features
print('Computing neg features...')
t0 = time.time()
neg_anno = build_anno_features(neg, neg_hits, neg_query_peaks, lib_peaks_cache)
neg_ev = build_evidence(neg, neg_anno)
print(f'  {time.time()-t0:.0f}s')

# Save IK cache
with open(INCHIKEY_CACHE, 'w') as f:
    json.dump(ik_cache, f)

FEATURE_COLS = [
    'delta_rt_abs','spectral_entropy','anno_entropy_sim','anno_forward','anno_reverse',
    'anno_delta_mda','sim_gap','anno_rank','max_deviation',
    'is_isf_adduct','isf_no_mh','is_dubious_adduct','n_compound_adducts','is_nacetyl',
]

train_mask = neg_ev['tier'].isin(['regular_tp','fp'])
df_train = neg_ev[train_mask].copy()
y_train = (df_train['tier']=='regular_tp').astype(int)
train_medians = df_train[FEATURE_COLS].median()
X_train = df_train[FEATURE_COLS].fillna(train_medians).values

gb = GradientBoostingClassifier(n_estimators=100, max_depth=3, learning_rate=0.1, random_state=42)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_probs = cross_val_predict(gb, X_train, y_train, cv=skf, method='predict_proba')[:,1]
cv_auc = roc_auc_score(y_train, cv_probs)
gb.fit(X_train, y_train)

print(f'Neg GBM: train={len(df_train)} (TP={int(y_train.sum())}, FP={int((~y_train.astype(bool)).sum())})')
print(f'5-fold CV AUC: {cv_auc:.3f}')
print('\nFeature importance:')
for col, imp in sorted(zip(FEATURE_COLS, gb.feature_importances_), key=lambda x: -x[1]):
    print(f'  {col:25s}  {imp:.3f}')

results = {}

for mode, csv_path, hits_path, peaks_path in [
    ('neg', TTOF_NEG_CSV, TTOF_NEG_HITS, TTOF_NEG_PEAKS),
    ('pos', TTOF_POS_CSV, TTOF_POS_HITS, TTOF_POS_PEAKS),
]:
    print(f'\n{"="*60}')
    print(f'TTOF {mode.upper()}')
    print(f'{"="*60}')

    # Check data exists
    if not os.path.exists(hits_path):
        print(f'  Hits file not found: {hits_path}')
        print(f'  Run: python code/fetch_ttof_hits_and_peaks.py <TOKEN>')
        continue

    # Load + label
    raw = pd.read_csv(csv_path, low_memory=False)
    names = raw['name'].fillna('').astype(str)
    raw['label'] = 'unlabeled'
    raw.loc[raw['is_manual_annotated'].astype(bool) & ~names.str.startswith('yy_') & ~names.str.startswith('zz_') & (names.str.strip()!=''), 'label'] = 'TP'
    raw.loc[names.str.startswith('yy_'), 'label'] = 'FP'
    raw.loc[names.str.startswith('zz_'), 'label'] = 'TN'
    raw['tier'] = 'unlabeled'
    raw.loc[raw['label']=='TP', 'tier'] = 'regular_tp'
    raw.loc[raw['label']=='FP', 'tier'] = 'fp'
    df = raw[raw['label'].isin(['TP','FP'])].copy().reset_index(drop=True)
    df = apply_adduct_features(df)
    print(f'  Annotated: {len(df)} (TP={df["label"].eq("TP").sum()}, FP={df["label"].eq("FP").sum()})')

    # Load hits + peaks
    hits = pd.read_csv(hits_path, low_memory=False)
    hits = hits.rename(columns={'id':'library_id','lib_name':'name'})
    with open(peaks_path) as f:
        qpeaks = json.load(f)
    print(f'  Hits: {len(hits):,}, Query peaks: {len(qpeaks):,}')

    # Fetch missing library peaks (no auth)
    lwids = hits.get('library_wiki_id', pd.Series()).dropna().unique().tolist()
    missing_lib = [l for l in lwids if l not in lib_peaks_cache]
    if missing_lib:
        import requests
        LIB_URL = 'https://masswiki.us-west-2.elasticbeanstalk.com/reference_library/get_spectra_data'
        print(f'  Fetching {len(missing_lib)} library peaks...')
        for i in range(0, len(missing_lib), 50):
            batch = missing_lib[i:i+50]
            try:
                r = requests.post(LIB_URL, json={'id_list': batch, 'get_details': False, 'include_fields': ['peaks']}, timeout=30)
                if r.status_code == 200:
                    for e in r.json():
                        wid = e.get('wiki_id'); pk = e.get('peaks', [])
                        if wid and pk: lib_peaks_cache[wid] = pk
            except: pass
            time.sleep(0.1)
        with open(LIB_PEAKS, 'w') as f:
            json.dump(lib_peaks_cache, f)

    # Compute features
    print(f'  Computing features...')
    t0 = time.time()
    anno_df = build_anno_features(df, hits, qpeaks, lib_peaks_cache)
    ev = build_evidence(df, anno_df)
    elapsed = time.time() - t0
    print(f'  Feature compute: {elapsed:.0f}s ({len(ev)} spectra)')

    # Coverage
    for ch in ['delta_rt_abs','anno_entropy_sim','anno_forward','sim_gap','max_deviation']:
        n = ev[ch].notna().sum()
        print(f'    {ch:20s}  {n}/{len(ev)} ({n/len(ev):.0%})')

    # Apply GBM
    X = ev[FEATURE_COLS].fillna(train_medians).values
    y = (ev['tier']=='regular_tp').astype(int)
    scores = gb.predict_proba(X)[:,1]
    ev['confidence_score'] = scores
    auc = roc_auc_score(y, scores) if y.nunique() > 1 else float('nan')

    tp_s = scores[y.astype(bool)]
    fp_s = scores[~y.astype(bool)]
    print(f'\n  AUC: {auc:.3f}  (neg ref: {cv_auc:.3f}, delta: {auc-cv_auc:+.3f})')
    print(f'  TP scores: mean={tp_s.mean():.3f}, median={np.median(tp_s):.3f}')
    print(f'  FP scores: mean={fp_s.mean():.3f}, median={np.median(fp_s):.3f}')

    # Export
    out = df.merge(
        ev[['wiki_id','confidence_score','has_anno_hit','anno_entropy_sim','anno_forward',
            'anno_reverse','sim_gap','anno_rank','max_deviation','anno_delta_mda',
            'is_isf_adduct','isf_no_mh','is_dubious_adduct','n_compound_adducts','is_nacetyl']],
        on='wiki_id', how='left'
    ).sort_values('confidence_score', ascending=False, na_position='last')

    xlsx = f'{OUT_DIR}/ttof_{mode}_gbm_scores.xlsx'
    out.to_excel(xlsx, index=False)
    out.to_csv(f'{OUT_DIR}/ttof_{mode}_gbm_scores.csv', index=False)
    print(f'  Saved: {xlsx} ({len(out)} rows)')

    results[mode] = {'auc': auc, 'n': len(ev), 'tp_mean': tp_s.mean(), 'fp_mean': fp_s.mean()}

# Save IK cache
with open(INCHIKEY_CACHE, 'w') as f:
    json.dump(ik_cache, f)

print('\n' + '='*60)
print('CROSS-PLATFORM GENERALIZATION SUMMARY')
print('='*60)
print(f'{"Dataset":<25s} {"AUC":>6s}  {"n":>5s}  {"TP mean":>8s}  {"FP mean":>8s}')
print('-'*60)
print(f'{"Orbitrap neg (CV)":<25s} {cv_auc:>6.3f}  {len(df_train):>5d}  {"—":>8s}  {"—":>8s}')
for mode, r in results.items():
    print(f'{"TTOF " + mode:<25s} {r["auc"]:>6.3f}  {r["n"]:>5d}  {r["tp_mean"]:>8.3f}  {r["fp_mean"]:>8.3f}')
print('-'*60)
print(f'\nOutput files:')
for mode in results:
    print(f'  results/ttof_gbm/ttof_{mode}_gbm_scores.xlsx')

print(f'\nInterpretation:')
for mode, r in results.items():
    delta = r['auc'] - cv_auc
    if abs(delta) < 0.05:
        verdict = 'GENERALIZES (within 0.05 of reference)'
    elif delta < -0.05:
        verdict = f'DEGRADED ({delta:+.3f} — platform shift hurts)'
    else:
        verdict = f'IMPROVED ({delta:+.3f})'
    print(f'  TTOF {mode}: {verdict}')