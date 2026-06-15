import pandas as pd
import glob
import os

# 1. Point to the folder with your CSV files
folder_path = "/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof_pos_rp_wiki_match"

# 2. Get list of CSV files
csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

# 3. Read and combine CSVs
dfs = [pd.read_csv(f) for f in csv_files]
combined_df = pd.concat(dfs, ignore_index=True)

# 4. Drop duplicate rows across all columns
combined_df = combined_df.drop_duplicates(keep='first')

# 5. Reset index for cleanliness
combined_df.reset_index(drop=True, inplace=True)

spectrum = pd.read_csv("/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof+pos+rp.csv")

# 0) Import the MVP modules I dropped in /mnt/data/mvp
import sys
sys.path.append('/Users/ellayoung/Desktop/metabolo_confi_score/code/0918mvp')

from mvp_ingest import join_hits_to_spectra, build_confusable_sets
from mvp_score import fit_channel_params, per_hit_logLR, aggregate_to_structure, local_rank_probability, compute_global_posterior
from mvp_export import make_assertion_table, make_top_calls

# 1) Normalize hits in memory (no file I/O)
hits = (
    combined_df
      .rename(columns={
          'id':'library_id',
          'precursor_mz':'lib_precursor_mz',
          'rt':'lib_rt'              # will be NaN in your sample, that’s fine
      })[
        ['wiki_id','library_id','name','adduct','lib_precursor_mz','lib_rt','entropy_similarity',
         'db','library_type','rank']
      ]
)

# 2) Minimal spectra view (fill optional fields if missing)
spec = spectrum[['wiki_id','rt','precursor_mz','is_manual_annotated']].copy()
spec['assay'] = 'unknown'
spec['polarity'] = 'unknown'

# 3) Join + confusable set (MS2 filter; RT passes through if lib_rt present)
joint = join_hits_to_spectra(hits, spec)
cset  = build_confusable_sets(joint, ms2_min=0.75, zrt_k=3.0, rt_sigma_sec=2.0)

# After you build `cset` (and before aggregate_to_structure)
# 1) normalize keys so NaNs don’t explode groups
cset['name']   = cset['name'].fillna('UNKNOWN')
cset['adduct'] = cset['adduct'].fillna('NA')

# 2) define the structure alias as NAME||ADDUCT
cset['library_id'] = cset['name'].astype(str) + '||' + cset['adduct'].astype(str)

# (optional but nice) keep a clean display name
cset['display_name'] = cset['name']  # you can use this in exports

# 4) Fit quick channel params (labels from is_manual_annotated) → per-hit logLR
# --- B) Clean hit-level labels using spectrum annotations ---
# 1) Bring in the annotation fields with controlled suffixes
ann = spectrum[['wiki_id','annotation-name','annotation-adduct']].copy()
lab = cset.merge(ann, on='wiki_id', how='left', suffixes=('', '_ann'))

# 2) Normalize fields for robust equality
lab['name_norm']       = lab['name'].fillna('').str.strip().str.lower()
lab['ann_name_norm']   = lab['annotation-name'].fillna('').str.strip().str.lower()
lab['adduct_norm']     = lab['adduct'].fillna('').str.replace(r'\s+', '', regex=True).str.upper()
lab['ann_adduct_norm'] = lab['annotation-adduct'].fillna('').str.replace(r'\s+', '', regex=True).str.upper()

# 3) Decide whether this spectrum is annotated
#    Prefer cset's column if present; otherwise infer from presence of annotation fields
if 'is_manual_annotated' in lab.columns:
    lab['is_annotated'] = lab['is_manual_annotated'].fillna(False)
elif 'is_manual_annotated' in spectrum.columns:
    lab = lab.merge(spectrum[['wiki_id','is_manual_annotated']], on='wiki_id', how='left', suffixes=('', '_spec'))
    lab['is_annotated'] = lab['is_manual_annotated'].fillna(lab['is_manual_annotated_spec']).fillna(False)
else:
    lab['is_annotated'] = lab['annotation-name'].notna() & lab['annotation-adduct'].notna()

# 4) Hit-level label: TP only if (a) spectrum is annotated AND (b) hit matches annotated name+adduct
lab['label'] = (
    lab['is_annotated'] &
    (lab['name_norm']   == lab['ann_name_norm']) &
    (lab['adduct_norm'] == lab['ann_adduct_norm'])
).astype(int)

# 5) Fit channel params using these clean labels
params = fit_channel_params(lab)

params["ppm"] = {
    "model": "t",     # TP ~ Student-t for Δppm
    "df": 3.0,
    "scale": 2.7017,  # from your diag
    "b": 3.3926       # FP Laplace b from your diag
}



cset['logLR_hit'] = per_hit_logLR(cset, params)

# --- sibling-import path shim (code/ root) ---
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from diagnostics_channels import diag_delta_ppm, diag_ms2
res_ppm = diag_delta_ppm(lab, 'label', 'delta_ppm'); print(res_ppm)
res_ms2 = diag_ms2(lab, 'label', 'entropy_similarity'); print(res_ms2)


# 5) Aggregate → structure score θ, local rank C_loc, global posterior with NoTA
struct  = aggregate_to_structure(cset, quality_cols=('entropy_similarity',), topK=None, lambda_var=0.0)
# right after: struct = aggregate_to_structure(...)
struct = struct.sort_values(['wiki_id','library_id','theta'],
                            ascending=[True, True, False]) \
               .drop_duplicates(['wiki_id','library_id'], keep='first')

struct  = local_rank_probability(struct)
post    = compute_global_posterior(struct, prior_policy='uniform', lr_nota=0.5, prior_nota=0.5)

# 6) Assertions + top calls
assertions = make_assertion_table(post, features_df=cset)
top_calls  = make_top_calls(assertions, nota_threshold=0.6)

# Peek
assertions.head(), top_calls.head()
print("Assertions,", assertions.shape)

# coverage: how many spectra now have at least one candidate?
n_spec = spectrum['wiki_id'].nunique()
n_cov  = assertions['wiki_id'].nunique()
print(f"Spectra total: {n_spec} | with candidates: {n_cov} ({n_cov/n_spec:.1%})")

# candidates per spectrum distribution
cand_per_spec = assertions.groupby('wiki_id').size()
print(cand_per_spec.describe())

# abstentions
print("top_calls shape:", top_calls.shape)
print("abstain %:", top_calls['abstain'].mean(), " count:", top_calls['abstain'].sum())

# sanity: unique rows = unique (wiki_id, library_id)
assert assertions[['wiki_id','library_id']].drop_duplicates().shape[0] == len(assertions)

# peek top confident candidates
cols = ['wiki_id','name','theta','C_loc','PEP_local','entropy_similarity','delta_ppm']
print(assertions.sort_values('theta', ascending=False)[cols].head(10))

# how many hits had missing lib_rt in the scoring stage?
print("RT missing rate in hits:", cset['rt_missing'].mean())

import numpy as np
import pandas as pd

# 1) Ensure per-hit LR exists
cset = cset.copy()
cset['LR_hit'] = np.exp(cset['logLR_hit'].astype(float))

# 2) Use a quality weight (same as aggregation; here entropy sim)
cset['_w'] = cset['entropy_similarity'].astype(float).clip(lower=1e-6).fillna(1e-6)

# 3) Conditional prob of a hit given its structure: P(h | structure)
cset['_w_norm'] = cset['_w'] / cset.groupby(['wiki_id','library_id'])['_w'].transform('sum')
cset['_num'] = cset['_w_norm'] * cset['LR_hit']
cset['_den'] = cset.groupby(['wiki_id','library_id'])['_num'].transform('sum').replace(0, np.nan)
cset['P_h_given_struct'] = cset['_num'] / cset['_den']

# 4) Bring in structure-level probabilities
#    - local (within confusable set): C_loc from `struct`
#    - global (with NoTA): post from `post`
cols_struct = ['wiki_id','library_id','C_loc']
if 'C_loc' not in struct.columns:  # safety
    raise RuntimeError("Run local_rank_probability(struct) first.")
hit_scores = cset.merge(struct[cols_struct], on=['wiki_id','library_id'], how='left')

cols_post = ['wiki_id','library_id','post','P_NoTA']
if not set(cols_post).issubset(post.columns):
    raise RuntimeError("Run compute_global_posterior(struct) first.")
hit_scores = hit_scores.merge(post[cols_post].drop_duplicates(), on=['wiki_id','library_id'], how='left')

# 5) Per-hit probabilities
#    Local (excludes NoTA): share the structure’s local mass to its hits
hit_scores['match_post_local']  = hit_scores['C_loc'] * hit_scores['P_h_given_struct']

#    Global (includes NoTA): share the global structure posterior to hits
hit_scores['match_post_global'] = hit_scores['post']  * hit_scores['P_h_given_struct']
hit_scores['PEP_hit']           = 1.0 - hit_scores['match_post_global']  # smaller is better

# 6) (Optional) naive per-hit softmax across ALL hits of a spectrum (duplicates inflate!)
hit_scores['C_loc_hit_naive'] = (
    hit_scores.groupby('wiki_id')['logLR_hit']
    .transform(lambda x: np.exp(x - x.max()) / np.exp(x - x.max()).sum())
)

# 7) Final columns to keep
hit_scores_out = hit_scores[[
    'wiki_id','library_id','name','adduct','entropy_similarity','delta_ppm',
    'logLR_hit','LR_hit','P_h_given_struct','C_loc','post','P_NoTA',
    'match_post_local','match_post_global','PEP_hit','C_loc_hit_naive'
]].sort_values(['wiki_id','match_post_global'], ascending=[True, False])

print("Hit-level scores:", hit_scores_out.head(10))

import os

OUTDIR = "/Users/ellayoung/Desktop/metabolo_confi_score/out_0918"
os.makedirs(OUTDIR, exist_ok=True)

assertions.to_csv(f"{OUTDIR}/assertions.csv", index=False)
top_calls.to_csv(f"{OUTDIR}/top_calls.csv", index=False)
hit_scores_out.to_csv(f"{OUTDIR}/hit_scores.csv", index=False)  # per-match scores

print("Wrote:",
      f"{OUTDIR}/assertions.csv",
      f"{OUTDIR}/top_calls.csv",
      f"{OUTDIR}/hit_scores.csv",
      sep="\n")


# After running your scoring pipeline:
# cset has: delta_ppm, entropy_similarity, zrt (NaN if no RT), is_manual_annotated
from diagnostics_channels import diag_delta_ppm, diag_ms2, diag_zrt

labeled = cset.assign(label=cset['is_manual_annotated'].astype(int))

res_ppm = diag_delta_ppm(labeled, label_col='label', value_col='delta_ppm', title="Δppm")
print(res_ppm)  # shows σ for TP, b for FP, and AICs

res_ms2 = diag_ms2(labeled, label_col='label', value_col='entropy_similarity', title="MS2 entropy")
print(res_ms2)  # shows μ/σ for TP/FP and AICs

# Only run if you actually have zRT (lib_rt present)
if 'zrt' in labeled.columns and np.isfinite(labeled['zrt']).any():
    res_zrt = diag_zrt(labeled, label_col='label', value_col='zrt', title="zRT")
    print(res_zrt)
