# --- 0) Imports & config ---------------------------------------------
import sys, os, numpy as np, pandas as pd
from scipy.stats import norm, laplace, t

import pandas as pd
import glob
import os

# # 1. Point to the folder with your CSV files
# folder_path = "/Users/ellayoung/Desktop/metabolo_confi_score/data/hilic_masswiki_reference_hits"

# # 2. Get list of CSV files
# csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

df = pd.read_csv("/Users/ellayoung/Desktop/metabolo_confi_score/data/hilic_masswiki_reference_hits.csv")
# 4. Drop duplicate rows across all columns
combined_df = df.drop_duplicates(keep='first')

# 5. Reset index for cleanliness
combined_df.reset_index(drop=True, inplace=True)

spectrum = pd.read_csv('/Users/ellayoung/Desktop/metabolo_confi_score/data/ttof+neg+hilic.csv')
# wikimatch_label_fit.py
import os, sys, numpy as np, pandas as pd
from scipy.stats import norm, laplace, t
# wikimatch_analysis_label_fit.py
import os, sys, numpy as np, pandas as pd
from scipy.stats import norm, laplace, t

# --- use your MVP modules exactly like the 0918 template ---
sys.path.append('/Users/ellayoung/Desktop/metabolo_confi_score/code/0918mvp')
from mvp_ingest import join_hits_to_spectra, build_confusable_sets
from mvp_score import aggregate_to_structure, local_rank_probability, compute_global_posterior
from mvp_export import make_assertion_table, make_top_calls

# Assumes combined_df and spectrum are already in memory (as in your 0918 script)
OUTDIR     = "/Users/ellayoung/Desktop/metabolo_confi_score/out_1017"
MS2_MIN    = 0.6          # gate for confusable set
PPM_MAX    = 20.0          # |Δppm| gate (ppm)
PRIOR_NOTA = 0.05          # softer NoTA so strong evidence can win
LR_NOTA    = 1.0

os.makedirs(OUTDIR, exist_ok=True)

# ---------------- 1) Normalize hits & spectra (template-compatible) ----------------
hits = (
    combined_df
      .rename(columns={'id':'library_id','precursor_mz':'lib_precursor_mz','rt':'lib_rt'})[
        ['wiki_id','library_id','name','adduct','lib_precursor_mz','lib_rt',
         'entropy_similarity','db','library_type','rank']
      ].copy()
)

# spec MUST have is_manual_annotated before join (like your 0918 script)
spec = spectrum[['wiki_id','rt','precursor_mz']].copy()
if 'is_manual_annotated' in spectrum.columns:
    spec['is_manual_annotated'] = spectrum['is_manual_annotated'].fillna(False).astype(bool)
else:
    # infer from presence of annotation fields if flag not present
    has_ann = spectrum.get('annotation-name') is not None and spectrum.get('annotation-adduct') is not None
    if has_ann:
        spec = spec.merge(
            spectrum[['wiki_id','annotation-name','annotation-adduct']],
            on='wiki_id', how='left'
        )
        spec['is_manual_annotated'] = (spec['annotation-name'].notna() & spec['annotation-adduct'].notna()).astype(bool)
    else:
        spec['is_manual_annotated'] = False
spec['assay'] = 'unknown'
spec['polarity'] = 'unknown'

# ---------------- 2) Join; compute Δppm; build confusable set ----------------
joint = join_hits_to_spectra(hits, spec)  # same API as your 0918 file

dmz = joint['precursor_mz'] - joint['lib_precursor_mz']
joint['delta_ppm'] = dmz

cset = build_confusable_sets(joint, ms2_min=MS2_MIN, zrt_k=3.0, rt_sigma_sec=2.0)
cset = cset[cset['delta_ppm'].abs() <= PPM_MAX].copy()

# collapse duplicates by (name, adduct) like we discussed
cset['name']   = cset['name'].fillna('UNKNOWN')
cset['adduct'] = cset['adduct'].fillna('NA')
cset['library_id'] = cset['name'].astype(str) + '||' + cset['adduct'].astype(str)

# ---------------- 3) Labeling-for-fit from annotations ----------------
# bring annotation fields (don’t rely on manual flag names)
ann_cols = [c for c in ['annotation-name','annotation-adduct'] if c in spectrum.columns]
lab = cset.merge(spectrum[['wiki_id'] + ann_cols], on='wiki_id', how='left')

def _norm_name(s: pd.Series) -> pd.Series:
    return s.fillna('').str.strip().str.lower()
def _norm_adduct(s: pd.Series) -> pd.Series:
    return s.fillna('').str.replace(r'\s+','', regex=True).str.upper()

lab['name_norm']       = _norm_name(lab['name'])
lab['adduct_norm']     = _norm_adduct(lab['adduct'])
lab['ann_name_norm']   = _norm_name(lab.get('annotation-name', pd.Series(index=lab.index)))
lab['ann_adduct_norm'] = _norm_adduct(lab.get('annotation-adduct', pd.Series(index=lab.index)))

lab['is_annotated'] = lab.get('annotation-name', pd.Series(index=lab.index)).notna() & \
                      lab.get('annotation-adduct', pd.Series(index=lab.index)).notna()

lab['label'] = (
    lab['is_annotated'] &
    (lab['name_norm']   == lab['ann_name_norm']) &
    (lab['adduct_norm'] == lab['ann_adduct_norm'])
).astype(int)

# ---------------- 4) Fit channel params from labels ----------------
# Δppm: TP ~ t(df=3), FP ~ Laplace
df_t   = 3.0
tp_ppm = lab.loc[lab['label']==1, 'delta_ppm'].dropna().astype(float)
fp_ppm = lab.loc[lab['label']==0, 'delta_ppm'].dropna().astype(float)
if len(tp_ppm) >= 5:
    var_tp  = float(tp_ppm.var(ddof=1))
    scale_t = (max(var_tp, 1e-12) * (df_t-2.0)/df_t) ** 0.5
else:
    scale_t = 3.0
b_lap = float(fp_ppm.abs().mean()) if len(fp_ppm) >= 5 else 8.0

# MS2: logit(upper-half) transform, fit Normals with guardrails
def _ms2_transform(arr, eps=1e-6):
    s = np.maximum(np.asarray(arr, float), 0.5)
    z = (s - 0.5) / 0.5
    z = np.clip(z, eps, 1-eps)
    return np.log(z/(1-z))

tp_ms2 = _ms2_transform(lab.loc[lab['label']==1, 'entropy_similarity'].dropna())
fp_ms2 = _ms2_transform(lab.loc[lab['label']==0, 'entropy_similarity'].dropna())
if len(tp_ms2) >= 20 and len(fp_ms2) >= 20:
    mu_tp, sd_tp = float(tp_ms2.mean()), float(max(tp_ms2.std(ddof=1), 0.3))
    mu_fp, sd_fp = float(fp_ms2.mean()), float(max(fp_ms2.std(ddof=1), 0.4))
    if not (mu_tp > mu_fp):  # ensure high sim helps
        mu_tp, sd_tp, mu_fp, sd_fp = 1.2, 0.6, -0.2, 1.0
else:
    mu_tp, sd_tp, mu_fp, sd_fp = 1.2, 0.6, -0.2, 1.0

params = {
    'ppm': {'model':'t', 'df':df_t, 'scale':scale_t, 'b':b_lap},
    'zrt': {'sigma_tp':1.0, 'sigma_fp':2.0},  # placeholder if RT missing
    'ms2': {'mu_tp':mu_tp, 'sd_tp':sd_tp, 'mu_fp':mu_fp, 'sd_fp':sd_fp},
}

# ensure Δppm at 0 favors TP (center density check)
t0   = 2.0 / (params['ppm']['scale'] * (3.0*np.pi)**0.5)
lap0 = 1.0 / (2.0 * params['ppm']['b'])
if lap0 >= t0:
    params['ppm']['b'] = max(params['ppm']['b'], 1.8*params['ppm']['scale'])

# ---------------- 5) Per-hit logLR (inline—no module dependency) ----------------
def _ms2_scalar(s, eps=1e-6):
    s = 0.5 + max(0.0, float(s) - 0.5)
    z = (s - 0.5) / 0.5
    z = min(max(z, eps), 1.0 - eps)
    return np.log(z/(1 - z))

def per_hit_logLR_inline(df: pd.DataFrame, params: dict) -> pd.Series:
    out = np.zeros(len(df), dtype=float)
    # Δppm: TP=t, FP=Laplace
    x = df['delta_ppm'].astype(float).to_numpy()
    out += t.logpdf(x, df=params['ppm']['df'], loc=0.0, scale=params['ppm']['scale']) \
         - laplace.logpdf(x, loc=0.0, scale=params['ppm']['b'])
    # zRT (skip NaNs)
    if 'zrt' in df.columns:
        r = df['zrt'].astype(float).to_numpy()
        add = norm.logpdf(r, 0.0, params['zrt']['sigma_tp']) - norm.logpdf(r, 0.0, params['zrt']['sigma_fp'])
        add[~np.isfinite(r)] = 0.0
        out += add
    # MS2
    ms2x = df['entropy_similarity'].astype(float).map(_ms2_scalar).to_numpy()
    out += norm.logpdf(ms2x, params['ms2']['mu_tp'], params['ms2']['sd_tp']) \
         - norm.logpdf(ms2x, params['ms2']['mu_fp'], params['ms2']['sd_fp'])
    return pd.Series(out, index=df.index, name='logLR_hit')

cset['logLR_hit'] = per_hit_logLR_inline(cset, params)
cset['LR_hit']    = np.exp(cset['logLR_hit'])

# ---------------- 6) Aggregate → structure; local; global posterior ----------------
struct = aggregate_to_structure(cset, quality_cols=('entropy_similarity',), topK=None, lambda_var=0.0)
struct = local_rank_probability(struct)
post   = compute_global_posterior(struct, prior_policy='uniform', lr_nota=LR_NOTA, prior_nota=PRIOR_NOTA)

# ---------------- 7) Assertions + top calls ----------------
assertions = make_assertion_table(post, features_df=cset)
top_calls  = make_top_calls(assertions, nota_threshold=0.60)

# ---------------- 8) Per-hit probabilities (split structure mass) ----------------
hitsc = cset.merge(struct[['wiki_id','library_id','C_loc']], on=['wiki_id','library_id'], how='left')
hitsc = hitsc.merge(post[['wiki_id','library_id','post','P_NoTA']].drop_duplicates(),
                    on=['wiki_id','library_id'], how='left')

hitsc['_w'] = hitsc['entropy_similarity'].astype(float).clip(lower=1e-6).fillna(1e-6)
hitsc['_w'] = hitsc['_w'] / hitsc.groupby(['wiki_id','library_id'])['_w'].transform('sum')

# num = hitsc['_w'] * hitsc['LR_hit']
# den = hitsc.groupby(['wiki_id','library_id'])[num.name].transform('sum').replace(0, np.nan)
# hitsc['P_h_given_struct']  = num / den

hitsc['_num'] = hitsc['_w'] * hitsc['LR_hit']
hitsc['_den'] = hitsc.groupby(['wiki_id','library_id'])['_num'].transform('sum').replace(0, np.nan)
hitsc['P_h_given_struct'] = hitsc['_num'] / hitsc['_den']


import numpy as np

def _logsumexp(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))

def structure_level_posts(hitsc, structure_key="library_id", nota_col="P_NoTA"):
    """
    hitsc: per-hit table containing at least:
        ['wiki_id', structure_key, 'theta_hit'] and spectrum-level nota_col
      (If you only have LR_hit, set theta_hit = np.log(LR_hit).)
    Returns: per-(wiki_id, structure) posterior without hit-level dilution.
    """
    df = hitsc.copy()
    if "theta_hit" not in df.columns:
        if "LR_hit" not in df.columns:
            raise ValueError("Need either 'theta_hit' or 'LR_hit' in hitsc.")
        df["theta_hit"] = np.log(np.clip(df["LR_hit"].astype(float), 1e-300, None))

    # 1) Collapse duplicates: combine multiple hits of the SAME structure
    #    Option A (recommended): independent-evidence assumption → log-sum-exp over hit thetas
    theta_struct = (df
        .groupby(["wiki_id", structure_key])["theta_hit"]
        .apply(lambda s: _logsumexp(s.to_numpy()))
        .reset_index()
        .rename(columns={"theta_hit":"theta_struct"}))

    # 2) Local softmax across *structures only* (no hit-level entries in the denominator)
    def _softmax_group(g):
        x = g["theta_struct"].to_numpy()
        m = np.max(x)
        ex = np.exp(x - m)
        g["C_loc_struct"] = ex / ex.sum()
        return g

    theta_struct = (theta_struct
        .groupby("wiki_id", as_index=False, group_keys=False)
        .apply(_softmax_group))

    # 3) Bring in spectrum-level NoTA and compute structure posterior
    nota = (df[["wiki_id", nota_col]].drop_duplicates("wiki_id"))
    out = theta_struct.merge(nota, on="wiki_id", how="left")
    out["post_struct"] = (1.0 - out[nota_col].astype(float)) * out["C_loc_struct"]

    # Optional QC aggregations on the structure bin
    # e.g., max entropy, median delta_ppm, member count
    qc_aggs = {}
    for c, fn in [("entropy_similarity","max"), ("delta_ppm","median"),
                  ("delta_rt","median"), ("zrt","median")]:
        if c in df.columns:
            qc_aggs[c] = fn
    qc = (df.groupby(["wiki_id", structure_key])
            .agg({**qc_aggs, structure_key:"size"})
            .rename(columns={structure_key:"member_hits"})
            .reset_index())
    out = out.merge(qc, on=["wiki_id", structure_key], how="left")
    return out  # one row per (wiki_id, structure)

# Collapse to structure-level posteriors (no dilution by duplicate hits)
struct_posts = structure_level_posts(hitsc, structure_key="library_id", nota_col="P_NoTA")

# Use this for assertions instead of hit-level diluted posts
assertions = struct_posts.rename(columns={
    "theta_struct":"theta",
    "C_loc_struct":"C_loc",
    "post_struct":"post"
})
assertions["PEP"] = 1.0 - assertions["post"]

hitsc['match_post_local']  = hitsc['C_loc'] * hitsc['P_h_given_struct']
hitsc['match_post_global'] = hitsc['post']  * hitsc['P_h_given_struct']
hitsc['PEP_hit']           = 1.0 - hitsc['match_post_global']

hit_scores_out = hitsc[[
    'wiki_id','library_id','name','adduct','entropy_similarity','delta_ppm',
    'logLR_hit','LR_hit','P_h_given_struct','C_loc','post','P_NoTA',
    'match_post_local','match_post_global','PEP_hit'
]].sort_values(['wiki_id','match_post_global'], ascending=[True, False])

# ---------------- 9) Save + quick prints ----------------
assertions.to_csv(f"{OUTDIR}/assertions.csv", index=False)
top_calls.to_csv(f"{OUTDIR}/top_calls.csv", index=False)
hit_scores_out.to_csv(f"{OUTDIR}/hit_scores.csv", index=False)

print("Assertions,", assertions.shape)
print("Spectra total:", spectrum['wiki_id'].nunique(),
      "| with candidates:", assertions['wiki_id'].nunique())
print("top_calls shape:", top_calls.shape,
      "| abstain %:", float(top_calls['abstain'].mean() or 0))
print("Wrote:\n", f"{OUTDIR}/assertions.csv\n{OUTDIR}/top_calls.csv\n{OUTDIR}/hit_scores.csv")
