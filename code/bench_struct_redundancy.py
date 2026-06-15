"""bench_struct_redundancy.py — does within-bin MolRex structural similarity
add ΔAUC when entropy_similarity / sim_gap are removed?

The question: the original struct_sim_gap bench was null when added on top of
the full 18-feature GBM. Univariate analysis (2026-05-29) showed cos_to_gbm_top1
discriminates TP/FP at AUC 0.7175. So the signal exists — is it being absorbed
by entropy_similarity / sim_gap, or is the GBM just not finding it?

Ablation arms (5-fold IK14 GroupKFold, same setup as bench_harness):
  1. baseline_full              — production 18-feature GBM
  2. baseline_full + struct     — add max_cos_to_siblings (control; expect ≈0)
  3. minus_entropy              — drop entropy_similarity
  4. minus_entropy + struct     — does struct recover the lost signal?
  5. minus_entropy_simgap       — drop both entropy_similarity AND sim_gap
  6. minus_entropy_simgap + struct  — does struct recover both?

If arms 4 and 6 recover ΔAUC over their reduced baselines, the redundancy
hypothesis is confirmed. If they don't, the signal lives somewhere else.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FT_PATH = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
CAND_PATH = os.path.join(ROOT, 'data', 'candidate_scores_v2.csv')
EMB_PATH = os.path.join(ROOT, 'data', 'molrex_embeddings.npz')
OUT = os.path.join(ROOT, 'data', 'bench_struct_redundancy_summary.csv')

BASE_FULL = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
]
CAT = ['hit_adduct_cat', 'db', 'polarity']
XGB_P = {
    'objective': 'binary:logistic', 'eval_metric': 'auc',
    'tree_method': 'hist', 'max_depth': 5, 'learning_rate': 0.05,
    'subsample': 0.85, 'colsample_bytree': 0.85, 'min_child_weight': 5,
    'reg_alpha': 0.1, 'reg_lambda': 1.0, 'seed': 42, 'verbosity': 0,
}
N_ROUND = 500


def compute_sibling_stats():
    """For each (bin, candidate), compute MolRex cosine to siblings within bin."""
    print('Loading MolRex embeddings...')
    e = np.load(EMB_PATH, allow_pickle=True)
    iks = e['ik14']
    V = e['embedding']
    Vn = V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)
    ik2idx = {ik: i for i, ik in enumerate(iks)}

    print('Loading candidates...')
    cand = pd.read_csv(CAND_PATH, usecols=['wiki_id', 'hit_ik14', 'gbm_cal'])
    cand = cand[cand['hit_ik14'].notna() & (cand['hit_ik14'] != '')].copy()
    cand['mol_idx'] = cand['hit_ik14'].map(ik2idx)

    out_rows = []
    g = cand.groupby('wiki_id')
    for wid, grp in g:
        idx = grp['mol_idx'].values
        valid = ~pd.isna(idx)
        if valid.sum() < 2:
            for _, r in grp.iterrows():
                out_rows.append({
                    'wiki_id': wid, 'hit_ik14': r['hit_ik14'],
                    'max_cos_to_siblings': np.nan,
                    'mean_cos_to_siblings': np.nan,
                    'cos_to_gbm_top1': np.nan,
                })
            continue
        idx_int = np.array([int(x) if not pd.isna(x) else -1 for x in idx])
        V_g = Vn[idx_int[idx_int >= 0]]
        C = V_g @ V_g.T
        np.fill_diagonal(C, np.nan)
        gbm_vals = grp['gbm_cal'].values
        valid_pos = np.where(idx_int >= 0)[0]
        cur = 0
        for k_local, (_, r) in enumerate(grp.iterrows()):
            if idx_int[k_local] < 0:
                out_rows.append({
                    'wiki_id': wid, 'hit_ik14': r['hit_ik14'],
                    'max_cos_to_siblings': np.nan,
                    'mean_cos_to_siblings': np.nan,
                    'cos_to_gbm_top1': np.nan,
                })
                continue
            sib = C[cur]
            sib_valid = sib[~np.isnan(sib)]
            if sib_valid.size == 0:
                out_rows.append({
                    'wiki_id': wid, 'hit_ik14': r['hit_ik14'],
                    'max_cos_to_siblings': np.nan,
                    'mean_cos_to_siblings': np.nan,
                    'cos_to_gbm_top1': np.nan,
                })
                cur += 1
                continue
            gbm_excl = gbm_vals.copy().astype(float)
            gbm_excl[k_local] = -np.inf
            top1_full = int(np.argmax(gbm_excl))
            if idx_int[top1_full] < 0:
                cos_top1 = np.nan
            else:
                top1_local = int(np.where(valid_pos == top1_full)[0][0])
                cos_top1 = float(C[cur, top1_local])
            out_rows.append({
                'wiki_id': wid, 'hit_ik14': r['hit_ik14'],
                'max_cos_to_siblings': float(np.nanmax(sib)),
                'mean_cos_to_siblings': float(np.nanmean(sib)),
                'cos_to_gbm_top1': cos_top1,
            })
            cur += 1
    sib_df = pd.DataFrame(out_rows)
    print(f'  sibling stats computed for {len(sib_df):,} (bin, candidate) rows')
    print(f'  with usable mean_cos_to_siblings: {sib_df["mean_cos_to_siblings"].notna().sum():,}')
    return sib_df


def prep_features(df, feats_num, feats_cat):
    X = df.copy()
    for c in feats_num:
        if c not in X.columns:
            X[c] = np.nan
        X[c] = pd.to_numeric(X[c], errors='coerce')
    for c in feats_cat:
        s = X[c].astype(str).fillna('missing')
        cats = sorted(s.unique().tolist())
        idx = {v: i for i, v in enumerate(cats)}
        X[c] = s.map(idx).astype('int32')
    return X[feats_num + feats_cat]


def ece(p, y, nbins=10):
    edges = np.linspace(0, 1, nbins + 1)
    bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
    total, n = 0.0, 0
    for b in range(nbins):
        m = bi == b
        if m.sum() == 0: continue
        total += m.sum() * abs(p[m].mean() - y[m].mean())
        n += m.sum()
    return total / max(n, 1)


def run_arm(top1, feats_num, groups, tag):
    X = prep_features(top1, feats_num, CAT)
    y = top1['hit_label'].values
    n_num, n_cat = len(feats_num), len(CAT)
    oof = np.full(len(top1), np.nan)
    fold_aucs = []
    gkf = GroupKFold(n_splits=5)
    for fold, (tr, te) in enumerate(gkf.split(top1, y, groups)):
        dtr = xgb.DMatrix(X.iloc[tr], label=y[tr], enable_categorical=True,
                          feature_types=['q'] * n_num + ['c'] * n_cat)
        m = xgb.train(XGB_P, dtr, num_boost_round=N_ROUND)
        dte = xgb.DMatrix(X.iloc[te], enable_categorical=True,
                          feature_types=['q'] * n_num + ['c'] * n_cat)
        oof[te] = m.predict(dte)
        fold_aucs.append(roc_auc_score(y[te], oof[te]))
    auc_oof = roc_auc_score(y, oof)
    iso = IsotonicRegression(out_of_bounds='clip').fit(oof, y)
    oc = iso.transform(oof)
    print(f'  [{tag:42s}] AUC={auc_oof:.4f}  Brier_cal={brier_score_loss(y, oc):.4f}  '
          f'ECE_cal={ece(oc, y):.4f}  features={n_num+n_cat}')
    return {'tag': tag, 'auc_oof': auc_oof, 'fold_aucs': fold_aucs,
            'brier_cal': brier_score_loss(y, oc), 'ece_cal': ece(oc, y)}


def main():
    sib = compute_sibling_stats()
    print(f'\nLoading {FT_PATH}')
    ft = pd.read_csv(FT_PATH, low_memory=False)
    print(f'  {len(ft):,} rows, {ft["wiki_id"].nunique():,} bins')
    labeled = ft[ft['spectrum_label'].isin(['TP', 'FP'])].copy()
    top1_idx = labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1 = ft.loc[top1_idx].reset_index(drop=True)
    print(f'  top-1 per bin in labeled set: {len(top1):,}')

    top1 = top1.merge(sib, on=['wiki_id', 'hit_ik14'], how='left')
    print(f'  merged sibling stats; valid mean_cos: {top1["mean_cos_to_siblings"].notna().sum():,}')

    groups = top1['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    minus_e = [c for c in BASE_FULL if c != 'entropy_similarity']
    minus_e_s = [c for c in BASE_FULL if c not in ('entropy_similarity', 'sim_gap')]

    arms = []
    print('\n=== Arms ===')
    arms.append(run_arm(top1, BASE_FULL, groups,
                        'A1 baseline_full'))
    arms.append(run_arm(top1, BASE_FULL + ['max_cos_to_siblings', 'cos_to_gbm_top1'], groups,
                        'A2 baseline_full + struct'))
    arms.append(run_arm(top1, minus_e, groups,
                        'A3 minus_entropy_sim'))
    arms.append(run_arm(top1, minus_e + ['max_cos_to_siblings', 'cos_to_gbm_top1'], groups,
                        'A4 minus_entropy_sim + struct'))
    arms.append(run_arm(top1, minus_e_s, groups,
                        'A5 minus_entropy_sim_AND_sim_gap'))
    arms.append(run_arm(top1, minus_e_s + ['max_cos_to_siblings', 'cos_to_gbm_top1'], groups,
                        'A6 minus_entropy_sim_AND_sim_gap + struct'))

    base = arms[0]['auc_oof']
    for a in arms:
        a['delta_vs_A1'] = a['auc_oof'] - base
    df = pd.DataFrame(arms)
    df.to_csv(OUT, index=False)

    print(f'\n=== Summary (Δ vs A1 baseline) ===')
    for a in arms:
        print(f'  {a["tag"]:50s}  AUC={a["auc_oof"]:.4f}  Δ={a["delta_vs_A1"]:+.4f}')

    print(f'\nWrote {OUT}')

    print('\n=== Critical comparisons ===')
    by = {a['tag']: a['auc_oof'] for a in arms}
    print(f'  Does struct add when entropy IS in model?  A2 - A1 = {by[arms[1]["tag"]] - by[arms[0]["tag"]]:+.4f}')
    print(f'  Does struct recover from minus_entropy?    A4 - A3 = {by[arms[3]["tag"]] - by[arms[2]["tag"]]:+.4f}')
    print(f'  Does struct recover from minus_both?       A6 - A5 = {by[arms[5]["tag"]] - by[arms[4]["tag"]]:+.4f}')


if __name__ == '__main__':
    main()
