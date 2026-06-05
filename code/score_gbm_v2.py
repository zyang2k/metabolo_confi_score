"""
score_gbm_v2.py — Full switch to GBM (XGBoost) on the v2 feature table.

Trains on labeled top-1-by-entropy-similarity rows (TP + FP only), excluding
blanks. 5-fold GroupKFold by `anno_ik14` for OOF scoring, then isotonic
calibration applied on top so the output is a defensible 0–100% probability.

Features used:
  Per-hit numeric:  entropy_similarity, sim_gap, signed_delta_rt, delta_mda,
                    forward_cosine, reverse_cosine, cov_count, cov_int
  Per-spectrum:     spectral_entropy, n_candidates
  Reinforcement:    n_candidate_adducts, compound_has_ok_adduct
  Flags:            hit_is_isf, hit_is_dubious, hit_isf_no_ok
  Categorical:      hit_adduct_cat, db, polarity

Excluded:
  hit_theoretical_mz, precursor_mz — polarity confound (pos m/z systematically
  larger); would let the model cheat via polarity.
  rank — consequence of entropy_similarity ranking, not independent.
  Identifiers (wiki_id, hit_ik14, anno_ik14, etc.).

Compared to calibrated Bayesian on:
  • OOF AUC / Brier / ECE
  • Deck case confidences (PEP, trig [Cat]+, deoxycarnitine [Cat]+, trig baseline)
  • Min ground-truth match rate
  • NoTA behavior on blanks

Outputs:
  data/deliverable_scores_gbm.csv
  data/candidate_scores_gbm.csv
  figures/gbm_vs_bayesian_20260423.png
"""

import os
import sys
import json
import shutil
import argparse
import datetime
import warnings
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # so `python code/score_gbm_v2.py` resolves siblings
from relational_graph import Graph, REL_COLS  # gate-2 frozen-graph confusability features

warnings.filterwarnings('ignore')
matplotlib_backend = None
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib_backend = 'Agg'
except ImportError:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEATURE_TABLE       = os.path.join(ROOT, 'data', 'feature_table_v2.csv')
QUERY_PEAKS         = os.path.join(ROOT, 'data', 'query_peaks_cache_v2.json')
# Primary GBM output — this is the production scorer as of 2026-04-23.
# Bayesian output is kept as a companion at deliverable_scores_bayesian_v2.csv.
DELIVERABLE_GBM_OUT = os.path.join(ROOT, 'data', 'deliverable_scores_v2.csv')
CANDIDATES_GBM_OUT  = os.path.join(ROOT, 'data', 'candidate_scores_v2.csv')
# Keep the _gbm.csv legacy paths as aliases so older notebooks still resolve.
DELIVERABLE_GBM_ALIAS = os.path.join(ROOT, 'data', 'deliverable_scores_gbm.csv')
CANDIDATES_GBM_ALIAS  = os.path.join(ROOT, 'data', 'candidate_scores_gbm.csv')
FIG_OUT             = os.path.join(ROOT, 'figures', 'gbm_vs_bayesian_20260423.png')
BAYESIAN_COMPARISON = os.path.join(ROOT, 'data', 'deliverable_scores_bayesian_v2.csv')
MIN_VERIFIED        = os.path.join(ROOT, 'benchmark', '_internal', 'min_verified_entries.csv')

NUMERIC_FEATURES = [
    'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
    'forward_cosine', 'reverse_cosine', 'cov_count', 'cov_int',
    'spectral_entropy', 'n_candidates', 'n_candidate_adducts',
    'compound_has_ok_adduct', 'hit_is_isf', 'hit_is_dubious', 'hit_isf_no_ok',
    # gate-2 relational confusability (frozen-graph; see relational_graph.py)
    'rel_nn_sim', 'rel_n_confirmed_nbr',
]
# per-bin frame columns the relational graph needs (one row per wiki_id)
REL_BIN_COLS = ['wiki_id', 'precursor_mz', 'polarity', 'anno_ik14', 'anno_name_lower', 'hit_label']
CATEGORICAL_FEATURES = ['hit_adduct_cat', 'db', 'polarity']
ALL_FEATURES = NUMERIC_FEATURES + CATEGORICAL_FEATURES

# Bootstrap-ensemble config — used to compute ensemble_sd column for curator triage.
# K bootstrap copies of the GBM, group-blocked by anno_ik14. Production confidence stays
# single-model; ensemble_sd is shipped alongside as a per-row epistemic uncertainty.
K_BOOTSTRAP = 10
BOOTSTRAP_SEED = 42


def block_bootstrap_indices(group_keys, rng):
    """Resample IK14 groups with replacement; empty/NaN treated as unique groups."""
    keys = np.asarray(group_keys, dtype=object).copy()
    for i in range(len(keys)):
        if keys[i] is None or (isinstance(keys[i], float) and np.isnan(keys[i])) or keys[i] == '':
            keys[i] = f'__no_ik14_{i}'
    unique_keys, inv = np.unique(keys, return_inverse=True)
    sampled_groups = rng.choice(len(unique_keys), size=len(unique_keys), replace=True)
    parts = [np.where(inv == g)[0] for g in sampled_groups]
    return np.concatenate(parts)


def prep_features(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce numeric columns + label-encode categoricals to integer codes."""
    X = df.copy()
    for c in NUMERIC_FEATURES:
        if c not in X.columns:
            X[c] = np.nan
        X[c] = pd.to_numeric(X[c], errors='coerce')
    # Categorical: label-encode with stable mapping
    cat_maps = {}
    for c in CATEGORICAL_FEATURES:
        if c not in X.columns:
            X[c] = 'missing'
        s = X[c].astype(str).fillna('missing')
        categories = sorted(s.unique().tolist())
        mapping = {v: i for i, v in enumerate(categories)}
        X[c] = s.map(mapping).astype('int32')
        cat_maps[c] = mapping
    return X[ALL_FEATURES], cat_maps


def train_xgb(X_train, y_train, cat_feature_names):
    """Train XGBoost with reasonable defaults on training subset."""
    dtrain = xgb.DMatrix(X_train, label=y_train,
                         enable_categorical=True,
                         feature_types=['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES))
    params = {
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'tree_method': 'hist',
        'max_depth': 5,
        'learning_rate': 0.05,
        'subsample': 0.85,
        'colsample_bytree': 0.85,
        'min_child_weight': 5,
        'reg_alpha': 0.1,
        'reg_lambda': 1.0,
        'seed': 42,
        'verbosity': 0,
    }
    num_rounds = 500
    return xgb.train(params, dtrain, num_boost_round=num_rounds)


def _isotonic_to_dict(iso: IsotonicRegression) -> dict:
    """Serialize a fitted IsotonicRegression to a portable JSON-able dict.

    Only the breakpoint arrays are needed: predict() linearly interpolates over
    (X_thresholds_, y_thresholds_) and, with out_of_bounds='clip', clamps to the
    end values — which is exactly numpy.interp's default behaviour. The inference
    module reconstructs the calibrator with np.interp(x, x_thresholds, y_thresholds),
    avoiding any sklearn/pickle version coupling. selftest.py asserts bit-parity.
    """
    return {
        'x_thresholds': np.asarray(iso.X_thresholds_, dtype=float).tolist(),
        'y_thresholds': np.asarray(iso.y_thresholds_, dtype=float).tolist(),
        'increasing': bool(iso.increasing_),
        'out_of_bounds': 'clip',
    }


def freeze_artifacts(out_dir, final_model, final_bootstrap, iso, iso_ens,
                     cat_maps, metrics):
    """Serialize the in-memory production artifacts to `out_dir`.

    Writes exactly what this run trained — the final XGBoost booster, the K
    bootstrap boosters, both isotonic calibrators, the categorical encodings,
    the feature contract, and provenance metadata. Because these are the same
    objects used to produce deliverable_scores_v2.csv, the frozen package cannot
    drift from production. Boosters use XGBoost's native JSON (portable across
    2.x). Everything else is plain JSON.
    """
    os.makedirs(out_dir, exist_ok=True)

    final_model.save_model(os.path.join(out_dir, 'gbm_final.json'))
    for k, mk in enumerate(final_bootstrap):
        mk.save_model(os.path.join(out_dir, f'gbm_bootstrap_{k:02d}.json'))

    with open(os.path.join(out_dir, 'isotonic_main.json'), 'w') as f:
        json.dump(_isotonic_to_dict(iso), f, indent=2)
    with open(os.path.join(out_dir, 'isotonic_ensemble.json'), 'w') as f:
        json.dump(_isotonic_to_dict(iso_ens), f, indent=2)

    meta = {
        'model': 'metabolo_confi_score GBM v2 (XGBoost + isotonic)',
        'frozen_at': datetime.datetime.now().isoformat(timespec='seconds'),
        'source_script': 'code/score_gbm_v2.py',
        'numeric_features': NUMERIC_FEATURES,
        'categorical_features': CATEGORICAL_FEATURES,
        'all_features': ALL_FEATURES,
        'feature_types': ['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES),
        'cat_maps': cat_maps,
        'k_bootstrap': K_BOOTSTRAP,
        'bootstrap_seed': BOOTSTRAP_SEED,
        'final_model_file': 'gbm_final.json',
        'bootstrap_model_files': [f'gbm_bootstrap_{k:02d}.json' for k in range(len(final_bootstrap))],
        'isotonic_main_file': 'isotonic_main.json',
        'isotonic_ensemble_file': 'isotonic_ensemble.json',
        'adduct_taxonomy_file': 'adduct_taxonomy_oliver.csv',
        'relational_graph_file': 'relational_graph.json',
        'versions': {
            'python': sys.version.split()[0],
            'xgboost': xgb.__version__,
            'sklearn': sklearn.__version__,
            'numpy': np.__version__,
            'pandas': pd.__version__,
        },
        'metrics': metrics,
    }
    with open(os.path.join(out_dir, 'model_meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    # The feature engineering needs the adduct taxonomy; bundle the exact copy used.
    shutil.copy(os.path.join(ROOT, 'data', 'adduct_taxonomy_oliver.csv'),
                os.path.join(out_dir, 'adduct_taxonomy_oliver.csv'))

    print(f'\nFroze model artifacts → {out_dir}')
    print(f'  final booster + {len(final_bootstrap)} bootstrap boosters, 2 isotonic calibrators, meta + taxonomy')


def main(freeze_dir=None):
    print(f'Reading {FEATURE_TABLE}')
    ft = pd.read_csv(FEATURE_TABLE, low_memory=False)
    print(f'  {len(ft):,} rows x {len(ft.columns)} cols, {ft["wiki_id"].nunique():,} spectra')

    # Exclude candidates without a resolved library structure (empty hit_ik14). These mostly
    # come from NIST23 entries with no SMILES; GNPS/MassBank empty-IK14 hits were already
    # dropped upstream in build_features_v2.py. A scored confidence on an unresolved structure
    # is meaningless to a curator, so we drop them everywhere — training and scoring alike.
    # Removes name-fallback labeling as a side effect (was 42 training rows, 15 hit_label=1).
    before_smiles = len(ft)
    ft = ft[ft['hit_ik14'].fillna('').ne('')].reset_index(drop=True)
    print(f'  Dropped {before_smiles - len(ft):,} rows with empty hit_ik14; '
          f'kept {len(ft):,} candidate rows with resolved SMILES')

    # Apply Oliver-reviewed label overrides (data/oliver_extra_labels_*.csv). These promote
    # explicitly-reviewed blanks into training, flip mis-curated TPs to FPs where Oliver
    # disagreed, and confirm existing labels for trustworthy-iso fitting. The override sets
    # spectrum_label, hit_label, and anno_ik14 (to hit_ik14 of the top-1 candidate so the
    # downstream anno_ik14 filter accepts it). Tracked as `oliver_override_source`.
    overrides_path = os.path.join(ROOT, 'data', 'oliver_extra_labels_20260520.csv')
    if os.path.exists(overrides_path):
        overrides = pd.read_csv(overrides_path)
        ft['oliver_override_source'] = ''
        n_applied = 0
        for _, ov in overrides.iterrows():
            mask = ft['wiki_id'] == ov['wiki_id']
            if not mask.any():
                continue
            sub = ft[mask]
            top1_pos = sub['entropy_similarity'].idxmax()
            ft.loc[mask, 'spectrum_label'] = ov['spectrum_label']
            ft.loc[top1_pos, 'hit_label'] = int(ov['hit_label'])
            ft.loc[top1_pos, 'anno_ik14'] = ft.loc[top1_pos, 'hit_ik14']
            ft.loc[mask, 'oliver_override_source'] = ov['source']
            n_applied += 1
        print(f'  Applied {n_applied}/{len(overrides)} Oliver overrides:')
        applied_top1 = ft[(ft['oliver_override_source'] != '') &
                          (ft.groupby('wiki_id')['entropy_similarity'].transform('idxmax') == ft.index)]
        applied_summary = applied_top1.groupby(['oliver_override_source','spectrum_label']).size()
        for (src, lbl), n in applied_summary.items():
            print(f'    {src:25s} → {lbl}: {n}')
    else:
        ft['oliver_override_source'] = ''
        print(f'  No overrides file found ({overrides_path}); proceeding without.')

    # Training = top-1 rows of labeled spectra (TP + FP only); blanks scored but not trained.
    # Also require anno_ik14 (curator structure) for labeled rows — should be ~universal but
    # guard the 1 anomalous TP row without anno_ik14.
    labeled_mask = ft['spectrum_label'].isin(['TP', 'FP']) & ft['anno_ik14'].fillna('').ne('')
    ft_labeled = ft[labeled_mask]
    top1_idx = ft_labeled.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1_train = ft.loc[top1_idx].reset_index(drop=True)
    labels = top1_train['hit_label'].values
    prior_tp = labels.mean()
    print(f'\nTraining set: {len(top1_train):,} top-1 rows (TP+FP)   prior={prior_tp:.3f}')

    # --- Gate-2 relational confusability features (frozen training graph) ---
    # The graph is built ONLY from labeled training bins, so scoring-time look-ups can
    # never read a label that wasn't in training (project_relational_kg gate 2). OOF folds
    # below rebuild the graph per-fold (train-only) so the reported AUC stays honest.
    print('\nComputing relational confusability features (frozen training graph)...')
    rel_peaks = {k: np.asarray(v, dtype=np.float64)
                 for k, v in json.load(open(QUERY_PEAKS)).items() if v}
    train_bins = top1_train[REL_BIN_COLS].copy()
    full_graph = Graph(train_bins, rel_peaks)                 # frozen graph = all training bins
    rel_train_full = full_graph.features(train_bins).set_index('wiki_id')   # final-fit features
    ft_bins = ft[REL_BIN_COLS].drop_duplicates('wiki_id')
    rel_all = full_graph.features(ft_bins).set_index('wiki_id')             # scoring features
    for c in REL_COLS:
        top1_train[c] = top1_train['wiki_id'].map(rel_train_full[c]).fillna(0.0).values
        ft[c] = ft['wiki_id'].map(rel_all[c]).fillna(0.0).values
    n_graph_nodes = sum(len(g['wid']) for g in full_graph.by_pol.values())
    print(f'  graph nodes: {n_graph_nodes:,}   bins scored: {len(rel_all):,}   '
          f'rel_nn_sim mean(train)={top1_train["rel_nn_sim"].mean():.3f}  '
          f'n_confirmed mean={top1_train["rel_n_confirmed_nbr"].mean():.3f}')

    # is_oliver_reviewed: FPs (yy_) are all explicitly flagged; among TPs, only those whose
    # anno_name appears in ≥2 spectra were reviewed for RT conflict. Singleton-name TPs were
    # auto-annotated by BinBase and not verified. Rows touched by the override CSV are also
    # forced to is_oliver_reviewed=True since Oliver explicitly looked at them. Used to fit
    # isotonic on a clean label slice.
    tp_mask = top1_train['spectrum_label'] == 'TP'
    fp_mask = top1_train['spectrum_label'] == 'FP'
    name_counts = top1_train[tp_mask].groupby('anno_name_lower')['wiki_id'].nunique()
    dup_names = set(name_counts[name_counts >= 2].index) - {''}
    is_dup_name = top1_train['anno_name_lower'].fillna('').isin(dup_names)
    has_override = top1_train.get('oliver_override_source', pd.Series('', index=top1_train.index)).fillna('').ne('')
    top1_train['is_oliver_reviewed'] = fp_mask | (tp_mask & is_dup_name) | has_override
    n_rev = int(top1_train['is_oliver_reviewed'].sum())
    n_override_trusted = int(has_override.sum())
    print(f'  Oliver-reviewed (trustworthy): {n_rev:,}  '
          f'({100 * n_rev / len(top1_train):.1f}%)  '
          f'[+{n_override_trusted} from override CSV]  — used for isotonic fit')

    # Prepare features for training and scoring
    X_train_df, cat_maps = prep_features(top1_train)
    X_all_df, _ = prep_features(ft)
    # Reuse train's cat mapping for consistency
    for c in CATEGORICAL_FEATURES:
        s = ft[c].astype(str).fillna('missing')
        mapping = cat_maps[c]
        max_code = max(mapping.values()) + 1 if mapping else 0
        X_all_df[c] = s.map(lambda v: mapping.get(v, max_code)).astype('int32')

    # --- 5-fold GroupKFold OOF on training set ---
    print('\nFitting OOF folds...')
    groups = top1_train['anno_ik14'].fillna('').values.copy()
    for i in range(len(groups)):
        if groups[i] == '':
            groups[i] = f'__no_ik14_{i}'

    gkf = GroupKFold(n_splits=5)
    oof = np.full(len(top1_train), np.nan)
    oof_per_k = np.full((K_BOOTSTRAP, len(top1_train)), np.nan)
    fold_aucs = []
    for fold, (tr, te) in enumerate(gkf.split(top1_train, labels, groups)):
        X_tr = X_train_df.iloc[tr].copy()
        X_te = X_train_df.iloc[te].copy()
        y_tr = labels[tr]
        # Honest OOF: rebuild the relational graph from THIS fold's train bins only, then
        # recompute rel features for both train and held-out rows against it. (The full-graph
        # values in X_train_df would let a test bin's different-compound neighbours from its
        # own fold leak in — train-only graph removes that optimism.)
        fold_graph = Graph(top1_train.iloc[tr][REL_BIN_COLS], rel_peaks)
        rel_tr = fold_graph.features(top1_train.iloc[tr][REL_BIN_COLS]).set_index('wiki_id')
        rel_te = fold_graph.features(top1_train.iloc[te][REL_BIN_COLS]).set_index('wiki_id')
        for c in REL_COLS:
            X_tr[c] = pd.Series(top1_train.iloc[tr]['wiki_id'].values).map(rel_tr[c]).fillna(0.0).values
            X_te[c] = pd.Series(top1_train.iloc[te]['wiki_id'].values).map(rel_te[c]).fillna(0.0).values
        model = train_xgb(X_tr, y_tr, CATEGORICAL_FEATURES)
        dte = xgb.DMatrix(X_te,
                          enable_categorical=True,
                          feature_types=['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES))
        oof[te] = model.predict(dte)
        fold_auc = roc_auc_score(labels[te], oof[te])
        fold_aucs.append(fold_auc)

        # K=10 IK14-group-block bootstrap models — used to estimate ensemble_sd
        tr_group_keys = top1_train.iloc[tr]['anno_ik14'].fillna('').values
        rng = np.random.default_rng(BOOTSTRAP_SEED + fold)
        for k in range(K_BOOTSTRAP):
            bs_idx = block_bootstrap_indices(tr_group_keys, rng)
            mk = train_xgb(X_tr.iloc[bs_idx], y_tr[bs_idx], CATEGORICAL_FEATURES)
            oof_per_k[k, te] = mk.predict(dte)

        print(f'  Fold {fold}: AUC={fold_auc:.4f}  (n_test={len(te):,})  '
              f'+ {K_BOOTSTRAP} bootstrap fits')

    valid = ~np.isnan(oof)
    auc_oof = roc_auc_score(labels[valid], oof[valid])
    brier_raw = brier_score_loss(labels[valid], oof[valid])
    print(f'\nGBM OOF AUC: {auc_oof:.4f}   folds: [{min(fold_aucs):.4f}, {max(fold_aucs):.4f}]')
    print(f'GBM OOF Brier (raw): {brier_raw:.4f}')

    # --- Isotonic calibration on TRUSTWORTHY OOF only ---
    # Fitting on the full OOF set calibrates to a population where 50% of positives are
    # unreviewed singleton-TPs. Fitting on is_oliver_reviewed=True rows targets the
    # curator-verified slice (yy_ FPs + duplicate-name TPs) — the labels that are actually
    # ground truth. Ships better calibration where it matters at the cost of slightly worse
    # Brier on the full set. See code/bench_bootstrap_ensemble.ipynb for the comparison.
    trust_oof = top1_train['is_oliver_reviewed'].values & valid
    oof_ens = oof_per_k.mean(axis=0)

    iso = IsotonicRegression(out_of_bounds='clip')
    iso.fit(oof[trust_oof], labels[trust_oof])
    iso_ens = IsotonicRegression(out_of_bounds='clip')
    iso_ens.fit(oof_ens[trust_oof], labels[trust_oof])

    oof_cal = iso.transform(oof[valid])
    brier_cal = brier_score_loss(labels[valid], oof_cal)
    print(f'GBM OOF Brier (trustworthy-isotonic, full eval):       {brier_cal:.4f}')
    print(f'GBM OOF Brier (trustworthy-isotonic, trustworthy eval): '
          f'{brier_score_loss(labels[trust_oof], iso.transform(oof[trust_oof])):.4f}')

    # ECE computation (same 10-bin scheme as Bayesian pipeline)
    nbins = 10
    edges = np.linspace(0, 1, nbins + 1)
    def ece(p, y):
        bi = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
        total, n = 0.0, 0
        for b in range(nbins):
            m = bi == b
            if m.sum() == 0: continue
            total += m.sum() * abs(p[m].mean() - y[m].mean())
            n += m.sum()
        return total / n
    ece_raw = ece(oof[valid], labels[valid])
    ece_cal = ece(oof_cal, labels[valid])
    print(f'GBM OOF ECE (full eval): raw={ece_raw:.4f}  calibrated={ece_cal:.4f}')
    # Calibration evaluated on the labels we actually trust:
    ece_trust = ece(iso.transform(oof[trust_oof]), labels[trust_oof])
    auc_trust = roc_auc_score(labels[trust_oof], oof[trust_oof])
    print(f'GBM OOF on trustworthy subset (n={trust_oof.sum():,}): '
          f'AUC={auc_trust:.4f}  ECE_cal={ece_trust:.4f}')

    # --- Train final model on ALL labeled top-1 rows ---
    print('\nTraining final model on full labeled top-1 set...')
    final_model = train_xgb(X_train_df, labels, CATEGORICAL_FEATURES)

    # --- K=10 bootstrap final models for ensemble_sd at inference ---
    print(f'Training {K_BOOTSTRAP} bootstrap final models on full labeled set...')
    all_group_keys = top1_train['anno_ik14'].fillna('').values
    rng = np.random.default_rng(BOOTSTRAP_SEED + 1000)
    final_bootstrap = []
    for k in range(K_BOOTSTRAP):
        bs_idx = block_bootstrap_indices(all_group_keys, rng)
        mk = train_xgb(X_train_df.iloc[bs_idx], labels[bs_idx], CATEGORICAL_FEATURES)
        final_bootstrap.append(mk)
    print(f'  ...done.')

    # --- Freeze artifacts for deployment (opt-in) ---
    # Serializes the exact in-memory objects just trained, so the deployable
    # package is provably the production model. See deploy/masswiki_gbm/.
    if freeze_dir:
        metrics = {
            'n_train_top1': int(len(top1_train)),
            'prior_tp': float(prior_tp),
            'n_trustworthy': int(trust_oof.sum()),
            'oof_auc': float(auc_oof),
            'oof_brier_raw': float(brier_raw),
            'oof_brier_cal': float(brier_cal),
            'oof_ece_raw': float(ece_raw),
            'oof_ece_cal': float(ece_cal),
            'trustworthy_auc': float(auc_trust),
            'trustworthy_ece_cal': float(ece_trust),
        }
        freeze_artifacts(freeze_dir, final_model, final_bootstrap, iso, iso_ens,
                         cat_maps, metrics)
        n_nodes = full_graph.freeze(os.path.join(freeze_dir, 'relational_graph.json'))
        print(f'  froze relational graph: {n_nodes:,} nodes '
              f'(deploy inference module must compute REL_COLS against it via relational_graph.load_graph)')

    # Feature importances
    imp = final_model.get_score(importance_type='gain')
    imp_sorted = sorted(imp.items(), key=lambda x: -x[1])
    print('\nFeature importance (gain):')
    for feat, g in imp_sorted[:15]:
        print(f'  {feat:28s}  {g:>12.2f}')

    # --- Score every candidate in the full feature table ---
    print('\nScoring every candidate...')
    dall = xgb.DMatrix(X_all_df,
                       enable_categorical=True,
                       feature_types=['q'] * len(NUMERIC_FEATURES) + ['c'] * len(CATEGORICAL_FEATURES))
    raw_all = final_model.predict(dall)
    cal_all = iso.transform(raw_all)

    # Bootstrap ensemble: K predictions per candidate, push each through ensemble isotonic,
    # then sd across K = epistemic uncertainty for the curator triage column.
    bs_preds_cal = np.zeros((K_BOOTSTRAP, len(X_all_df)))
    for k, mk in enumerate(final_bootstrap):
        bs_preds_cal[k] = iso_ens.transform(mk.predict(dall))
    ensemble_sd_cal = bs_preds_cal.std(axis=0)
    ensemble_mean_cal = bs_preds_cal.mean(axis=0)

    ft_out = ft.copy()
    ft_out['gbm_raw'] = raw_all
    ft_out['gbm_cal'] = cal_all
    ft_out['ensemble_sd'] = ensemble_sd_cal
    ft_out['ensemble_mean_cal'] = ensemble_mean_cal

    # --- Per-spectrum deliverable using same pick_scored_rows semantics ---
    # Reuse the branch logic: curator_ik14 > curator_name_fallback > top1_by_sim
    ft_out['hit_ik14'] = ft_out['hit_ik14'].fillna('')
    ft_out['anno_ik14'] = ft_out['anno_ik14'].fillna('')
    has_curator = (ft_out['spectrum_label'] == 'TP') & (ft_out['anno_ik14'] != '')

    ik14_match = ft_out[has_curator & (ft_out['hit_ik14'] == ft_out['anno_ik14'])]
    ik14_picks_idx = ik14_match.groupby('wiki_id')['entropy_similarity'].idxmax()
    ik14_picks = ft_out.loc[ik14_picks_idx].copy()
    ik14_picks['pick_source'] = 'curator_ik14'

    covered = set(ik14_picks['wiki_id'])
    hit_name_lower = ft_out['name'].fillna('').str.strip().str.lower()
    anno_name_lower = ft_out['anno_name_lower'].fillna('').str.strip().str.lower()
    name_match = ft_out[
        has_curator & (~ft_out['wiki_id'].isin(covered)) &
        (ft_out['hit_ik14'] == '') & (anno_name_lower != '') &
        (hit_name_lower == anno_name_lower)
    ]
    name_picks_idx = name_match.groupby('wiki_id')['entropy_similarity'].idxmax()
    name_picks = ft_out.loc[name_picks_idx].copy()
    name_picks['pick_source'] = 'curator_name_fallback'

    covered.update(name_picks['wiki_id'])
    remaining = ft_out[~ft_out['wiki_id'].isin(covered)]
    top1_pick_idx = remaining.groupby('wiki_id')['entropy_similarity'].idxmax()
    top1_picks = ft_out.loc[top1_pick_idx].copy()
    top1_picks['pick_source'] = 'top1_by_sim'

    picks = pd.concat([ik14_picks, name_picks, top1_picks], ignore_index=True)
    assert picks['wiki_id'].nunique() == len(picks)

    picks['confidence'] = picks['gbm_cal']
    picks['confidence_pct'] = (picks['confidence'] * 100).round(1)
    picks['confidence_raw'] = picks['gbm_raw']
    # ensemble_sd already carried through from ft_out; expose a rounded % version for readers.
    picks['ensemble_sd_pct'] = (picks['ensemble_sd'] * 100).round(1)

    out_cols = [
        'wiki_id', 'spectrum_label', 'hit_label', 'pick_source',
        'name', 'adduct', 'hit_ik14', 'anno_ik14',
        'confidence', 'confidence_pct', 'confidence_raw',
        'ensemble_sd', 'ensemble_sd_pct',
        'entropy_similarity', 'sim_gap', 'signed_delta_rt', 'delta_mda',
        'spectral_entropy', 'hit_adduct_cat', 'db',
    ]
    deliverable = picks[[c for c in out_cols if c in picks.columns]].rename(columns={'name': 'annotation'})
    deliverable.to_csv(DELIVERABLE_GBM_OUT, index=False)
    deliverable.to_csv(DELIVERABLE_GBM_ALIAS, index=False)  # legacy alias
    print(f'\nWrote {DELIVERABLE_GBM_OUT}: {len(deliverable):,} rows')
    print(f'  (also mirrored to {DELIVERABLE_GBM_ALIAS})')

    # Candidate-level dump
    cand_cols = ['wiki_id', 'hit_ik14', 'name', 'adduct', 'library_wiki_id',
                 'entropy_similarity', 'sim_gap', 'signed_delta_rt',
                 'gbm_raw', 'gbm_cal', 'ensemble_sd', 'ensemble_mean_cal',
                 'hit_label']
    candidate_out = ft_out[[c for c in cand_cols if c in ft_out.columns]]
    candidate_out.to_csv(CANDIDATES_GBM_OUT, index=False)
    candidate_out.to_csv(CANDIDATES_GBM_ALIAS, index=False)  # legacy alias
    print(f'Wrote {CANDIDATES_GBM_OUT}: {len(ft_out):,} rows')

    # --- Compare to Bayesian deliverable ---
    print('\n=== Comparison with Bayesian ===')
    try:
        bay = pd.read_csv(BAYESIAN_COMPARISON)[['wiki_id','confidence','confidence_raw','spectrum_label']].rename(
            columns={'confidence':'bay_cal','confidence_raw':'bay_raw','spectrum_label':'spec_label'})
        merged = deliverable[['wiki_id','confidence','confidence_raw','spectrum_label']].merge(bay, on='wiki_id', how='inner')
        for grp in ['TP', 'FP', 'blank']:
            sub = merged[merged['spectrum_label']==grp]
            print(f'\n  {grp} (n={len(sub):,})')
            print(f'    Bayesian median: {sub["bay_cal"].median()*100:5.1f}%   GBM median: {sub["confidence"].median()*100:5.1f}%')
            for thr in [0.5, 0.7, 0.9]:
                b_n = (sub['bay_cal']>=thr).sum(); g_n = (sub['confidence']>=thr).sum()
                print(f'    ≥ {thr:.1f}: Bayesian {100*b_n/max(len(sub),1):5.1f}%   GBM {100*g_n/max(len(sub),1):5.1f}%')
    except FileNotFoundError:
        print(f'  {BAYESIAN_COMPARISON} not found — skipping comparison')

    # --- Deck cases ---
    print('\n=== Deck cases (GBM confidence) ===')
    for wid, label in [('aPUDE1U/2244', 'PEP slide 4'),
                       ('aEKJ9AS/10123', 'trig [Cat]+ slide 7'),
                       ('aEKJ9AS/994', 'deoxycarnitine [Cat]+ slide 8'),
                       ('aEKJ9AS/1674', 'trig baseline')]:
        r = deliverable[deliverable['wiki_id'] == wid]
        if len(r):
            r = r.iloc[0]
            print(f'  {wid:22s} ({label:28s}): '
                  f'annotation="{str(r["annotation"])[:28]:28s}"  '
                  f'conf={r["confidence_pct"]:5.1f}%  raw={r["confidence_raw"]*100:5.1f}%  '
                  f'source={r["pick_source"]}')

    # --- Min validation ---
    try:
        min_df = pd.read_csv(MIN_VERIFIED)[['wiki_id','host_ik14','name']]
        mer = min_df.merge(deliverable[['wiki_id','hit_ik14','confidence','confidence_pct']], on='wiki_id', how='inner')
        mer['matches'] = (mer['hit_ik14'] == mer['host_ik14']).astype(int)
        print(f'\n=== Min validation (GBM) ===')
        print(f'  Pipeline matches Min: {mer["matches"].sum()}/{len(mer)} ({100*mer["matches"].mean():.1f}%)')
        print(f'  Median conf on Min-correct: {mer[mer["matches"]==1]["confidence"].median()*100:.1f}%')
        print(f'  ≥ 90% among Min-correct: {(mer[mer["matches"]==1]["confidence"]>=0.9).sum()}/{mer["matches"].sum()}')
    except FileNotFoundError:
        pass

    # --- Reliability plot ---
    if matplotlib_backend:
        print('\nPlotting GBM vs Bayesian reliability...')
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        for ax, (name, p) in zip(axes, [('GBM raw', oof[valid]), ('GBM calibrated', oof_cal)]):
            bin_idx = np.clip(np.digitize(p, edges[1:-1]), 0, nbins - 1)
            bc, bf, bn = [], [], []
            for b in range(nbins):
                m = bin_idx == b
                if m.sum() == 0: continue
                bc.append(p[m].mean()); bf.append(labels[valid][m].mean()); bn.append(m.sum())
            bc, bf, bn = np.array(bc), np.array(bf), np.array(bn)
            ax.plot([0,1],[0,1],'k--',alpha=0.4)
            ax.plot(bc, bf, 'o-', lw=2, ms=8, color='purple' if 'raw' in name else 'darkgreen')
            for c_, f_, n_ in zip(bc, bf, bn):
                ax.annotate(f'n={n_}', (c_, f_), textcoords='offset points', xytext=(5,-10), fontsize=7)
            my_ece = sum(bn*abs(bf-bc))/bn.sum()
            ax.set_title(f'{name}\nECE={my_ece:.3f}')
            ax.set_xlabel('Predicted confidence'); ax.set_ylabel('Observed TP rate')
            ax.set_xlim(0,1); ax.set_ylim(0,1); ax.grid(alpha=0.3)
        plt.tight_layout()
        os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True)
        fig.savefig(FIG_OUT, dpi=140, bbox_inches='tight')
        print(f'Saved {FIG_OUT}')

    return deliverable


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--freeze', metavar='DIR', default=None,
                    help='Also serialize the trained production artifacts to DIR '
                         '(for deployment; see deploy/masswiki_gbm/). Normal outputs '
                         'are unchanged.')
    cli = ap.parse_args()
    main(freeze_dir=cli.freeze)
